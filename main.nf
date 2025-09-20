#!/usr/local/bin/nextflow

params.sample_sheet = null
params.outdir = null
params.top_coverage_regions = 50
params.probe_tiling_gap = 25
params.padding = 50
params.coverage_threshold = 500
params.test_dir = "${params.outdir}/test_probes"

include { TEST_PROBES } from './subworkflows/test_probes'

if (!params.sample_sheet) {
    error "Please provide a samplesheet using --sample_sheet argument"
}

if (!params.outdir) {
    error "Please provide a output dir using --outdir argument"
}

workflow {
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample_id.trim()
            def read1 = file(row.read1.trim())
            def read2 = row.containsKey('read2') && row.read2?.trim() ? file(row.read2.trim()) : null
            def is_pe = read2 ? true: false
            tuple(sample_id, read1, read2, is_pe)
        }
        .set { samples_ch }

    //MERGE_PAIRED_READS(samples_ch)
    RUN_SORTMERNA(samples_ch, "/app/idx", "${params.cpus}")
    RUN_SORTMERNA.out.map { sample_id, bam_file, read1, read2 ->
        tuple(sample_id, bam_file, read1, read2)
    }
    .set { sortmerna_out }
    GENOME_COVERAGE_BED(sortmerna_out, "/app/resources/Genome/allFasta.fasta")
    IDENTIFY_HIGH_COVERAGE_BLOCKS(GENOME_COVERAGE_BED.out, params.coverage_threshold)
    MERGE_CLOSE_BY_BLOCKS(IDENTIFY_HIGH_COVERAGE_BLOCKS.out)
    ADD_READ_PERCENT(MERGE_CLOSE_BY_BLOCKS.out)
    ADD_READ_PERCENT.out.collect().set { all_cov_sorted_bed_percent_added }
    GET_TOP_COVERAGE_REGIONS(
        params.top_coverage_regions, 
        all_cov_sorted_bed_percent_added
    )
    MERGE_OVERLAPPING_REGIONS(params.top_coverage_regions, GET_TOP_COVERAGE_REGIONS.out)
    BED_TO_FASTA(MERGE_OVERLAPPING_REGIONS.out, "/app/resources/Genome/allFasta.fasta", params.top_coverage_regions)
    RUN_BLAST(BED_TO_FASTA.out)
    IDENTIFY_HEATMAP_BLOCKS(RUN_BLAST.out.blast_hit_json, RUN_BLAST.out.cov_fasta, params.top_coverage_regions)
    MAKE_BED_25BP_GAP(IDENTIFY_HEATMAP_BLOCKS.out.top_regions_80_perc_only_fai, params.top_coverage_regions, params.probe_tiling_gap)
    GET_FASTA(MAKE_BED_25BP_GAP.out, IDENTIFY_HEATMAP_BLOCKS.out.top_regions_80_perc_only_fasta, params.top_coverage_regions)
    
    TEST_PROBES(
        samples_ch,
        GENOME_COVERAGE_BED.out,
        "/app/resources/Genome/allFasta.fasta",
        GET_FASTA.out.probes_fasta,
        GET_FASTA.out.probes_summary,
        params.top_coverage_regions
    )
}

process GET_FASTA {
    label 'medium'

    publishDir "${params.outdir}/probes/"

    input:
    path(bed_file)
    path(fasta_file)
    val(top_coverage_regions)

    output:
    path("top_${top_coverage_regions}_additional_probe_80perc_only.fasta"), emit: probes_fasta
    path("probes_summary.csv"), emit: probes_summary
    
    script:
    """
    bedtools getfasta -s -fi $fasta_file -bed $bed_file -fo top_${top_coverage_regions}_additional_probe_80perc_only.fasta
    # Create probes table
    python /app/bin/summarize_probes.py top_${top_coverage_regions}_additional_probe_80perc_only.fasta probes_summary.csv
    """
}

process MAKE_BED_25BP_GAP {
    label 'small'

    publishDir "${params.outdir}/probes/"

    input:
    path(fasta_index)
    val(top_coverage_regions)
    val(probe_tiling_gap)

    output:
    path("top_${top_coverage_regions}_additional_probe_80perc_only.bed")

    script:
    """
    /app/bin/make_bed_25bp_gap_v2.py $fasta_index "top_${top_coverage_regions}_additional_probe_80perc_only.bed" --gap $probe_tiling_gap
    """
}

process IDENTIFY_HEATMAP_BLOCKS {
    label 'small'

    publishDir "${params.outdir}/probes/"

    input:
    path(blast_hit_json)
    path(cov_fasta)
    val(top_coverage_regions)

    output:
    path("top_${top_coverage_regions}_region_to_make_probes_80perc_only.fasta"), emit: top_regions_80_perc_only_fasta
    path("top_${top_coverage_regions}_region_to_make_probes_80perc_only.fasta.fai"), emit: top_regions_80_perc_only_fai
    '*'
    
    script:
    """
    /app/bin/identify_heatmap_blocks_v2.py -t ${top_coverage_regions} -b $blast_hit_json -f $cov_fasta  -o "top_${top_coverage_regions}_region_to_make_probes_80perc_only.fasta"
    samtools faidx "top_${top_coverage_regions}_region_to_make_probes_80perc_only.fasta"
    """
}

process RUN_BLAST {
    label 'medium'

    publishDir "${params.outdir}/probes"

    input:
    path(cov_fasta)
    
    output:
    path("blastHitResult.json"), emit: blast_hit_json
    path(cov_fasta), emit: cov_fasta

    script:
    """
    makeblastdb -in $cov_fasta -dbtype 'nucl' -out db
    /app/bin/blast.py -f $cov_fasta -d db -o ./
    """
}

process BED_TO_FASTA {
    label 'medium'

    publishDir "${params.outdir}/probes/"

    input:
    path(top_cov_region_bed)
    path(all_fasta)
    val(top_coverage_regions)

    output:
    path("top${top_coverage_regions}_cov_region_seq.fasta")

    script:
    """
    fastaFromBed -fi $all_fasta -bed $top_cov_region_bed -fo top${top_coverage_regions}_cov_region_seq.fasta
    """
}

process MERGE_OVERLAPPING_REGIONS {
    label 'small'

    publishDir "${params.outdir}/probes/"

    input:
    val(top_coverage_regions)
    path(top_coverage_regions_bed)

    output:
    path("top${top_coverage_regions}_cov_region.bed")

    script:
    """
    /app/bin/merge_overlapping_regions_top.py -t ${top_coverage_regions} -i ${top_coverage_regions_bed} -o top${top_coverage_regions}_cov_region.bed
    """
}

process GET_TOP_COVERAGE_REGIONS {
    label 'small'

    publishDir "${params.outdir}/probes/"

    input:
    val(top_coverage_regions)
    val(bed_file_list)

    output:
    path("high_coverage_blocks_gap_merged_cov_sorted_percentage_added_top${top_coverage_regions}.bed")

    script:
    outfile = "high_coverage_blocks_gap_merged_cov_sorted_percentage_added_top${top_coverage_regions}.bed"
    """
    for bed_file in ${bed_file_list.join(' ')}; do
        head -n ${top_coverage_regions} \$bed_file >> $outfile
    done
    """
}

process ADD_READ_PERCENT {
    label 'small'

    publishDir "${params.outdir}/$sample_id"

    tag "$sample_id"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(cov_sorted_bed), path(read1), path(read2)

    output:
    path("${sample_id}_high_coverage_blocks_gap_merged_cov_sorted_percentage_added.bed")

    script:
    """
    read_count=`wc -l $read1 | awk '{print \$1/4;}'`
    /app/bin/add_read_percent.py -s ${sample_id} -c $cov_sorted_bed -n \$read_count
    """
}

process MERGE_CLOSE_BY_BLOCKS {
    label 'small'

    publishDir "${params.outdir}/$sample_id"

    tag "$sample_id"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(high_cov_blocks), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_cov_blocks_merged_sorted.bed"), path(read1), path(read2)

    script:
    """
    /app/bin/merge_close_by_blocks.py -s ${sample_id} -c $high_cov_blocks
    tail -n +2 ${sample_id}_high_coverage_blocks_gap_merged.tsv | sort -k 4 -nr > ${sample_id}_cov_blocks_merged_sorted.bed
    """
}

process IDENTIFY_HIGH_COVERAGE_BLOCKS {
    label 'small'

    publishDir "${params.outdir}/$sample_id"

    tag "$sample_id"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(genome_cov_bed), path(read1), path(read2)
    val(coverage_threshold)

    output:
    tuple val(sample_id), path("${sample_id}_high_coverage_blocks.tsv"), path(read1), path(read2)

    script:
    """
    /app/bin/identify_blocks.py -s ${sample_id} -c $genome_cov_bed --high ${coverage_threshold}
    """
}

process GENOME_COVERAGE_BED {
    label 'medium'

    publishDir "${params.outdir}/$sample_id"

    tag "$sample_id"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(bam_file), path(read1), path(read2)
    path(all_fasta)

    output:
    tuple val(sample_id), path("${sample_id}_genomeCoverage.bed"), path(read1), path(read2)

    script:
    """
    genomeCoverageBed -bga -ibam $bam_file -g $all_fasta > ${sample_id}_genomeCoverage.bed
    """
}

process SORT_AND_INDEX {
    label 'medium'
    
    publishDir "${params.outdir}/$sample_id"

    tag "$sample_id"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(bam_file), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_SortMeRna.sorted.bam"), path(read1), path(read2)

    script:
    """
    samtools sort $bam_file -o ${sample_id}_SortMeRna.sorted.bam
    samtools index ${sample_id}_SortMeRna.sorted.bam
    """
}

process MERGE_PAIRED_READS {
    label 'medium'

    publishDir "${params.outdir}/$sample_id"

    tag "$sample_id"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_Merged.fastq")

    script:
    def read2_arg = read2 ? "\"$read2\"" : ""

    """
    if [[ $read2_arg == "" ]]; then
        if [[ "$read1" == *.gz ]]; then
            cp "$read1" "${sample_id}_Merged.fastq.gz"
        else
            gzip -c "$read1" > "${sample_id}_Merged.fastq.gz"
        fi
    else
        reformat.sh in1="$read1" in2="$read2_arg" out="${sample_id}_Merged.fastq" threads=6
    fi

    """
}

process RUN_SORTMERNA {
    label 'high'
    tag "$sample_id"
    errorStrategy 'terminate'

    publishDir "${params.outdir}/$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2), val(is_pe)
    path(index_files)
    val(cpus)

    output:
    tuple val(sample_id), path("${sample_id}_SortMeRna.sorted.bam"), path(read1), path(read2)

    script:
    def ref_base = "${projectDir}/resources/rRNA_databases"
    def read2_opt = read2 ? "--reads ${read2}" : ""
    """
    sortmerna \
      --workdir './' \
      --ref ${ref_base}/silva-arc-23s-id98.fasta \
      --ref ${ref_base}/silva-bac-23s-id98.fasta \
      --ref ${ref_base}/silva-bac-16s-id90.fasta \
      --ref ${ref_base}/rfam-5.8s-database-id98.fasta \
      --ref ${ref_base}/silva-euk-18s-id95.fasta \
      --ref ${ref_base}/rfam-5s-database-id98.fasta \
      --ref ${ref_base}/silva-arc-16s-id95.fasta \
      --ref ${ref_base}/silva-euk-28s-id98.fasta \
      --reads ${read1} \
      ${read2_opt} \
      --aligned ${sample_id}_SortMeRna \
      --sam \
      --threads ${cpus} \
      --SQ \
      --num_alignments 0
    
    sam_file="${sample_id}_SortMeRna.sam"
    [[ ! -f "\$sam_file" ]] && sam_file="${sample_id}_SortMeRna.sam.gz"

    samtools view -@ ${cpus} -b "\$sam_file" | \
        samtools sort -@ ${cpus} -o "${sample_id}_SortMeRna.sorted.bam"
    samtools index -@ ${cpus} "${sample_id}_SortMeRna.sorted.bam"

    rm -rf "\$sam_file"
    """
}
