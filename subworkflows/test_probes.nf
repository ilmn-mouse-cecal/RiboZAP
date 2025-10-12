#!/usr/local/bin/nextflow

workflow TEST_PROBES {
    take:
        samples_ch
        merged_fastq
        ref_fasta
        additional_probe_80_percent_fasta
        probes_summary
        top_coverage_regions

    main:
        RUN_SORTMERNA_BEST_HIT(merged_fastq, "/app/idx", "${params.cpus}")
        RUN_BLAST(ref_fasta, additional_probe_80_percent_fasta, top_coverage_regions)
        FILTER_AND_ADD_PADDING(RUN_BLAST.out, ref_fasta, top_coverage_regions, params.padding)
        MERGE_CAN_DEPLETE_REGIONS(FILTER_AND_ADD_PADDING.out, top_coverage_regions)
        GET_NEAR_PROBE_READS(
            params.analysis_name.replaceAll("_", "-"),
            RUN_SORTMERNA_BEST_HIT.out,
            MERGE_CAN_DEPLETE_REGIONS.out.can_deplete_regions_merged,
            top_coverage_regions
        )

        samples_ch.collect(flat: false).set {all_samples}
        RUN_SORTMERNA_BEST_HIT.out.collect(flat: false).set {sortmerna_bam}
        GET_NEAR_PROBE_READS.out.near_probe_bam.collect(flat: false).set {near_probe_reads}

        CALCULATE_STATS(
            all_samples,
            sortmerna_bam,
            near_probe_reads
        )
        new File("NO_SUMMARY").text = ""
        probes_summary = probes_summary ? probes_summary: Channel.fromPath("NO_SUMMARY")
        GENERATE_REPORTS(
            CALCULATE_STATS.out,
            additional_probe_80_percent_fasta,
            probes_summary
        )
}

process GENERATE_REPORTS {
    label 'small'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path "top_coverage_result.csv"
    path probes_fasta
    path probes_summary

    output:
    path "./reports"

    script:
    """
    cp /app/bin/probe_download.md .
    cp top_coverage_result.csv test_probes_rrna_reduction.csv
    cp top_coverage_result.csv test_probes_composition.csv
    cut -f1,7,8,9 -d, top_coverage_result.csv >test_probes_heatmap.csv
    rm -rf reports &&  multiqc . -o reports --config /app/config/multiqc_custom.yaml
    cp $probes_fasta reports/probes.fasta
    """
}

process CALCULATE_STATS {
    label 'small'

    publishDir "${params.test_dir}"
    errorStrategy 'finish'

    input:
    val all_samples
    val sortmerna_bam_files
    val near_probe_bam_files

    output:
    val "${params.test_dir}/top_coverage_result.csv"

    exec:
    def samples_map = all_samples.collectEntries { [it[0], it[1]] }
    def bam_map = sortmerna_bam_files.collectEntries { [it[0], it[1]] }
    def near_probe_bam_map = near_probe_bam_files.collectEntries { [it[0], it[1]] }
    def is_pe = all_samples[0][3]

    def resultFile = new File(params.test_dir.toString() + "/top_coverage_result.csv")
    resultFile.text = "Sample ID,Total Reads,Total Mapped,Unmapped,Depleted,Remaining Mapped,rRNA Mapped Percent,Remaining rRNA Percent,rRNA Depletion Percent\n"

    samples_map.keySet().each { sample_id ->
        def fastq_path = samples_map[sample_id]
        def bam_path = bam_map[sample_id]
        def near_probe_bam_path = near_probe_bam_map[sample_id]

        def depleted = "samtools view -c ${near_probe_bam_path}".execute().text.trim().toInteger()
        def mappedReads = ["bash", "-c", "samtools view -F 4 ${bam_path}| wc -l"].execute().text.trim().toInteger()
        def cmd = fastq_path.toString().endsWith(".gz") ?
            "zcat ${fastq_path} | wc -l" :
            "cat ${fastq_path} | wc -l"
        def totalfastq_lines = ["bash", "-c", cmd].execute().text.trim().toInteger()
        def totalfastq = is_pe ? (totalfastq_lines / 4) *2 : (totalfastq_lines / 4)
        def unmapped = totalfastq - mappedReads
        def mapped_percent = (mappedReads / totalfastq) * 100
        def remaining_rRNA = mappedReads - depleted
        def remaining_rRNA_percent = (remaining_rRNA / mappedReads) * 100
        def rrna_depletion_percent = (depleted / mappedReads) * 100
        resultFile.append("${sample_id},${totalfastq},${mappedReads},${unmapped},${depleted},${remaining_rRNA},${String.format('%.2f', mapped_percent)},${String.format('%.2f', remaining_rRNA_percent)},${String.format('%.2f', rrna_depletion_percent)}\n")
    }
}

/*process CALCULATE_STATS {
    label 'small'
    tag "$sample_id"
    publishDir "${params.test_dir}"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(read1), path(read2)
    tuple val(sample_id), path(sorted_bam)
    path(near_probe_reads_sam)
    path(top_coverage_result)

    output:
    path(top_coverage_result)

    script:
    """
    depleted=`cat $near_probe_reads_sam  | wc -l`
    totalmapped=`samtools view -F 4 $sorted_bam | wc -l`
    totalfastq=`wc -l $read1 | awk '{print \$1/4*2;}'`
    echo ${sample_id}","\${totalmapped}","\${depleted}","\${totalfastq} | awk 'BEGIN { FS=",";OFS = "\t";}{print \$1,\$2,\$3,\$4,(\$2/\$4*100)"%",((\$2-\$3)/(\$4-\$3)*100)"%",((\$2/\$4*100)-((\$2-\$3)/(\$4-\$3)*100))"%";}' >> $top_coverage_result
    """
}*/

process GET_NEAR_PROBE_READS {
    label 'medium'

    tag "$sample_id"

    publishDir "${params.test_dir}/$sample_id"
    errorStrategy 'ignore'

    input:
    val analysis_name
    tuple val(sample_id), path(sorted_bam)
    path(can_deplete_regions_bed)
    val(top_coverage_regions)

    output:
    tuple val(sample_id), path("top_${top_coverage_regions}_additional_probe_80perc_only_near_probe_reads.bam"), emit: near_probe_bam
    tuple val(sample_id), path("top_${top_coverage_regions}_additional_probe_80perc_only_NOT_near_probe_reads.bam"), emit: not_near_probe_bam
    tuple val(sample_id), path("${sample_id}_${analysis_name}-residual-rRNA_S1_L001_R1_001.fastq.gz"), emit: rRNA_fastq

    script:
    
    """
    samtools view -b $sorted_bam -L $can_deplete_regions_bed > top_${top_coverage_regions}_additional_probe_80perc_only_near_probe_reads.bam
    samtools view -b -L $can_deplete_regions_bed -U top_${top_coverage_regions}_additional_probe_80perc_only_NOT_near_probe_reads.bam $sorted_bam > /dev/null
    samtools fastq -0 /dev/stdout -s /dev/null -n top_${top_coverage_regions}_additional_probe_80perc_only_NOT_near_probe_reads.bam | gzip > ${sample_id}-${analysis_name}-residual-rRNA_S1_L001_R1_001.fastq.gz
    """
}

process MERGE_CLOSE_BY_BLOCKS {
    label 'small'

    publishDir "${params.test_dir}/$sample_id"

    tag "$sample_id"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(high_cov_blocks)

    output:
    tuple val(sample_id), path("${sample_id}_cov_blocks_merged_sorted.bed")

    script:
    """
    /app/bin/merge_close_by_blocks.py -s ${sample_id} -c $high_cov_blocks
    tail -n +2 ${sample_id}_high_coverage_blocks_gap_merged.tsv | sort -k 4 -nr > ${sample_id}_cov_blocks_merged_sorted.bed
    """
}

process GENOME_COVERAGE_BED {
    label 'medium'

    publishDir "${params.test_dir}/$sample_id"

    tag "$sample_id"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(bam_file)
    path(ref_fasta)

    output:
    tuple val(sample_id), path("${sample_id}_genomeCoverage.bed")

    script:
    """
    genomeCoverageBed -bga -ibam $bam_file -g $ref_fasta > ${sample_id}_genomeCoverage.bed
    """
}

process IDENTIFY_ALL_COVERAGE_BLOCKS {
    label 'small'

    publishDir "${params.test_dir}/$sample_id"

    tag "$sample_id"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(genome_cov_bed)

    output:
    tuple val(sample_id), path("${sample_id}_all_coverage_blocks.tsv")

    script:
    """
    /app/bin/identify_blocks.py -s ${sample_id} -c $genome_cov_bed --high 1
    mv "${sample_id}_high_coverage_blocks.tsv" "${sample_id}_all_coverage_blocks.tsv"
    """
}

process RUN_SORTMERNA_BEST_HIT {
    label 'high'

    tag "$sample_id"

    publishDir "${params.test_dir}/$sample_id"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(merged_fastq)
    path(index_files)
    val(cpus)

    output:
    tuple val(sample_id), path("${sample_id}_SortMeRna.sorted.bam")

    script:
    def ref_base = "/app/resources/rRNA_databases"
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
      --reads ${merged_fastq} \
      --aligned ${sample_id}_SortMeRna \
      --threads ${cpus} \
      --sam \
      --SQ \
      --num_alignments 1

    samtools view -Sb "${sample_id}_SortMeRna.sam" > "${sample_id}_SortMeRna.bam"
    samtools sort "${sample_id}_SortMeRna.bam" -o "${sample_id}_SortMeRna.sorted.bam"
    rm -rf "${sample_id}_SortMeRna.sam"
    rm -rf "${sample_id}_SortMeRna.bam"
    """
}

process MERGE_CAN_DEPLETE_REGIONS {
    label 'small'

    publishDir "${params.test_dir}"

    input:
    path(can_deplete_regions)
    val(top_coverage_regions)

    output:
    path("top_${top_coverage_regions}_additional_probe_80perc_only_can_deplete_regions_merged.bed"), emit: can_deplete_regions_merged
    
    script:
    """
    /app/bin/merge_candeplete_regions.py -i $can_deplete_regions -o top_${top_coverage_regions}_additional_probe_80perc_only_can_deplete_regions_merged.bed
    """
}

process FILTER_AND_ADD_PADDING {
    label 'medium'

    publishDir "${params.test_dir}"

    input:
    path(blast_result_txt)
    path(ref_fasta)
    val(top_coverage_regions)
    val(padding)

    output:
    path("top_${top_coverage_regions}_additional_probe_80perc_only_can_deplete_regions_sorted.txt")
    
    script:
    """
    /app/bin/filter_add_padding.py -i $blast_result_txt -f $ref_fasta -o top_${top_coverage_regions}_additional_probe_80perc_only_can_deplete_regions.txt -p ${padding}
    sortBed -i top_${top_coverage_regions}_additional_probe_80perc_only_can_deplete_regions.txt > top_${top_coverage_regions}_additional_probe_80perc_only_can_deplete_regions_sorted.txt
    """
}

process RUN_BLAST {
    label 'medium'

    publishDir "${params.test_dir}"

    input:
    path(ref_fasta)
    path(additional_probe_80_percent_fasta)
    val(top_coverage_regions)

    output:
    path("top_${top_coverage_regions}_additional_probe_80perc_only_blast_result.txt")

    script:
    """
    makeblastdb -dbtype nucl -in $ref_fasta -out db
    blastn -db db -query $additional_probe_80_percent_fasta -evalue 10 -outfmt 6 -out top_${top_coverage_regions}_additional_probe_80perc_only_blast_result.txt
    """
}