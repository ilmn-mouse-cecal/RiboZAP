#!/usr/local/bin/nextflow

workflow TEST_PROBES {
    take:
        merged_reads
        genome_cov_bed
        ref_fasta
        additional_probe_80_percent_fasta
        probes_summary
        top_coverage_regions

    main:
        RUN_SORTMERNA_BEST_HIT(merged_reads, "/app/idx", "${params.cpus}")
        RUN_BLAST(ref_fasta, additional_probe_80_percent_fasta, top_coverage_regions)
        FILTER_AND_ADD_PADDING(RUN_BLAST.out, ref_fasta, top_coverage_regions, params.padding)
        MERGE_CAN_DEPLETE_REGIONS(FILTER_AND_ADD_PADDING.out, top_coverage_regions)
        //GENOME_COVERAGE_BED(RUN_SORTMERNA_BEST_HIT.out, ref_fasta)
        //IDENTIFY_ALL_COVERAGE_BLOCKS(GENOME_COVERAGE_BED.out)
        //MERGE_CLOSE_BY_BLOCKS(IDENTIFY_ALL_COVERAGE_BLOCKS.out)
        GET_NEAR_PROBE_READS(
            RUN_SORTMERNA_BEST_HIT.out,
            MERGE_CAN_DEPLETE_REGIONS.out.can_deplete_regions_merged,
            top_coverage_regions
        )

        merged_reads.collect(flat: false).set {all_samples}
        RUN_SORTMERNA_BEST_HIT.out.collect(flat: false).set {srotmerna_bam}
        GET_NEAR_PROBE_READS.out.collect(flat: false).set {near_probe_reads}

        CALCULATE_STATS(
            all_samples,
            srotmerna_bam,
            near_probe_reads
        )
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
    val near_probe_sam_files

    output:
    val "${params.test_dir}/top_coverage_result.csv"

    exec:
    def fastq_map = all_samples.collectEntries { [it[0], it[1]] }
    def bam_map = sortmerna_bam_files.collectEntries { [it[0], it[1]] }
    def sam_map = near_probe_sam_files.collectEntries { [it[0], it[1]] }

    def resultFile = new File(params.test_dir.toString() + "/top_coverage_result.csv")
    resultFile.text = "Sample ID,Total Reads,Total Mapped,Remaining Mapped,Unmapped,Depleted,Mapped Percent,Depleted Mapped Percent,rRNA Depletion Percent\n"

    fastq_map.keySet().each { sample_id ->
        def merged_fastq_path = fastq_map[sample_id]
        def bam_path = bam_map[sample_id]
        def sam_path = sam_map[sample_id]

        def depleted = "grep -vc '^@' ${sam_path}".execute().text.trim().toInteger()
        def totalmapped = "samtools view -c -F 4 ${bam_path}".execute().text.trim().toInteger()
        def totalfastq_lines = new File(merged_fastq_path.toString()).readLines().size()
        def totalfastq = totalfastq_lines / 4

        def mapped_percent = (totalmapped / totalfastq) * 100
        def depleted_mapped_percent = ((totalmapped - depleted) / (totalfastq - depleted)) * 100
        def diff_percent = mapped_percent - depleted_mapped_percent
        def remaining_mapped = totalmapped - depleted
        def unmapped = totalfastq - totalmapped
        resultFile.append("${sample_id},${totalfastq},${totalmapped},${remaining_mapped},${unmapped},${depleted},${String.format('%.2f', mapped_percent)},${String.format('%.2f', depleted_mapped_percent)},${String.format('%.2f', diff_percent)}\n")
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
    tuple val(sample_id), path(sorted_bam)
    path(can_deplete_regions_bed)
    val(top_coverage_regions)

    output:
    tuple val(sample_id), path("top_${top_coverage_regions}_additional_probe_80perc_only_near_probe_reads.sam")

    script:
    """
    samtools view $sorted_bam -L $can_deplete_regions_bed > top_${top_coverage_regions}_additional_probe_80perc_only_near_probe_reads.sam
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

    samtools view -b "${sample_id}_SortMeRna.sam" | \
    samtools sort -o "${sample_id}_SortMeRna.sorted.bam"
    samtools index "${sample_id}_SortMeRna.sorted.bam"
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
    path("top_${top_coverage_regions}_result.txt"), emit: top_coverage_result
    
    script:
    """
    /app/bin/merge_candeplete_regions.py -i $can_deplete_regions -o top_${top_coverage_regions}_additional_probe_80perc_only_can_deplete_regions_merged.bed
    echo -e "SampleNumber\tTotal # of Reads Mapped to SIVLA\tReads Overlaps with Possible Depleted Region\tTotal Number of Reads (R1+R2)\tOriginal rRNA Contents\tEstimate rRNA Contents After Extra Probes\trRNA Reduction" > top_${top_coverage_regions}_result.txt
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