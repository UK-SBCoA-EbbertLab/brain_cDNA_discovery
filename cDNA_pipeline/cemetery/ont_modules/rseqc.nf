process RSEQC {

    publishDir 'results/RseQC/', mode: 'copy', overwrite: false
 
    label "medium_small"

    input:
    val(id)
    path(bam)
    path(bai)
    path(housekeep)
    path(bed)

    output:
    path "*"

    script:
    """
    bam_stat.py -i $bam > "${id}_bam_stat"
    #geneBody_coverage.py -r $housekeep -i $bam  -o "${id}_gene_body_coverage"
    junction_annotation.py -i $bam -o "junction_annotation" -r $bed 2> "${id}_junction_annotation_multiqc"
    junction_saturation.py -i $bam -r $bed -o "${id}_junction_saturation"
    read_GC.py -i $bam -o "${id}_read_gc"
    read_NVC.py -i $bam -o "${id}_read_nvc"
    read_quality.py -i $bam -o "${id}_read_quality"
    read_distribution.py  -i $bam -r $bed > "${id}_read_distribution"
    read_duplication.py -i $bam -o "${id}_read_duplication"
    """
}
