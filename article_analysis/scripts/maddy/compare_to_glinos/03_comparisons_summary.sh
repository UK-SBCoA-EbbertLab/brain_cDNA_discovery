#!/bin/bash

o_t_gtf="extended_annotation_grch38_only_novel.sorted.just_transcripts.gtf"
o_e_gtf="extended_annotation_grch38_only_novel.sorted.just_exons.gtf"

g_t_gtf="flair_filter_transcripts_no_chr.sorted.just_transcripts.gtf"
g_e_gtf="flair_filter_transcripts_no_chr.sorted.just_exons.gtf"

o_v_e_t_intersect="ours_vs_ensembl/ours_vs_ensembl_just_transcripts.intersect.gtf"
o_v_e_5_t_intersect="ours_vs_ensembl/ours_vs_ensembl_just_transcripts_0_5.intersect.gtf"
o_v_e_9_t_intersect="ours_vs_ensembl/ours_vs_ensembl_just_transcripts_0_9.intersect.gtf"

o_v_e_t_no_overlap="ours_vs_ensembl/ours_vs_ensembl_just_transcripts.no_overlap.gtf"
o_v_e_5_t_no_overlap="ours_vs_ensembl/ours_vs_ensembl_just_transcripts_0_5.no_overlap.gtf"
o_v_e_9_t_no_overlap="ours_vs_ensembl/ours_vs_ensembl_just_transcripts_0_9.no_overlap.gtf"

o_v_e_e_intersect="ours_vs_ensembl/ours_vs_ensembl_just_exons.intersect.gtf"
o_v_e_5_e_intersect="ours_vs_ensembl/ours_vs_ensembl_just_exons_0_5.intersect.gtf"
o_v_e_9_e_intersect="ours_vs_ensembl/ours_vs_ensembl_just_exons_0_9.intersect.gtf"

o_v_e_e_no_overlap="ours_vs_ensembl/ours_vs_ensembl_just_exons.no_overlap.gtf"
o_v_e_5_e_no_overlap="ours_vs_ensembl/ours_vs_ensembl_just_exons_0_5.no_overlap.gtf"
o_v_e_9_e_no_overlap="ours_vs_ensembl/ours_vs_ensembl_just_exons_0_9.no_overlap.gtf"

g_v_e_t_intersect="glinos_vs_ensembl/glinos_vs_ensembl_just_transcripts.intersect.gtf"
g_v_e_5_t_intersect="glinos_vs_ensembl/glinos_vs_ensembl_just_transcripts_0_5.intersect.gtf"
g_v_e_9_t_intersect="glinos_vs_ensembl/glinos_vs_ensembl_just_transcripts_0_9.intersect.gtf"

g_v_e_t_no_overlap="glinos_vs_ensembl/glinos_vs_ensembl_just_transcripts.no_overlap.gtf"
g_v_e_5_t_no_overlap="glinos_vs_ensembl/glinos_vs_ensembl_just_transcripts_0_5.no_overlap.gtf"
g_v_e_9_t_no_overlap="glinos_vs_ensembl/glinos_vs_ensembl_just_transcripts_0_9.no_overlap.gtf"

g_v_e_e_intersect="glinos_vs_ensembl/glinos_vs_ensembl_just_exons.intersect.gtf"
g_v_e_5_e_intersect="glinos_vs_ensembl/glinos_vs_ensembl_just_exons_0_5.intersect.gtf"
g_v_e_9_e_intersect="glinos_vs_ensembl/glinos_vs_ensembl_just_exons_0_9.intersect.gtf"

g_v_e_e_no_overlap="glinos_vs_ensembl/glinos_vs_ensembl_just_exons.no_overlap.gtf"
g_v_e_5_e_no_overlap="glinos_vs_ensembl/glinos_vs_ensembl_just_exons_0_5.no_overlap.gtf"
g_v_e_9_e_no_overlap="glinos_vs_ensembl/glinos_vs_ensembl_just_exons_0_9.no_overlap.gtf"

o_v_g_t_intersect="ours_vs_glinos/ours_vs_glinos_just_transcripts.intersect.gtf"
o_v_g_5_t_intersect="ours_vs_glinos/ours_vs_glinos_just_transcripts_0_5.intersect.gtf"
o_v_g_9_t_intersect="ours_vs_glinos/ours_vs_glinos_just_transcripts_0_9.intersect.gtf"

o_v_g_t_no_overlap="ours_vs_glinos/ours_vs_glinos_just_transcripts.no_overlap.gtf"
o_v_g_5_t_no_overlap="ours_vs_glinos/ours_vs_glinos_just_transcripts_0_5.no_overlap.gtf"
o_v_g_9_t_no_overlap="ours_vs_glinos/ours_vs_glinos_just_transcripts_0_9.no_overlap.gtf"

o_v_g_e_intersect="ours_vs_glinos/ours_vs_glinos_just_exons.intersect.gtf"
o_v_g_5_e_intersect="ours_vs_glinos/ours_vs_glinos_just_exons_0_5.intersect.gtf"
o_v_g_9_e_intersect="ours_vs_glinos/ours_vs_glinos_just_exons_0_9.intersect.gtf"

o_v_g_e_no_overlap="ours_vs_glinos/ours_vs_glinos_just_exons.no_overlap.gtf"
o_v_g_5_e_no_overlap="ours_vs_glinos/ours_vs_glinos_just_exons_0_5.no_overlap.gtf"
o_v_g_9_e_no_overlap="ours_vs_glinos/ours_vs_glinos_just_exons_0_9.no_overlap.gtf"

g_v_o_t_intersect="glinos_vs_ours/glinos_vs_ours_just_transcripts.intersect.gtf"
g_v_o_5_t_intersect="glinos_vs_ours/glinos_vs_ours_just_transcripts_0_5.intersect.gtf"
g_v_o_9_t_intersect="glinos_vs_ours/glinos_vs_ours_just_transcripts_0_9.intersect.gtf"

g_v_o_t_no_overlap="glinos_vs_ours/glinos_vs_ours_just_transcripts.no_overlap.gtf"
g_v_o_5_t_no_overlap="glinos_vs_ours/glinos_vs_ours_just_transcripts_0_5.no_overlap.gtf"
g_v_o_9_t_no_overlap="glinos_vs_ours/glinos_vs_ours_just_transcripts_0_9.no_overlap.gtf"

g_v_o_e_intersect="glinos_vs_ours/glinos_vs_ours_just_exons.intersect.gtf"
g_v_o_5_e_intersect="glinos_vs_ours/glinos_vs_ours_just_exons_0_5.intersect.gtf"
g_v_o_9_e_intersect="glinos_vs_ours/glinos_vs_ours_just_exons_0_9.intersect.gtf"

g_v_o_e_no_overlap="glinos_vs_ours/glinos_vs_ours_just_exons.no_overlap.gtf"
g_v_o_5_e_no_overlap="glinos_vs_ours/glinos_vs_ours_just_exons_0_5.no_overlap.gtf"
g_v_o_9_e_no_overlap="glinos_vs_ours/glinos_vs_ours_just_exons_0_9.no_overlap.gtf"

o_t_gtf_c=$(wc -l < $o_t_gtf)
o_e_gtf_c=$(wc -l < $o_e_gtf)

g_t_gtf_c=$(wc -l < $g_t_gtf)
g_e_gtf_c=$(wc -l < $g_e_gtf)

ove_t_i=$(wc -l < "$o_v_e_t_intersect")
ove5_t_i=$(wc -l < "$o_v_e_5_t_intersect")
ove9_t_i=$(wc -l < "$o_v_e_9_t_intersect")

ove_t_no=$(wc -l < "$o_v_e_t_no_overlap")
ove5_t_no=$(wc -l < "$o_v_e_5_t_no_overlap")
ove9_t_no=$(wc -l < "$o_v_e_9_t_no_overlap")

ove_e_i=$(wc -l < "$o_v_e_e_intersect")
ove5_e_i=$(wc -l < "$o_v_e_5_e_intersect")
ove9_e_i=$(wc -l < "$o_v_e_9_e_intersect")

ove_e_no=$(wc -l < $o_v_e_e_no_overlap)
ove5_e_no=$(wc -l < $o_v_e_5_e_no_overlap)
ove9_e_no=$(wc -l < $o_v_e_9_e_no_overlap)


gve_t_i=$(wc -l < $g_v_e_t_intersect)
gve5_t_i=$(wc -l < $g_v_e_5_t_intersect)
gve9_t_i=$(wc -l < $g_v_e_9_t_intersect)

gve_t_no=$(wc -l < $g_v_e_t_no_overlap)
gve5_t_no=$(wc -l < $g_v_e_5_t_no_overlap)
gve9_t_no=$(wc -l < $g_v_e_9_t_no_overlap)

gve_e_i=$(wc -l < $g_v_e_e_intersect)
gve5_e_i=$(wc -l < $g_v_e_5_e_intersect)
gve9_e_i=$(wc -l < $g_v_e_9_e_intersect)

gve_e_no=$(wc -l < $g_v_e_e_no_overlap)
gve5_e_no=$(wc -l < $g_v_e_5_e_no_overlap)
gve9_e_no=$(wc -l < $g_v_e_9_e_no_overlap)


ovg_t_i=$(wc -l < $o_v_g_t_intersect)
ovg5_t_i=$(wc -l < $o_v_g_5_t_intersect)
ovg9_t_i=$(wc -l < $o_v_g_9_t_intersect)

gvo_t_i=$(wc -l < $g_v_o_t_intersect)
gvo5_t_i=$(wc -l < $g_v_o_5_t_intersect)
gvo9_t_i=$(wc -l < $g_v_o_9_t_intersect)

ovg_t_no=$(wc -l < $o_v_g_t_no_overlap)
ovg5_t_no=$(wc -l < $o_v_g_5_t_no_overlap)
ovg9_t_no=$(wc -l < $o_v_g_9_t_no_overlap)

gvo_t_no=$(wc -l < $g_v_o_t_no_overlap)
gvo5_t_no=$(wc -l < $g_v_o_5_t_no_overlap)
gvo9_t_no=$(wc -l < $g_v_o_9_t_no_overlap)


ovg_e_i=$(wc -l < $o_v_g_e_intersect)
ovg5_e_i=$(wc -l < $o_v_g_5_e_intersect)
ovg9_e_i=$(wc -l < $o_v_g_9_e_intersect)

gvo_e_i=$(wc -l < $g_v_o_e_intersect)
gvo5_e_i=$(wc -l < $g_v_o_5_e_intersect)
gvo9_e_i=$(wc -l < $g_v_o_9_e_intersect)

ovg_e_no=$(wc -l < $o_v_g_e_no_overlap)
ovg5_e_no=$(wc -l < $o_v_g_5_e_no_overlap)
ovg9_e_no=$(wc -l < $o_v_g_9_e_no_overlap)

gvo_e_no=$(wc -l < $g_v_o_e_no_overlap)
gvo5_e_no=$(wc -l < $g_v_o_5_e_no_overlap)
gvo9_e_no=$(wc -l < $g_v_o_9_e_no_overlap)

i_com_o="intersect_comparison_summary.tsv"
no_com_o="no_overlap_comparison_summary.tsv"

echo -e "comparisons\tintersect\tintersect at least 50% overlap\tintersect at least 90% overlap" > ${i_com_o}
echo -e "ours v ensembl transcripts\t${ove_t_i}/${o_t_gtf_c}\t${ove5_t_i}/${o_t_gtf_c}\t${ove9_t_i}/${o_t_gtf_c}" >> ${i_com_o}
echo -e "ours v ensembl exons\t${ove_e_i}/${o_e_gtf_c}\t${ove5_e_i}/${o_e_gtf_c}\t${ove9_e_i}/${o_e_gtf_c}" >> ${i_com_o}
echo -e "glinos v ensembl transcripts\t${gve_t_i}/${g_t_gtf_c}\t${gve5_t_i}/${g_t_gtf_c}\t${gve9_t_i}/${g_t_gtf_c}" >> ${i_com_o}
echo -e "glinos v ensembl exons\t${gve_e_i}/${g_e_gtf_c}\t${gve5_e_i}/${g_e_gtf_c}\t${gve9_e_i}/${g_e_gtf_c}" >> ${i_com_o}
#echo -e "ours v glinos transcripts\t${ovg_t_i}\t${ovg5_t_i}\t${ovg9_t_i}" >> ${i_com_o}
#echo -e "ours v glinos exons\t${ovg_e_i}\t${ovg5_e_i}\t${ovg9_e_i}" >> ${i_com_o}
#echo -e "glinos v ours transcripts\t${gvo_t_i}\t${gvo5_t_i}\t${gvo9_t_i}" >> ${i_com_o}
#echo -e "glinos v ours exons\t${gvo_e_i}\t${gvo5_e_i}\t${gvo9_e_i}" >> ${i_com_o}

echo -e "comparisons\tno overlap\tno overlap (intersect at least 50% overlap)\tno overlap (intersect at least 90% overlap)" > ${no_com_o}
echo -e "ours v ensembl transcripts\t${ove_t_no}/${o_t_gtf_c}\t${ove5_t_no}/${o_t_gtf_c}\t${ove9_t_no}"/${o_t_gtf_c} >> ${no_com_o}
echo -e "ours v ensembl exons\t${ove_e_no}/${o_e_gtf_c}\t${ove5_e_no}/${o_e_gtf_c}\t${ove9_e_no}/${o_e_gtf_c}" >> ${no_com_o}
echo -e "glinos v ensembl transcripts\t${gve_t_no}/${g_t_gtf_c}\t${gve5_t_no}/${g_t_gtf_c}\t${gve9_t_no}/${g_t_gtf_c}" >> ${no_com_o}
echo -e "glinos v ensembl exons\t${gve_e_no}/${g_e_gtf_c}\t${gve5_e_no}/${g_e_gtf_c}\t${gve9_e_no}/${g_e_gtf_c}" >> ${no_com_o}
#echo -e "ours v glinos transcripts\t${ovg_t_no}\t${ovg5_t_no}\t${ovg9_t_no}" >> ${no_com_o}
#echo -e "ours v glinos exons\t${ovg_e_no}\t${ovg5_e_no}\t${ovg9_e_no}" >> ${no_com_o}
#echo -e "glinos v ours transcripts\t${gvo_t_no}\t${gvo5_t_no}\t${gvo9_t_no}" >> ${no_com_o}
#echo -e "glinos v ours exons\t${gvo_e_no}\t${gvo5_e_no}\t${gvo9_e_no}" >> ${no_com_o}
