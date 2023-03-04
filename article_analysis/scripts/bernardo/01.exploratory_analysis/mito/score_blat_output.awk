#!/usr/bin/awk -f
BEGIN {
	OFS="\t"
        {print "transcript_id","t_align_start","t_align_end","chr","strand","chr_start","chr_end","percent_identity"}
}
{
	qSize = $13 - $12
	tSize = $17 - $16
	totalSize = qSize > tSize ? qSize : tSize
	score = $1 / totalSize
	if (score > .90) {
		print $14,$16,$17,$10,$9,$12,$13,score
	}
}
