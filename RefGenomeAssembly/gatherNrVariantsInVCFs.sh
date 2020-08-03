#find /research/btc_bioinformatic/operations/scratch/ -type f -name "*_HaplotypeCaller.vcf.gz*" > hg19ref_haplotypeCaller.txt
declare -a samples
dirs=( scratch )
for wdir in "${dirs[@]}"; do
    samples=()
    for i in $(find /research/btc_bioinformatic/operations/${wdir}/ -type f -name "*_HaplotypeCaller.vcf.gz"); do
	sample=$(basename $i|cut -d'_' -f1)
	statsfile="variantStats_${sample}_${wdir}.csv"
	echo "${sample}${wdir:0:1} chromo" > $statsfile #"$(basename $i)"
	samples+=($sample)    
	zcat $i|grep -v '^#' | cut -d$'\t' -f 1|uniq -c > $statsfile
    done
done
