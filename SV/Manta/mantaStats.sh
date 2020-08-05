while read p; do echo "$p"; zcat "$p" |grep -v '^#' |cut -d$'\t' -f 3 |cut -d':' -f 1| sort | uniq -c; done < diploidSV.txt > mantaStats.txt
