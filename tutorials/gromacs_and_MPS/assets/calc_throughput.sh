grep -H "Performance" $@ | sed -e "s/\s\+/\t/g" | cut -f 1-2
grep -H "Performance" $@ | sed -e "s/\s\+/\t/g" | cut -f 2 | \
	awk '{sum+=$1} END {print "Total throughput:",sum}'
