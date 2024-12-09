# Prepend a value to every line in a list
for i in $(seq 1 8); do
	echo $i;
done | xargs -L 1 echo prepend


# Can also run a function
function times2 {
	echo "$1 * 2 =" $(( $1 * 2 ))
	# Simulates processing time
	sleep 1
}
# Export function
export -f times2
time for i in $(seq 1 8); do
	echo $i;
# -I takes a whole line (-L 1 implied) and substitutes it as the matching symbol
done | xargs -I {} bash -c 'times2 {}'


time for i in $(seq 1 8); do
	echo $i;
# -I takes a whole line (-L 1 implied) and substitutes it as the matching symbol
# Can also run in parallel
done | xargs -P 4 -I {} bash -c 'times2 {}'
