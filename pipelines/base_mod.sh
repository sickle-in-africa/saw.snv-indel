
custom_call() {
	local argv=("$@") \
		msg \
		routine

	msg=${argv[0]}
	routine=${argv[1]}

	echo $msg
	$routine || { echo "...failed"; return 1; }
}