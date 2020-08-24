#!/usr/bin/env bash

usage="$(basename "$0") [-h] [-e s] [-n n] -- program to sort collapsed events table.
	where:
		-h show this help text
		-e collapsed current intensity table
		-n number of cpus used for sorting the input (default: 6)"
num_cpus=6

while getopts ':he:n:' option; do 
	case "$option" in
	h) echo "$usage"
		exit
		;;
	e) events=$OPTARG
		;;
	n) num_cpus=$OPTARG
		;;
	:) printf "missing argument for -%e\n" "$OPTARG" >&2
		echo "$usage" >&2
		exit 1 
		;;
	:) printf "missing argument for -%n\n" "$OPTARG" >&2
		echo "$usage" >&2
		exit 1 
		;;
	esac 
done 
	


events=$1
num_cpus=$2

out=${events}.sort

if [ ${events: -3} == '.gz' ] ; then
	pigz -cd ${events} | sort -S 50% --parallel=8 > $out
    else
	sort -S 50% --parallel=8 ${events} > $out
fi




