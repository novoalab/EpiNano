#!/usr/bin/env bash

usage="sh $(basename "$0") [-m -u -o]"


if [[ "$#" -eq 0 ]]; then
	echo "$usage"
	exit 0
fi 


while getopts ":hm:u:o:" option; do 
	case $option in 
		h) 
			echo -e " 
		Usage: 
			$usage \n
			-m modified sample feature table [either slided or not]
			-u unmodified sample feature table [either slided or not]
			-o output file name 
			"
			exit;;
		m) mod=${OPTARG} ;;
		u) unm=${OPTARG} ;;
		o) out=${OPTARG};;
		
		:) printf "missing argment for -%s\n" "$OPTARG" >&2
			printf "$usage" >&2
			exit 1 
		;;
		
		\?) # incrorrect option
		echo "Error: invalid option: -%s\n" "$OPTARG" >&2
		exit 1 ;;
	esac 
done
	

awk '{if (/^#/) print $0",sample"; else print $0",mod";}' ${mod}  > ${mod}.labeled

awk '{if (/^#/) print $0",sample"; else print $0",unm";}' ${unm} > ${unm}.labeled 

cat ${mod}.labeled ${unm}.labeled | awk '!(NR>1 && /^#/)' > ${out} && rm ${mod}.labeled ${unm}.labeled

