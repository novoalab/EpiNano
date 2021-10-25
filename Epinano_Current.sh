#!/usr/bin/env bash 

usage="$(basename "$0") [-h] [-b bam -r reads -f genome/transriptome reference -d fast5dir -t threads -m bam_file_mapping_type] 
	
	it runs nanopolish eventalign; aggreagets current intensity values associated with single positions 
	
	nanopolish, samtools have to be installed and added to environmental paths!!

	-h [this help text]
	-b [bam file]
	-d [fast5 folder]
	-r [fastq reads file]
	-f [reference file used to produce bam file]
	-t [number of threads]
	-m [t: reads mapping performed using reference transcriptome; g: reads mapping performed with reference genome]
"

# example usage:
# sh ~/mygit/epinano_devel/v1.2/to_commit/misc/Epinano_Current.sh -b test_epinano_current_wrapper.bam -r test_epinano_current_wrapper.fastq -f ~/DATA_BACKUPS/ReferenceSequnces/Curlcakes/curlcake_constructs.fasta  -t 6  -m g -d /users/enovoa/hliu/basecalling/guppy3/rna_mod_vitro/RNAAB089716_unmodified/ 

if [[ "$#" -eq 0 ]]; then
	echo "$usage"
	exit 0
fi 

while getopts b:r:f:d:t:m:d:h option ; do 
	case "${option}" in	
		h) echo "$usage"
			exit ;;
		b) bam=${OPTARG};;
		r) reads=${OPTARG};;
		f) reference=${OPTARG};;
		d) fast5dir=${OPTARG};;
		t) threads=${OPTARG};;
		m) mapping_type=${OPTARG};;
		:) printf "missing argment for -%s\n" "$OPTARG" >&2
				printf "$usage" >&2
				exit 1 ;;
		
		\?) # incrorrect option
		echo "Error: invalid option: -%s\n" "$OPTARG" >&2
		exit 1 ;;
	
esac
done 
#shift $((OPTIND - 1))

#require: 
evaln=${bam/bam/eventalign.tsv}.gz 
src_d=$(dirname "$0") 

nanopolish index -d ${fast5dir} ${reads}
nanopolish eventalign --reads ${reads} --bam ${bam} --genome ${reference} --scale-events -t ${threads} -n --samples --signal-index  | pigz > ${evaln}

if [[ $mapping_type =~ 'g' ]]; then 
	samtools view -F3860 ${bam} | cut -f1 > ${bam}.forward.reads
	python ${src_d}/misc/eventalign_strandedness.py ${bam}.forward.reads ${evaln}
	fwrd=${evaln}.forward_strand.gz
	bwrd=${evaln}.reverse_strand.gz
	python ${src_d}/misc/Epinano_Current.py --infile ${evaln}.forward_strand.gz --reference ${reference} --outdir ${evaln}.forward  --threads ${threads} --strand + 
	python ${src_d}/misc/Epinano_Current.py --infile ${evaln}.reverse_strand.gz --reference ${reference} --outdir ${evaln}.reverse  --threads ${threads} --strand -
	python ${src_d}/misc/concat_events.py ${evaln}.forward_events.collapsed
	python ${src_d}/misc/concat_events.py ${evaln}.reverse_events.collapsed
	python ${src_d}/misc/Slide_Intensity.py --intensity_table ${evaln}.reverse_events.collapsed/Intensity.collapsed.tsv --window 5
	python ${src_d}/misc/Slide_Intensity.py --intensity_table ${evaln}.forward_events.collapsed/Intensity.collapsed.tsv --window 5 
fi

if [[ $mapping_type =~ 't' ]]; then 
	python ${src_d}/misc/Epinano_Current.py --infile ${evaln} --reference ${reference} --outdir ${evaln}.forward  --threads ${threads} --strand + 
	python ${src_d}/misc/concat_events.py ${evaln}.forward_events.collapsed 
	python ${src_d}/misc/Slide_Intensity.py --intensity_table ${evaln}.forward_events.collapsed/Intensity.collapsed.tsv --window 5 
fi 
