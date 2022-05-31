## run with docker container 
docker run -it -d --rm --name epivar2 -v "$PWD/xiao_data/":/project/ epi12 python3 /usr/local/bin/EpiNano/Epinano_Variants.py -R /project/ref.fa -b /project/reads.bam -s /usr/local/bin/EpiNano/misc/sam2tsv.jar -n 2

## use singularity and qsub on cluster
singularity pull docker://huanleliu/epi12

singularity exec -e /full/path/to/epi12_latest.sif python3 /usr/local/bin/EpiNano/Epinano_Variants.py -R /path/to/ref.fa -b /path/to/reads.bam -s /usr/local/bin/EpiNano/misc/sam2tsv.jar -n 2

