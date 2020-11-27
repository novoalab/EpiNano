#!/usr/bin/env bash

if [ ! -f "../../misc/sam2tsv.jar" ]; then
 echo "../../misc/sam2tsv.jar not found!"
 echo "please place/link sam2tsv.jar to that place"
 echo "reference file shoule be indexed with samtools faidx"
 echo "reference file shoule have a sequence dictionary file created with picard.jar"
 exit 1;
fi


echo "Basecalled reads and reads mapping commands are shown but are not going to be executed to save time"

echo "basecalling with guppy3.1.5"
echo "CMDs:"
echo "guppy_basecaller -c rna_r9.4.1_70bps_hac.cfg --compress_fastq -i ko_raw_fast5/ -r -s ko_fast5 --fast5_out --cpu_threads_per_caller 12"
echo "guppy_basecaller -c rna_r9.4.1_70bps_hac.cfg --compress_fastq -i wt_raw_fast5/ -r -s wt_fast5 --fast5_out --cpu_threads_per_caller 12"

echo "reads mapping either with minimap2 or graphmap2" 
echo "minimap2: minimap2 --MD -t 6 -ax map-ont ref.fa ko.fastq | samtools view -hbS -F 3844 - | samtools sort -@ 6 - ko"
echo "graphmap2:  graphmap align -r ref.fa -d ko.fastq -o ko.sam -v 1 -K fastq"
echo "graphmap2 with higher sensitivity: graphmap align -r ref.fa -d ko.fastq -o ko.sam -v 1 -K fastq  --rebuild-index --double-index --mapq -1 -x sensitive -z -1 --min-read-len 0 -A 7 -k 5"

echo "reads can also be mapped to reference genome with minimap2"


echo "compute varitants/error frequencies from bam file"
python ../../Epinano_Variants.py -R ref.fa -b wt.bam -n 6 -T t -s ../../misc/sam2tsv.jar
python ../../Epinano_Variants.py -R ref.fa -b ko.bam -n 6 -T t -s ../../misc/sam2tsv.jar

echo "predict based on deviance of mis" 
echo "AND"
echo "predict based on linear regression model residuals, using mis feature"
echo "CMD: Rscript ../Epinano_DiffErr.R -k ko.plus_strand.per.site.var.csv  -w wt.plus_strand.per.site.var.csv -t 5 -o HL -c 30 -f mis  -d 0.1 -p"
Rscript ../../Epinano_DiffErr.R -k ko.plus_strand.per.site.var.csv  -w wt.plus_strand.per.site.var.csv -t 5 -o HL -c 30 -f mis  -d 0.1 -p

echo "similarly we can use the same method but with sum_err"
echo "generate sum_err (mis, ins, del)"
python ../../misc/Epinano_sumErr.py --file ko.plus_strand.per.site.var.csv --out ko.sum_err.csv --kmer 0
python ../../misc/Epinano_sumErr.py --file wt.plus_strand.per.site.var.csv --out wt.sum_err.csv --kmer 0
Rscript ../../Epinano_DiffErr.R -k ko.sum_err.csv -w wt.sum_err.csv -d 0.1 -t 3 -p -o sumErr -f sum_err
echo "Switching on Epinano_sumErr.py -q will include quality score to compute sum_err"

echo "predict using pretrained SVM models"
echo "using q3,mis3,del3 features"
python ../../Epinano_Predict.py -o SVM_Predict -M ../models/rrach.q3.mis3.del3.linear.dump -p wt.plus_strand.per_site.5mer.csv -cl 8,13,23


echo "generate delta-features"
python ../../misc/Epinano_make_delta.py wt.plus_strand.per_site.5mer.csv ko.plus_strand.per_site.5mer.csv 5 > wt_ko_delta.5mer.csv
echo "predict using pretrained SVM models with delta features"
python ../../Epinano_Predict.py -o SVM_Predict_delta_features -M ../models/rrach.deltaQ3.deltaMis3.deltaDel3.linear.dump -p wt_ko_delta.5mer.csv -cl 7,12,22

echo "plot SVM-based prediction p-values"
Rscript ../../Epinano_Plot.R SVM_Predict_delta_features.mis3.del3.q3.MODEL.rrach.deltaQ3.deltaMis3.deltaDel3.linear.dump.csv
echo "../misc/Epinano_Plot.R can also generate plots from ../misc/Epinano_DiffErr.R outputs"

echo "Oh Wait! there is the current intensity values that you can incoporate into the above analyses"
echo "Extract and collapse current intensity values"
sh ../../Epinano_Current.sh -b ko.bam -r ko.fastq -f ref.fa -t 12 -m t -d ko_fast5/
sh ../../Epinano_Current.sh -b wt.bam -r wt.fastq -f ref.fa -t 12 -m t -d wt_fast5/

echo "create Intensity feature table contains Intensity and duration values"
python ../../misc/Slide_Intensity.py --intensity_table wt.eventalign.tsv.gz.forward_events.collapsed/Intensity.collapsed.tsv --window 5
python ../../misc/Slide_Intensity.py --intensity_table ko.eventalign.tsv.gz.forward_events.collapsed/Intensity.collapsed.tsv --window 5

echo "join current intensity features with error/variants features so that you can train model with interesting combinations of feature!"
python ../../misc/Join_variants_currents.py --variants wt.plus_strand.per_site.5mer.csv --intensity wt.eventalign.tsv.gz.forward_events.collapsed/Intensity.collapsed.tsv.5mer.csv --outfile wt.err.and.intensity.5mer.csv
python ../../misc/Join_variants_currents.py --variants ko.plus_strand.per_site.5mer.csv --intensity ko.eventalign.tsv.gz.forward_events.collapsed/Intensity.collapsed.tsv.5mer.csv --outfile ko.err.and.intensity.5mer.csv

echo "now it is possible to use both current intensity and Error features to do training or apply certain statistical tests to detect modifications!!"
