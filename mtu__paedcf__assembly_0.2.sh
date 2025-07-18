## assembly v0.3 :: The Morrigan

## setup 
  
  # separate setup and firing

# scp -i $key -oProxyJump=$jgary /Volumes/actimool/hpc-backup/ms__paedcf_raw.tar.gz $jdaed:/mnt/workspace2/jamie/


## getting off the ground!   ===================

mamba create -n climber -c bioconda -c conda-forge trimmomatic fastqc multiqc bowtie2 samtools hostile nonpareil -y   # 346MB
conda activate climber
# mamba activate k2 on geoff
# mamba install -c bioconda bracken kraken2 krakentools -y # >300MB ?
mamba create -n ccm ccmetagen -c bioconda -c conda-forge -y

## neo
# git bug - needed to do this in daed and pass it along
cd /claesson/jamie/bin ; git clone https://bitbucket.org/genomicepidemiology/kma.git
cd kma && make
# pip install ccm

# curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/GCF_000001405.40/download?filename=GCF_000001405.40.zip" -H "Accept: application/zip"
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5ss.fa.gz
# parallel gzip {} ::: $RAW/ref/*fna
# 
# time bowtie2-build --large-index --threads 14 \
#   $INOUT/GCF_000001405.40_GRCh38.p14_genomic.fna.gz,\
#   $INOUT/hs37d5ss.fa.gz \
#   $INOUT/mult_ref #> $INOUT/mult_ref.buildlog # 2>&1

# ccm_db=$DATA/db/ccmetagen    # neo
ccm_db=$WRK2/db/ccmetagen
mkdir -p $ccm_db ; cd $ccm_db ; curl "https://mediaflux.researchsoftware.unimelb.edu.au:443/mflux/share.mfjp?_token=Lqaic1pBmpDdqX8ofv1C1128247855&browser=true&filename=RefSeq_bf.zip" -d browser=false -o RefSeq_bf.zip

## one off, for below
# env_parallel --install


## =============================================================================

## dada daedalus
RAW=/mnt/workspace2/jamie/ms__paedcf_raw
WRK=/mnt/workspace2/jamie/ms__paedcf
## an morrigan
# RAW=/srv/jamie/raw_data/ms__paedcf_raw
# WRK=/srv/jamie/mtu__paedcf
## GEOFFREY
# RAW=~/data/ms__paedcf_raw
# WRK=~/data/ms__paedcf
## adapt
WRK=$SPN/mtu__paedcf

MAT=$WRK/Materials

DB=$WRK2/db
QC=$WRK/1__qc
FILT=$WRK/2__filt
HOST=$WRK/3__hostile
KRAK=$WRK/4__krak2
CCM=$WRK/4__ccmetagen

TEST=SC144C1-BAL_S45
TEST2=SC121C5-NS_S74
TEST3=/mnt/workspace2/jamie/db/ncbi_refseq__reference_genome/005/GCF_900451005.1_59123_A01/GCF_900451005.1_59123_A01_genomic.fna.gz


mkdir -p $RAW $WRK $QC $FILT $HOST $KRAK $MAT $CCM
mkdir -p $QC/mtu__paedcf_raw $QC/mtu__paedcf_filt $QC/mtu__paedcf_raw_multi $QC/mtu__paedcf_filt_multi

ls $RAW/*gz | sed -r 's/.*\/(.*)_R.*/\1/g' | sort -u > $MAT/mtu__paedcf__samples.txt


## wrangle   ==============================================================

cat $MAT/mtu__paedcf__samples.txt | parallel -j 12 "zcat $RAW/{}_L00*_R1_001.fastq.gz | gzip > $RAW/{}_R1.fastq.gz && rm $RAW/{}_L00*_R1_001.fastq.gz"
cat $MAT/mtu__paedcf__samples.txt | parallel -j 12 "zcat $RAW/{}_L00*_R2_001.fastq.gz | gzip > $RAW/{}_R2.fastq.gz && rm $RAW/{}_L00*_R2_001.fastq.gz"


## F / M Q C    ================================================================

fastqc -t 20 $RAW/*fastq.gz -o $QC/mtu__paedcf_raw && multiqc $QC/mtu__paedcf_raw -o $QC/mtu__paedcf_raw_multi
  

# delete negative controls if empty   ------------------------------------------
ls $RAW/*gz | sed -r 's/.*\/(.*)_R.*/\1/g' | sort -u  | grep -vE 'LIBNEG' > $MAT/mtu__paedcf__samples.txt
# adapt_GROVE ad hoc
mkdir $MAT
ls $HOST/*/*trimm.clean_1.fastq.gz | sed -r 's/.*\/(.*)_R.*/\1/g' | sort -u  | grep -vE 'LIBNEG' > $MAT/mtu__paedcf__samples.txt

  
##  t r i m m o   ==============================================================

echo '>Ampli_Tru_Seq_adapter__MOST_IMPORTANT_FEEEL_ME
CTGTCTCTTATACACATCT
>Transposase_Adap__for_tagmentation_1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Transposase_Adap__for_tagmentation_2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>PCR_primer_index_1
CAAGCAGAAGACGGCATACGAGATNNNNNNNGTCTCGTGGGCTCGG
>PCR_primer_index_2
AATGATACGGCGACCACCGAGATCTACACNNNNNTCGTCGGCAGCGTC
>TruSeq_single_index_LT_CD_HT__1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>TruSeq_single_index_LT_CD_HT__2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>PCR-Free_Prep__Tagm__additional_seq
ATGTGTATAAGAGACA
>Ampli_Tru_Seq_adapter__MOST_IMPORTANT_FEEEL_ME__mate_RC
AGATGTGTATAAGAGACAG
>Transposase_Adap__for_tagmentation_1.RC:
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>Transposase_Adap__for_tagmentation_2.RC:
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
>PCR_primer_index_1_RC:
CCGAGCCCACGAGACNNNNNNNATCTCGTATGCCGTCTTCTGCTTG
>PCR_primer_index_2_RC:
GACGCTGCCGACGANNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT
>TruSeq_single_index_LT_CD_HT__1_RC:
TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>TruSeq_single_index_LT_CD_HT__2_RC:
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PCR-Free_Prep__Tagm__additional_seq_RC:
TGTCTCTTATACACAT
>polG_just_from_PCF_concerns
GGGGGGGGGGGGGGGGGGGGGGGG
>polA_just_from_PCF_concerns
AAAAAAAAAAAAAAAAAAAAAAAA'> $MAT/fqc_trimmo_ill_ref.fa

  
##  + + +   MAKE SURE you parameterise via F/MQC  + + +   

# checking on TEST shows that removal of adapters also deals with a lot of 3' distortion. 
# Addition of AA/GGG to hopefully remove polys does little
## INSTEAD, ensure illumina first such that theres more pattern to be recognised - PDF ntoes this!
  
cat $MAT/mtu__paedcf__samples.txt | parallel -j 4 "trimmomatic PE \
  $RAW/{}_R1.fastq.gz \
  $RAW/{}_R2.fastq.gz \
  $FILT/{}_R1_trimm.fastq.gz \
  $FILT/{}_R1_trimm_unpaired.fastq.gz \
  $FILT/{}_R2_trimm.fastq.gz \
  $FILT/{}_R2_trimm_unpaired.fastq.gz \
  ILLUMINACLIP:$MAT/fqc_trimmo_ill_ref.fa:2:30:10:5 \
  SLIDINGWINDOW:6:15 \
  CROP:142 \
  HEADCROP:18 \
  MINLEN:120 \
  -threads 4 > $FILT/{}_trim.log 2>&1"

  # cat $FILT/*log | grep "Input Read Pairs" | sed -r "s/.*Read Pairs: ([0-9]*) Both Surviving: ([0-9]*) .*ward Only Surviving: ([0-9]*) .*erse Only Surviving: ([0-9]*).*pped: ([0-9]*) .*/\1\t\2\t\3\t\4\t\5/g"
## do with what yiou will  
  # ggplot(
  #   (reshape2::melt(read.table(file="input/mtu__paedcf__trimmo_counts.csv", header = FALSE, sep=","))),
  #   aes( y = value, x = variable)
  # ) + 
  # geom_line( )

mkdir $FILT/unpaired ; mv $FILT/*unpaired* $FILT/unpaired/
fastqc -t 12 $FILT/*_trimm.fastq.gz -o $QC/mtu__paedcf_filt ; multiqc $QC/mtu__paedcf_filt -o $QC/mtu__paedcf_filt_multi


##   dirty cleaning by hand - bt2 hum_deco   ==================================

bt2_db=$WRK2/db/bt2/hum_deco
deco=$WRK/3__bt2 ; mkdir $deco
bt2_threads=14

# for samp in $TEST ; done
while read samp ;
do  
  echo " + + +    $samp    -----------------------------------------------------------" ; 
  time bowtie2 -p $bt2_threads -x $bt2_db -1 $FILT/${samp}_R1_trimm.fastq.gz -2 $FILT/${samp}_R2_trimm.fastq.gz -S ${deco}/${samp}_bt2_refmapped.sam ;
  samtools view -bS ${deco}/${samp}_bt2_refmapped.sam > ${deco}/${samp}_bt2_refmapped.bam ;
  samtools view -b -f 12 -F 256 ${deco}/${samp}_bt2_refmapped.bam > ${deco}/${samp}_bt2_decontamd.bam ;
  samtools sort -n -m 5G -@ 2 ${deco}/${samp}_bt2_decontamd.bam -o ${deco}/${samp}_bt2_decontamd_sort.bam > ${deco}/${samp}_bam_sort.out ;
  samtools fastq -@ 8 ${deco}/${samp}_bt2_decontamd_sort.bam -1 ${deco}/${samp}_bt2decon_R1.fastq.gz -2 ${deco}/${samp}_bt2decon_R2.fastq.gz -0 /dev/null -s /dev/null -n && \
     echo " + + +   sample ${samp} completed decontam " >> ${deco}/bt2_deco_okay.log ;
  echo -e " + + +    --------------------------------------------------------------------\n" ; 

done< $MAT/mtu__paedcf__samples.txt >> $deco/deco.log

rm $deco/*am    # crucially...


##   hostile human+argos   =====================================================
  
time cat $MAT/mtu__paedcf__samples.txt | parallel -j 4 "hostile clean \
  --fastq1 $FILT/{}_R1_trimm.fastq.gz \
  --fastq2 $FILT/{}_R2_trimm.fastq.gz \
  --output $HOST/{}_hostless \
  --threads 4 \
  --index human-t2t-hla-argos985 > $HOST/{}_hostless.log 2>&1"


##  K r a k e n 2   ============================================================

## build kraken2   -------------------------------------------------------------
  #
  #  default is built already.. 
  #
  
## March2025 - bracken issuesl, get index online https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20241228.tar.gz
lang=$DB/kraken2_standard_langm ; mkdir $lang ; cd $lang
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20241228.tar.gz ; tar -xzcvf k2_pluspf_20241228.tar.gz ; rm k2_pluspf_20241228.tar.gz
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20241228.tar.gz ; tar -xzcvf k2_standard_20241228.tar.gz ; rm k2_standard_20241228.tar.gz

while read samp ;
do time kraken2 --db $lang \
    $HOST/${samp}_hostless/${samp}_R1_trimm.clean_1.fastq.gz \
    $HOST/${samp}_hostless/${samp}_R2_trimm.clean_2.fastq.gz \
    --paired \
    --threads 14 \
    --confidence 0.1 \
    --gzip-compressed \
    --report-zero-counts \
    --minimum-hit-groups 5 \
    --minimum-base-quality 20 \
    --report $KRAK/${samp}_kraken2_report \
    --unclassified-out $KRAK/${samp}_kraken_unclass# \
    --output $KRAK/${samp}_kraken_output > $KRAK/${samp}_krak2.log && echo " + + +   sample ${samp} completed task" >> $KRAK/kraken_okay.log ; 
done< $MAT/mtu__paedcf__samples.txt


kreport2mpa.py -r $KRAK/${TEST}_kraken_report -o $KRAK/mtu__paedcf__${TEST}_kraken_mpa
grep -h '|s_' $KRAK/mtu__paedcf__${TEST}_kraken_mpa | cut -f 1 |   sort | uniq | sed 's/|/\t/g' > $KRAK/mtu__paedcf__krakenStnd_taxonomy.tsv
less -S $KRAK/mtu__paedcf__krakenStnd_taxonomy.tsv


##   B r a c k e n    ==========================================================
  
# # mamba local git 0
# alias "bracken=$WRK2/bin/Bracken/bracken"
# alias "bracken-build=$WRK2/bin/Bracken/bracken-build"
  
## build Bracken   ----------------------------------
BR_kmer=35    # this is the default kmer length of the Kraken2 DB on the HPC
BR_leng=150   # beggars v. choosers - lenght of the Langmead kmer_dist file
## prebuilt langemad instances
# time bracken-build -d /mnt/workspace2/jamie/db/kraken2_standard -k 35 -l 124 -t 14
## segmentation fault persists... stil takes 45 mins, but perjhaps is actually related to the mem...

BR_r=$BR_leng
BR_l=S
BR_t=50   # counts! not threads
for i in $( cat $MAT/mtu__paedcf__samples.txt );
do 
  bracken -d $DB/kraken2_standard_langm/ -i $KRAK/${i}_kraken2_report -o $KRAK/${i}.bracken -r $BR_r -l $BR_l -t $BR_t ;
done > $KRAK/mtu__paedcf_krak2_bracken.log                                                               
combine_bracken_outputs.py --files $KRAK/*.bracken -o $KRAK/mtu__paedcf_krakenStnd_abundances.tsv >> $KRAK/mtu__paedcf_krak2_bracken.log                                                               

  
##  a c t F r u g a l l y   ---  salvage the unpaired reads for K2   ===========
  
  # could question quality & utility of this

# cat unpaired, uneven reads
cat $MAT/mtu__paedcf__samples.txt | parallel -j 16 "zcat $FILT/unpaired/{}*_trimm_unpaired.fastq.gz > $FILT/unpaired/{}_R12_trimm_unpaired.fastq && gzip $FILT/unpaired/{}_R12_trimm_unpaired.fastq"
  
time cat $MAT/mtu__paedcf__samples.txt | parallel -j 4 "hostile clean \
  --fastq1 $FILT/unpaired/{}_R12_trimm_unpaired.fastq.gz \
  --output $HOST/{}_hostless_unp \
  --threads 4 \
  --index human-t2t-hla-argos985 > $HOST/{}_hostless_unpaired.log"  2>&1

mkdir ${KRAK}_unp
cp -r $DB/kraken2_standard_langm/*k2d /dev/shm/ ; vmtouch -t /dev/shm/*.k2d 
while read samp ;
do time kraken2 --db /dev/shm \
    $HOST/${samp}_hostless/${samp}_R12_trimm_unpaired.clean.fastq.gz \
    --confidence 0.1 \
    --memory-mapping \
    --minimum-hit-groups 5 \
    --minimum-base-quality 20 \
    --threads 14 \
    --gzip-compressed \
    --report ${KRAK}_unp/${samp}_kraken2_unpaired_report \
    --unclassified-out ${KRAK}_unp/${samp}_kraken2_unpaired_unclass# \
    --output ${KRAK}_unp/${samp}_kraken2_unpaired_output > ${KRAK}_unp/${samp}_krak2_unpaired.log 2>&1
done< $MAT/mtu__paedcf__samples.txt
  
for i in $( cat $MAT/mtu__paedcf__samples.txt );
do 
  bracken -d $DB/kraken2_standard_langm/ -i ${KRAK}_unp/${i}_kraken2_unpaired_report -o ${KRAK}_unp/${i}.bracken -r $BR_r -l $BR_l -t $BR_t ;
done > ${KRAK}_unp/mtu__paedcf_krak2_bracken.log
combine_bracken_outputs.py --files ${KRAK}_unp/*.bracken -o ${KRAK}_unp/mtu__paedcf_krakenStnd_abundances_unpaired.tsv >> ${KRAK}_unp/mtu__paedcf_krak2_bracken.log                                                               
  
mkdir $MAT/output
cp ${KRAK}_unp/mtu__paedcf_krakenStnd_abundances_unpaired.tsv $KRAK/mtu__paedcf_krakenStnd_abundances.tsv ${KRAK}*/mtu__paedcf__krakenStnd_taxonomy.tsv $MAT/output/
scp -r -o ProxyJump=$jgary $jdaed://mnt/workspace2/jamie/ms__paedcf/Materials/output ~/


##  K a i j u   ================================================================
  # ... 
  # ..
  # ...    mooted through peer review. 
  # .
  # .


##  m u d d y p l a n   =========================================================

    # (metaphlan)
    
mamba create -n mpa metaphlan -c bioconda -c conda-forge -y 

mpa_db=$WRK2/db/metaphlan4 ; mpa=$WRK/4__metaphlan4
mkdir $mpa_db $mpa
# metaphlan --install --bowtie2db $mpa_db

while read samp ;
do time metaphlan \
  $HOST/${samp}_hostless/${samp}_R1_trimm.clean_1.fastq.gz,$HOST/${samp}_hostless/${samp}_R2_trimm.clean_2.fastq.gz \
  --bowtie2out $mpa/${samp}_mpa4crying_bt2out \
  -o $mpa/${samp}_mpa4crying_output \
  --input_type fastq \
  --bowtie2db $mpa_db \
  --nproc 14 && echo " + + +   sample ${samp} completed sad metaphlan task" >> $mpa/mpa_okay.log
  
done< $MAT/mtu__paedcf__samples.txt >> $mpa/mpa.log 2>&1

## well
  # >  WARNING: MetaPhlAn did not detect any microbial taxa in the sample.
  # >  [  . . .  ]
  # >  WARNING: The metagenome profile contains clades that represent multiple species merged into a single representant.
  # >  An additional column listing the merged species is added to the MetaPhlAn output.


##   n o n p a r e i l   =======================================================

## benefit of sequences run a-ground
# fastq is recommended for kmer algorithm
nonp=$WRK/5__nonpareil ; mkdir $nonp

while read samp ;
do 
  nonpareil -s ${HOST}/${samp}_hostless/${samp}_R1_trimm.clean_1.fastq.gz -T kmer -f fastq -b ${nonp}/${samp}_hostile_nonpareil -t 14
  # nonpareil -s ${HOST}/${samp}_hostless_unp/${samp}_R1_trimm.clean_1.fastq.gz -T kmer -f fastq -b ${nonp}/${samp}_hostile-unp_nonpareil
  # nonpareil -s ${deco}/${samp}_bt2decon_R1.fastq.gz -T kmer -f fastq -b ${nonp}/${samp}_bt2_nonpareil
done< $MAT/mtu__paedcf__samples.txt


