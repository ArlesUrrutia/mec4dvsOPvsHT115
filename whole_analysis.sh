#!/bin/bash

#define directories

maindir="$1"
sourcedir="$maindir/source"
readdir="$maindir/read"
bamdir="$maindir/bam"
tracksdir="$maindir/tracks"
countsdir="$maindir/counts"
resultsdir="$maindir/results"


#create directories if non existent
if [ ! -d $sourcedir ]; then
        mkdir -p $sourcedir
fi

if [ ! -d $readdir ]; then
        mkdir -p $readdir
fi

if [ ! -d $bamdir ]; then
        mkdir -p $bamdir
fi

if [ ! -d $tracksdir ]; then
        mkdir -p $tracksdir
fi

if [ ! -d $countsdir ]; then
        mkdir -p $countsdir
fi

if [ ! -d $resultsdir ]; then
        mkdir -p $resultsdir
fi


#set values for variables (only used for cuffdiff analysis)
name_base="IndexCHKPEI852160800"

#Choose which part of the analysis to carry
echo "Choose which part of the analysis to run:\n\n For the complete analysis type \"all\" \n For data download and quality control only type \"qc\" \n For the trimming only type \"trim\" \n For the mapping only type \"map\" \n For the genome browser tracks type \"tracks\" \n For the counts type \"count\" \n For the differential expression analysis and figures type \"de\" \n \n to run the whole analysis starting at different step type \"name_of_step-all\" (for example to do the whole analysis starting from the mapping type \"map-all\") \n" ;
read REPLY;

case $REPLY in
 *all )
follow="yes";;
 * )
follow="no";;
esac

case $REPLY in 
 all* )
start="1";;
 qc* )
start="1";;
 trim* )
start="2" ;;
 map* )
start="3";;
 tracks* )
start="4";;
 count* )
start="5";;
 de* )
start="6";;
esac


#get sample list from file

saved_IFS=$IFS
echo $IFS
IFS=$'\r\n' GLOBIGNORE='*' sample_list=$(cat name_samples.txt | awk '{print $1}')
IFS=$'\r\n' GLOBIGNORE=''
IFS=$saved_IFS

#version of genome, gtf and functional annotation (different version have different file names so just changing the version number might not work and some adjustment might be needed for the files to download)
vers="235"

#step 1 quatily control statistics
if [ $start =  1 ]
then

#Choose to download or give location of raw data
 echo "Choose to download the raw data or give the path to raw data directory by typing \"download\" or \"path/to/directory\") \n" ;
 read REPLY;
 
 case $REPLY in
  download )
 download="yes";;
  * )
 download="no"
 path_raw=$REPLY;;
 esac

#Choose to remove raw data
 echo "Do you want to remove raw data from working directory once processed? type \"yes\" or \"no\") \n" ;
 read REPLY;
 
 case $REPLY in
  yes )
 remove="yes";;
  * )
 remove="no";;
 esac

#change output file to avoid the already existing data to be written over
 mv nohup.out nohup.out_old

#get functional descripton of genes
 if [ ! -f "$sourcedir/c_elegans.PRJNA275000.WS254.functional_descriptions.txt" ]
 then
 wget -P $sourcedir ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/functional_descriptions/c_elegans.PRJNA275000.WS254.functional_descriptions.txt.gz
 gunzip $sourcedir/c_elegans.PRJNA275000.WS254.functional_descriptions.txt.gz
 fi
 
#get gtf
 if [ ! -f "$sourcedir/Caenorhabditis_elegans.WBcel${vers}.79.gtf" ]
 then
 wget -P $sourcedir http://ftp.ensembl.org/pub/release-79/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel${vers}.79.gtf.gz
 gunzip $sourcedir/Caenorhabditis_elegans.WBcel${vers}.79.gtf.gz
 fi 
#remove unnecessary lines that are creating issues in cufflinks
 cat $sourcedir/Caenorhabditis_elegans.WBcel${vers}.79.gtf | grep -P "ensembl\tCDS|ensembl\texon" > $sourcedir/c_elegans.gtf
#make simplified gtf table with gene length to be used by the R script
 cat $sourcedir/Caenorhabditis_elegans.WBcel${vers}.79.gtf | grep -P "gene\t" | awk '{print $10, $14, $7 , $4, $5 , $5-$4}' | sed 's/\"//g' | sed 's/\;//g' > $sourcedir/original_gtf_length.txt

#get genome sequence
#first remove old files
 if [ ! -f "$sourcedir/c_elegans.WS${vers}.genomic.fa" ]
 then
 rm $sourcedir/c_elegans.WS${vers}.genomic.fa.gz $sourcedir/c_elegans.WS${vers}.genomic.fa
 wget -P $sourcedir/ ftp://ftp.wormbase.org/pub/wormbase/releases/WS${vers}/species/c_elegans/c_elegans.WS${vers}.genomic.fa.gz
 gunzip $sourcedir/c_elegans.WS${vers}.genomic.fa.gz
 fi
#change fasta chromosome name format to match the one from the gtf and UCSC sequence
   sed -i 's/CHROMOSOME_// ' $sourcedir/c_elegans.WS${vers}.genomic.fa
fi

#get raw data
#downloading raw data
 if [ $download = yes ] 
 then
   for name in $sample_list
   do
    sample_nb=`echo $name | tail -c 3`
    sample1="$readdir/${name}_1.fq.gz"
    sample2="$readdir/${name}_2.fq.gz"
    dl_nb=`expr $sample_nb - 9`
 
    if [ ! -f $sample1 ]
    then
     wget  -P $readdir http://200.12.130.109/Secuenciacion/CGB1030/CGB1030-$dl_nb/${name}_1.fq.gz 
    fi 
    if [ ! -f $sample2 ]
    then
     wget -P $readdir http://200.12.130.109/Secuenciacion/CGB1030/CGB1030-$dl_nb/${name}_2.fq.gz 
    fi 

#do the quality control   
    echo "do quality control\n";
    nohup fastqc $sample1 
    nohup fastqc $sample2 
   done

 else
#copy data from internal source directory
   for name in $sample_list 
   do
    sample1="$readdir/${name}_1.fq.gz"
    sample2="$readdir/${name}_2.fq.gz"
    origin1="$path_raw/${name}_1.fq.gz"
    origin2="$path_raw/${name}_2.fq.gz"

    cp $origin1 $sample1
    cp $origin2 $sample2
    
#do the quality control   
    echo "do quality control\n";
    nohup fastqc $sample1 
    nohup fastqc $sample2 
   done
 fi
   
 start="done";
fi

#step 2 trimming
if  [ $start = done ] &&  [ $follow = yes  ] || [ $start =  2 ]  ;
then
echo "do trimming\n";
 for name in $sample_list 
 do
  read1="$readdir/${name}_1.fq.gz"
  read2="$readdir/${name}_2.fq.gz"
  trim1p="$readdir/${name}_1_p.fastq"  
  trim1u="$readdir/${name}_1_u.fastq"  
  trim2p="$readdir/${name}_2_p.fastq"  
  trim2u="$readdir/${name}_2_u.fastq"  

#Do trimming for paired end reads averaging across 4 base with an average quality required of 28 and the minimum lenght of 90bp for a read to be kept. phred + 33 quality score is used depending on the Illumina protocol 
  java -jar /home/carlos/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 $read1 $read2 $trim1p $trim1u $trim2p $trim2u SLIDINGWINDOW:4:28 MINLEN:90 
 
#redo quality control
   nohup fastqc $trim1p 
   nohup fastqc $trim2p 
   
 done
 
 mv nohup.out $readdir/fastqc_nohup_output

start="done";
fi

#step 3 mapping
if  [ $start = done ] &&  [ $follow = yes  ] || [ $start =  3 ]  ;
then
echo "do mapping\n";

#indexing of reference file with bowtie
 # ref="$sourcedir/ce_ws${vers}"
  ref="$sourcedir/c_elegans.WS${vers}.genomic"
  bowtie2-build -f ${ref}.fa $ref 

 for name in $sample_list 
 do
  gtf="$sourcedir/c_elegans.gtf"
  ref="$sourcedir/c_elegans.WS${vers}.genomic"
  trim1p="$readdir/${name}_1_p.fastq"  
  trim2p="$readdir/${name}_2_p.fastq"  
  output="$bamdir/${name}_thout"


#tophat mapping using 12 threads, a reference GTF file of known transcripts, the strand orientation of reads, and the indexed reference genome
   tophat2 -p 12 -G $gtf --library-type fr-firststrand -o $output $ref $trim1p $trim2p 
   if [ $remove = yes ]
   then
    rm $readdir/${name}_1.fq.gz  $readdir/${name}_2.fq.gz $trim1p $trim2p 
   fi

done


start="done";
fi

#step 4 tracks for data visualisation in UCSC genome browser
#for each sample, two track files are made with reads comming from the forward and reverse strand separately

if  [ $start = done ] &&  [ $follow = yes  ] || [ $start =  4 ]  ;
then
echo "do tracks\n";
 
 for name in $sample_list 
 do
  input="$bamdir/${name}_thout/accepted_hits.bam"
  ref="$sourcedir/c_elegans.WS${vers}.genomic.fa"
  tmpbam_f1="$bamdir/${name}_f1.bam"
  tmpbam_f2="$bamdir/${name}_f2.bam"
  tmpbam_r1="$bamdir/${name}_r1.bam"
  tmpbam_r2="$bamdir/${name}_r2.bam"
  tmpbam_r2="$bamdir/${name}_r2.bam"
  bam_fw="$bamdir/${name}_fw.bam"
  bam_re="$bamdir/${name}_re.bam"
  tmpbed_fw="$tracksdir/${name}_fw_tmp.bed"
  tmpbed_re="$tracksdir/${name}_re_tmp.bed"
  bed_fw="$tracksdir/${name}_fw.bed"
  bed_re="$tracksdir/${name}_re.bed"
 

  #first in pair, read reverse strand, transcript on forward strand in genome
   samtools view -bh -f 80 $input > $tmpbam_f1 

  #second in pair, read not on reverse strand, transcript on forward strand in genome
   samtools view -bh -f 128 -F 16 $input > $tmpbam_f2 
 
  #first in pair, read not on reverse strand, transcript on reverse strand in genome
   samtools view -bh -f 64 -F 16 $input > $tmpbam_r1 
 
  #second in pair, read on reverse strand, transcript on reverse strand in genome
   samtools view -bh -f 144 $input > $tmpbam_r2 
 
  #merge bam for reads mapping to transcripts on forward strand 
   samtools merge -f -h $tmpbam_f1 $bam_fw $tmpbam_f1 $tmpbam_f2 
   samtools index $bam_fw 
   rm $tmpbam_f1 $tmpbam_f2 
 
  #merge bam for reads mapping to transcripts on reverse strand 
   samtools merge -f -h $tmpbam_r1 $bam_re $tmpbam_r1 $tmpbam_r2 
   samtools index $bam_re 
   rm $tmpbam_r1 $tmpbam_r2 

  #make bedfiles 
   genomeCoverageBed -ibam $bam_fw -bg -g $ref > $tmpbed_fw 
   genomeCoverageBed -ibam $bam_re -bg -g $ref > $tmpbed_re 

  #add track specifications 
   echo "track type=\"bedGraph\" name=\"transcript forward $name\" color=0,0,128 visiblity=full" > $bed_fw 
   cat $tmpbed_fw >> $bed_fw
   sed -i '/^MtDNA/ d' $bed_fw
 
   echo "track type=\"bedGraph\" name=\"transcript reverse $name\" color=128,0,0 visiblity=full" > $bed_re 
   cat $tmpbed_re >> $bed_re
   sed -i '/^MtDNA/ d' $bed_re

   rm $tmpbed_fw $tmpbed_re $bam_fw $bam_re
 
 done

start="done";
fi

#step 5 counting
if  [ $start = done ] &&  [ $follow = yes  ] || [ $start =  5 ]  
then
 echo "do counting \n"
 rm $sourcedir/transcript_files.txt 
 touch $sourcedir/transcript_files.txt 
 for name in $sample_list 
 do
  input="$bamdir/${name}_thout/accepted_hits.bam"
  cuffout="$countsdir/${name}_clout"

#step 5a cufflinks counting and transcript reconstruction
  cufflinks -p 50 --library-type fr-firststrand -o $cuffout $input
  echo "$cuffout/transcripts.gtf \n" >> $sourcedir/transcript_files.txt 
 done

#prepare variable and sample ID/sample name correspondance for cuffdiff
 gtf="$sourcedir/c_elegans.gtf"
 ref="$sourcedir/c_elegans.WS${vers}.genomic.fa"
 OP50_12hA="$bamdir/${name_base}10_thout/accepted_hits.bam"
 OP50_12hB="$bamdir/${name_base}11_thout/accepted_hits.bam"
 OP50_12hC="$bamdir/${name_base}12_thout/accepted_hits.bam"
# HT115_12hA="$bamdir/${name_base}13_thout/accepted_hits.bam"
 HT115_12hB="$bamdir/${name_base}14_thout/accepted_hits.bam"
 HT115_12hC="$bamdir/${name_base}15_thout/accepted_hits.bam"
 OP50_24hA="$bamdir/${name_base}16_thout/accepted_hits.bam"
 OP50_24hB="$bamdir/${name_base}17_thout/accepted_hits.bam"
# OP50_24hC="$bamdir/${name_base}18_thout/accepted_hits.bam"
 HT115_24hA="$bamdir/${name_base}19_thout/accepted_hits.bam"
 HT115_24hB="$bamdir/${name_base}20_thout/accepted_hits.bam"
 HT115_24hC="$bamdir/${name_base}21_thout/accepted_hits.bam"
 OP50_48hA="$bamdir/${name_base}22_thout/accepted_hits.bam"
 OP50_48hB="$bamdir/${name_base}23_thout/accepted_hits.bam"
 OP50_48hC="$bamdir/${name_base}24_thout/accepted_hits.bam"
 HT115_48hA="$bamdir/${name_base}25_thout/accepted_hits.bam"
 HT115_48hB="$bamdir/${name_base}26_thout/accepted_hits.bam"
 HT115_48hC="$bamdir/${name_base}27_thout/accepted_hits.bam"

  cuffmerge  -g $gtf -o $countsdir/cuffmerge_out -p 50 -s $ref  $sourcedir/transcript_files.txt 

  cuffcompare -o $countsdir/cuffcompare -r $gtf $countsdir/cuffmerge_out/merged.gtf 
# cuffdiff 
  cuffdiff -o $resultsdir/cuff_12h -L HT115_12h,OP50_12h  $countsdir/cuffcompare.combined.gtf  $HT115_12hB,$HT115_12hC  $OP50_12hA,$OP50_12hB,$OP50_12hC 
  cuffdiff -o $resultsdir/cuff_24h -L HT115_24h,OP50_24h  $countsdir/cuffcompare.combined.gtf  $HT115_24hA,$HT115_24hB,$HT115_24hC  $OP50_24hA,$OP50_24hB 
  cuffdiff -o $resultsdir/cuff_48h -L HT115_48h,OP50_48h  $countsdir/cuffcompare.combined.gtf  $HT115_48hA,$HT115_48hB,$HT115_48hC  $OP50_48hA,$OP50_48hB,$OP50_48hC 

#step 5b HTseq counting using cufflinks transcripts
 #prepare the gtf file for HTseq count using the perl script to make an exon-based gtf file
  perl make_HT_gtf.pl $countsdir

 #do HTseq count
 for name in $sample_list 
 do
  input="$bamdir/${name}_thout/accepted_hits.bam"
  gtf="$countsdir/HTgtf.gtf"
  gtf_or="$sourcedir/Caenorhabditis_elegans.WBcel${vers}.79.gtf"
  htout="$countsdir/${name}.count"
  htout_or="$countsdir/${name}_original.count"
#do HTseq count with stranded paired end where first in pair is reverse to transcript orientation, with the more relaxed parameters, using original gtf file and the gtf from the cufflinks analysis
  htseq-count -f bam -s reverse -m intersection-nonempty --order=pos $input $gtf > $htout
  htseq-count -f bam -s reverse -m intersection-nonempty --order=pos $input $gtf_or > $htout_or

 done


 start="done";
fi

#step 6 R analysis
if  [ $start = done ] &&  [ $follow = yes  ] || [ $start =  6 ]  ;
then
 echo "do differential expression analysis\n"; 
 #make files for gene ID and sample names
 grep -vP "^#" $sourcedir/c_elegans.PRJNA275000.WS254.functional_descriptions.txt | head -1 | sed 's/ /\t/g' > $resultsdir/c_elegans.functional_descriptions.txt ; cat $sourcedir/c_elegans.PRJNA275000.WS254.functional_descriptions.txt | grep -P "^W" >> $resultsdir/c_elegans.functional_descriptions.txt
 cat $countsdir/cuffcompare.combined.gtf | awk '{print $16, $12, $12 "_" $14 , $4, $5}' | sed 's/\"//g' | sed 's/\;//g' > $resultsdir/cufflinks_transcripts.txt
#make file for gene ID sample names start end and length for TPM calculation for the exon-based gtf.
cat $countsdir/HTgtf.gtf | awk '{print $10, $16, $7 , $4, $5 , $5-$4}' | sed 's/\"//g' | sed 's/\;//g' > $resultsdir/HT_gtf_length.txt
#copy the rest of the needed files to the result directory
 cp heatmap_script.r name_samples.txt food_change_down_HT115.csv food_change_up_HT115.csv ortholist.csv c_elegans.PRJNA13758.WS257.geneIDs.txt regeneration_gene_nix.txt regeneration_gene_chen.txt  $resultsdir/
 cp $sourcedir/original_gtf_length.txt   $resultsdir/

#make summary for trimmed and mapped data
  out="$resultsdir/stat.txt"
  rm $out ; touch $out
for name in $sample_list
 do
  pre_trim1="$readdir/${name}_1.fq_fastqc/fastqc_data.txt"
  pre_trim2="$readdir/${name}_2.fq_fastqc/fastqc_data.txt"
  post_trim1="$readdir/${name}_1_p_fastqc/fastqc_data.txt"
  post_trim2="$readdir/${name}_2_p_fastqc/fastqc_data.txt"
  mapped="$bamdir/${name}_thout/align_summary.txt"
#output all information in file
  printf "$name\t" >> $out ; printf "pre-trimmed 1\t" >> $out ; cat $pre_trim1 | grep "Total Sequences" | tr -d '\n' >> $out ; printf "pre-trimmed 2\t" >> $out ; cat $pre_trim2 | grep "Total Sequences" | tr -d '\n' >> $out ; printf "post_trimmed 1\t" >> $out ; cat $post_trim1  | grep "Total Sequences" | tr -d '\n' >> $out  ; printf "post_trimmed 2\t" >> $out ; cat $post_trim2  | grep "Total Sequences" | tr -d '\n' >> $out ; printf "mapped\t" >> $out ; cat $mapped | grep "Aligned pairs"  >> $out

 done

#launch R script specifying the needed input/output directories
 Rscript de_analysis.r "$countsdir/ $resultsdir/"
 

start="done";
fi

echo "done \n";
