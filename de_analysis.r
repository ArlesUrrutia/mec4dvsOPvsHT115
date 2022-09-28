#load libraries
library(genefilter) #for row standard deviation calculation
library(edgeR) #for the differential expression with edgeR
library( DESeq ) #for the differential expression with DESeq
library(xtable) #to make tables suitable for latex input


#define threshold for the differential expression analysis cutoffs
logFC_threshold <- 1
pVal_threshold <- 0.05


#----parse command line arguments and put input and output directory names in "dir" object
dir <-commandArgs(trailingOnly = TRUE)
dir <- unlist(strsplit (dir, " "))

#load script for heatmaps
path <- paste ( dir[2] , 'heatmap_script.r' , sep='' )
source(path)

#read gene table with gene length made earlier from the original gtf file
path <- paste ( dir[2] , 'original_gtf_length.txt' , sep='' )
original_gtf <- read.table( path , header = F)
colnames(original_gtf) <- c("gene", "g_name","orientation","start","end","length")
original_gtf$gene <- as.character(original_gtf$gene)
original_gtf$g_name <- as.character(original_gtf$g_name)

#read exon-based table with exon length made earlier from the HTseq gtf file 
path <- paste ( dir[2] , 'HT_gtf_length.txt' , sep='' )
HTseq_gtf <- read.table( path , header = F)
colnames(HTseq_gtf) <- c("gene", "g_name","orientation","start","end","length")
HTseq_gtf$gene <- as.character(HTseq_gtf$gene)
HTseq_gtf$g_name <- as.character(HTseq_gtf$g_name)

#make table of file and extentions for each of the count analysis (original genes and exons)
gtf_version <- data.frame(matrix(c("HTseq_gtf","original_gtf",".count","_original.count","","_original"), ncol =3, byrow=F))
gtf_version$X1 <- as.character(gtf_version$X1)
gtf_version$X2 <- as.character(gtf_version$X2)
gtf_version$X3 <- as.character(gtf_version$X3)

#get the sample lists (list for ID and list for names) and correspondance table from the file. Also get the number of samples
path <- paste ( dir[2] , 'name_samples.txt' , sep='' )
sample_table <- read.table(path, header = F)
#sample_I <- c("HT115_12h_a","HT115_12h_b","HT115_12h_c","OP50_12h_a","OP50_12h_b","OP50_12h_c","HT115_24h_a","HT115_24h_b","HT115_24h_c","OP50_24h_a","OP50_24h_b","OP50_24h_c","HT115_48h_a","HT115_48h_b","HT115_48h_c","OP50_48h_a","OP50_48h_b","OP50_48h_c")
samples_I <- as.character(as.vector(sample_table$V1))
samples_N <- as.character(as.vector(sample_table$V2))
smpl_nb <- length(samples_I)

#---for each count analysis (original genes and exons)
for (j in c(1,2))
{
#----read information from count file for each sample
 for (i in 1:smpl_nb )
 {
  #get the right HTseq output file and put data in table
  path <- paste (dir[1],samples_I[i] , gtf_version$X2[j], sep='' )
  table_s <- read.table (path , header = F)
  table_s <- table_s[1:(length(table_s$V1)-5),]
  colnames (table_s) [1] <- "gene"
  colnames (table_s) [2] <- "count"
  table_s$gene <- as.character(table_s$gene)

  #merge with gene table 
  gene_table <- get(gtf_version$X1[j])
  table_s <- merge(table_s,gene_table[,c(1,6)], by.table_s=gene , by.gene_table=gene)
 
#  #calculate TMP values (here the read length is set as 100)
#  T_smpl <- sum (( table_s$count * 100 ) / table_s$length )
#  table_s$TPM <- ( table_s$count * 100 * 1000000 )/ ( table_s$length * T_smpl )
#  N_smpl <- sum ( table_s$count )
# 
#----calculate TPM values for genes
 
  # calculate effective gene length (here all samples have a read length of 100)
  read_len <- 100
  table_s$eff_len <- table_s$length + read_len -1
 
  # calculate relative gene coverage
  T_smpl <- sum (( table_s$count ) / table_s$eff_len )
 
  # calculate total counts in library 
  N_smpl <- sum ( table_s$count )
 
  # calculate TPM  
  table_s$TPM <- ( table_s$count * 1000000 )/ ( table_s$eff_len * T_smpl )
 

 #put TPM and count results in different tables and replace the TPM/counts columns with sample name
  sub_TPM <- table_s[,c(1,5)]
  sub_count <- table_s[,c(1,2)]
  cur_sample <- paste (samples_N[i])
  colnames (sub_TPM) [2] <- paste ( cur_sample ,sep = "")
  colnames (sub_count) [2] <- paste ( cur_sample ,sep = "")
 
  if (i == 1)
  {
   #----create the global TPM and count tables using the first sample in the list
   table_count_tmp <- sub_count
   table_TPM_tmp <- sub_TPM
   assign(paste("table_count",gtf_version$X3[j],sep=""), table_count_tmp)
   assign(paste("table_TPM",gtf_version$X3[j],sep=""), table_TPM_tmp)
  }
  #---rename the sample tables with sample id
  assign(paste (cur_sample , "_count" , gtf_version$X3[j] , sep = ""), sub_count)
  assign(paste (cur_sample , "_TPM" , gtf_version$X3[j] , sep = ""), sub_TPM)
  assign(paste (cur_sample , "_N" , gtf_version$X3[j] , sep = ""), N_smpl)
 }
 
     #---------------------------------------TPM-----------------------------
 #----add the TMP value of the other samples in the list to the global TPM table
 #get the table for the right experiment (gene or exon based) that contains only the sample one so far into a temporary table
 table_TPM_tmp <- get(paste("table_TPM",gtf_version$X3[j], sep=""))
 #add the other samples to this temporary table
 for (i in 2:smpl_nb )
 {
  cur_sample <- paste (samples_N[i])
  current <- get (paste (cur_sample , "_TPM" , gtf_version$X3[j], sep = ""))
  table_TPM_tmp <- merge (table_TPM_tmp, current, by.table_TPM_tmp = gene, by.current  = gene )
 }
 #once all the samples are in, add gene names as row names and convert the temporary table in final table with proper name
 row.names(table_TPM_tmp) <- as.character(table_TPM_tmp$gene)
 table_TPM_tmp <- table_TPM_tmp[,c(2:(smpl_nb+1))]
 assign(paste("table_TPM",gtf_version$X3[j], sep=""), table_TPM_tmp)
 
     #---------------------------------------count-----------------------------
 #----add the count value of the other samples in the list to the global count table
 #get the table for the right experiment (gene or exon based) that contains only the sample one so far into a temporary table
 table_count_tmp <- get(paste("table_count",gtf_version$X3[j],sep=""))
 #add the other samples to this temporary table
 for (i in 2:smpl_nb )
 {
  cur_sample <- paste (samples_N[i])
  current <- get (paste (cur_sample , "_count" , gtf_version$X3[j] , sep = ""))
  table_count_tmp <- merge (table_count_tmp, current, by.table_count_tmp = gene, by.current  = gene )
 }
 #once all the samples are in, add gene names as row names and convert the temporary table in final table with proper name
 row.names(table_count_tmp) <- as.character(table_count_tmp$gene)
 table_count_tmp <- table_count_tmp[,c(2:(smpl_nb+1))]
 assign(paste("table_count",gtf_version$X3[j],sep=""), table_count_tmp)
 
}
 
#-----------save data in a file first checkpoint-----------
f <- paste ( dir[2] , 'savedData.RData' , sep='' )
save.image(f)


#make exon-based correlation heatmap for tpm for all samples
#make matrix with correlation coefficient
x <- cor(table_TPM)
#define color scale
colorscale <- colorRampPalette(c('darkblue','yellow'))(30)
#make plot
PlotName <- paste(dir[2] , "correlation_heatmap.pdf", sep = "")
pdf( PlotName , width=15 , height = 16)

par(mar=c(22,17,11,10) )
image(as.matrix(t(x)) , axes = FALSE , col = colorscale )
axis ( 1 , at = seq(0 , 1 , by = 1/(length(x[1,])-1) ) , tick = FALSE , labels = colnames(x) ,  las = 2 , cex.axis = 1.2 )
axis ( 2 , at = seq(0 , 1 , by = 1/(length(x[,1])-1) ) , tick = FALSE , labels = row.names(x) ,  las = 2 , cex.axis = 1.2)
for (i in c(1:15))
{
 axis ( 3 , at = 0.7+(i-1)/55 , labels = round(seq (min(x, na.rm = TRUE) , max(x, na.rm = TRUE) , by = (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))/14 ), digit = 2)[i] ,  las = 2 , cex.axis = 1.3 , line = 3 ,  col = colorscale[(2*i-1)] , lwd.ticks = 7)
}
 axis ( 3 , at = 0.825 , label = "correlation coefficient" ,  las = 1 , cex.axis = 1.3 , line = 0.5 ,  tick = F )
dev.off()


#make gene-based correlation heatmap for tpm for all samples 
#make matrix with correlation coefficient
x <- cor(table_TPM_original)
#define color scale
colorscale <- colorRampPalette(c('darkblue','yellow'))(30)
#make plot
PlotName <- paste(dir[2] , "correlation_heatmap_original.pdf", sep = "")
pdf( PlotName , width=15 , height = 16)

par(mar=c(22,17,11,10) )
image(as.matrix(t(x)) , axes = FALSE , col = colorscale )
axis ( 1 , at = seq(0 , 1 , by = 1/(length(x[1,])-1) ) , tick = FALSE , labels = colnames(x) ,  las = 2 , cex.axis = 1.2 )
axis ( 2 , at = seq(0 , 1 , by = 1/(length(x[,1])-1) ) , tick = FALSE , labels = row.names(x) ,  las = 2 , cex.axis = 1.2)
for (i in c(1:15))
{
 axis ( 3 , at = 0.7+(i-1)/55 , labels = round(seq (min(x, na.rm = TRUE) , max(x, na.rm = TRUE) , by = (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))/14 ), digit = 2)[i] ,  las = 2 , cex.axis = 1.3 , line = 3 ,  col = colorscale[(2*i-1)] , lwd.ticks = 7)
}
 axis ( 3 , at = 0.825 , label = "correlation coefficient" ,  las = 1 , cex.axis = 1.3 , line = 0.5 ,  tick = F )
dev.off()



##----------------Generate tracks---------------------------------
#
##get gene files 
##path <- paste (dir[1],samples_I[i] , gtf_version$X2[j], sep='' 
#path <- paste (dir[1], samples_s[1] , ".txt", sep='' )
#count_table <- read.table (path , header = F, sep="\t" ,  quote="")
#bed_genes <- count_table[,c(1,4,5,9,7)]
#bed_genes$V9 <- levels(count_table$V9)[as.numeric(count_table$V9)]
#bed_genes$V9 <- substr(count_table$V9,9,14)
#bed_genes$score <- "1000"
#bed_genes <- bed_genes[,c(1,2,3,4,6,5)]
#colnames(bed_genes) <- c("chrom","start","end","name","score","strand")
#bed_genes$strand <- as.character(bed_genes$strand)
##bed_genes$strand[bed_genes$strand == "+"] <- 1
##bed_genes$strand[bed_genes$strand == "-"] <- -1
#bed_genes <- as.data.frame(bed_genes)
#bed_genes$score <- as.integer(bed_genes$score)
##bed_genes$strand <- as.numeric(bed_genes$strand)
#bed_genes <- bed_genes[order(bed_genes$start),]
#
##get samples track 
#
#path <- paste (dir[4], "Embryos_0_2hr_r1_F.bg" , sep='' )
#sample1F <- read.table(path)
#path <- paste (dir[4], "Embryos_0_2hr_r1_R.bg" , sep='' )
#sample1R <- read.table(path)
#
#wsp_smpl1 <- sample1F[sample1F$V1 == "Chromosome" & sample1F$V2 >= 1023047 & sample1F$V3 <= 1023760,]
#
#total_wsp_smpl1 <- 0
#for (i in c(1:length(wsp_smpl1[,1])))
#{
#total_wsp_smpl1 <- total_wsp_smpl1+(wsp_smpl1[i,3]-wsp_smpl1[i,2])*wsp_smpl1[i,4]
#}
#mean_wsp_smpl1 <- total_wsp_smpl1/714
#mean_wsp_smpl1_tenth <- mean_wsp_smpl1/10
#
#plotGene(bed_genes$start, bed_genes$end, bed_genes$strand, bed_genes$chrom, bed_genes$name,"Chromosome" ,  928000, 940000, idCex = 3, geneIdPatDif=NULL, useAxis = FALSE)
#axis(side=1, at=seq(930000,938000, by=2000), labels=T, cex.axis=3, padj = 1)
#mtext("Wolbachia gene models",side=3, line=-4 , font=2, cex=2, adj = 0.04)
##title("B", cex.main = 4, adj = 0, line = -1)
##title("C", cex.main = 4, adj = 0, line = -1)
##labelgenome("",928000,940000, scaleadjust=-2)
#
#par(mar=c(4,4,3,4))
#plotBedgraph(sample1F,"Chromosome",928000,940000, transparency=0.50,range=c(0,150) , main="", cex.main=0.8, color="blue")
#plotBedgraph(sample1R,"Chromosome",928000,940000, transparency=0.50,range=c(0,150) , main="", cex.main=0.8, color="red", overlay = T)
#axis(side=2, at=seq(0,150,by=50),las=2,tcl=.2 , cex.axis=3)
#labelgenome("",928000,940000, scaleadjust=-2 , cex.axis = 3, padj = 1)
#mtext("RNA-seq depth\n ",side=2, line=3, cex=2)
#mtext("embryos 0-2 hours",side=3, line=-1 , font=2, cex=2, adj = 0.03)
#abline(mean_wsp_smpl1_tenth/2, 0 , lty=2 )
#


#read cuffdiff results

#read the cufflings gtf info
path = paste(dir[2], "cufflinks_transcripts.txt" , sep = "")
cuff_gtf <- read.table(path)
gene_transcript <- unique(cuff_gtf[,c(1,2)])
gene_transcript$V1 <- as.character(gene_transcript$V1)
gene_transcript$V2 <- as.character(gene_transcript$V2)

#read the cuffdiff output table and make a simplified table
path = paste(dir[2], "cuff_12h/isoform_exp.diff" , sep = "")
cuff_12_comp <- read.table(path, header = T)
simple_cuff_12_comp <- cuff_12_comp[,c(1,3,8,9,10,14)]
colnames(simple_cuff_12_comp)[3] <- "value_HT115"
colnames(simple_cuff_12_comp)[4] <- "value_OP50"

path = paste(dir[2], "cuff_24h/isoform_exp.diff" , sep = "")
cuff_24_comp <- read.table(path, header = T)
simple_cuff_24_comp <- cuff_24_comp[,c(1,3,8,9,10,14)]
colnames(simple_cuff_24_comp)[3] <- "value_HT115"
colnames(simple_cuff_24_comp)[4] <- "value_OP50"

path = paste(dir[2], "cuff_48h/isoform_exp.diff" , sep = "")
cuff_48_comp <- read.table(path, header = T)
simple_cuff_48_comp <- cuff_48_comp[,c(1,3,8,9,10,14)]
colnames(simple_cuff_48_comp)[3] <- "value_HT115"
colnames(simple_cuff_48_comp)[4] <- "value_OP50"

#make table column name test specific for further merging of the data
colnames(simple_cuff_12_comp)[c(3:6)] <- paste(colnames(simple_cuff_12_comp)[c(3:6)], "12h", sep="_")
colnames(simple_cuff_24_comp)[c(3:6)] <- paste(colnames(simple_cuff_24_comp)[c(3:6)], "24h", sep="_")
colnames(simple_cuff_48_comp)[c(3:6)] <- paste(colnames(simple_cuff_48_comp)[c(3:6)], "48h", sep="_")

#merge results for three tests in one table
cuffdiff_table <- merge(simple_cuff_12_comp, simple_cuff_24_comp, by.simple_cuff_12_comp = test_id, by.simple_cuff_24_comp = test_id )
cuffdiff_table <- merge(cuffdiff_table, simple_cuff_48_comp, by.cuffdiff_table = test_id, by.simple_cuff_48_comp = test_id )

#get genes DE in one of the time points
allDE_list <- cuffdiff_table[cuffdiff_table$significant_12h == "yes" | cuffdiff_table$significant_24h == "yes" | cuffdiff_table$significant_48h == "yes",1]
allDE_list <- as.vector(allDE_list)

#make table for the subset of DE genes
table_DEgenes_cuffdiff <- cuffdiff_table[cuffdiff_table$test_id %in% allDE_list,c(1,3,4,7,8,11,12)]
row.names(table_DEgenes_cuffdiff) <- table_DEgenes_cuffdiff$test_id
table_DEgenes_cuffdiff <- table_DEgenes_cuffdiff[,c(2:7)]

#make list of the DE genes for each time point
comp_12h <- cuffdiff_table[cuffdiff_table$significant_12h == "yes",1]
comp_12h <- as.vector(comp_12h)
comp_24h <- cuffdiff_table[cuffdiff_table$significant_24h == "yes",1]
comp_24h <- as.vector(comp_24h)
comp_48h <- cuffdiff_table[cuffdiff_table$significant_48h == "yes",1]
comp_48h <- as.vector(comp_48h)

#plot heatmap of cuffdiff results
PlotName <- paste( dir[2], "heatmap_cuffdiff_genes.pdf", sep = "" )
x <- zrow_funct(table_DEgenes_cuffdiff)
pdf(PlotName  , width=25, height= 150)
par(mar=c(5,5,7,40) ,oma = c(15, 50, 7 , 5 ) , xpd=TRUE)
image(as.matrix(t(x[c(length(x[,1]):1),])) , axes = FALSE  )
axis ( 1 , at = seq(0 , 1 , by = 1/(length(x[1,])-1) ) , tick = FALSE , labels = colnames(x) ,  las = 2 , cex.axis = 1.5 )
for (n in c(1:length(x[,1])))
{
axis ( 2 , at = 1 - ((n-1)/(length(x[,1]) -1))  , tick = FALSE , labels = row.names(x)[n] ,  las = 2 , cex.axis = 1.5 )
}
mtext(side = 4,  "12h" , line = 2.5 , las = 1  , at = 1.003 , cex=1.6)
mtext(side = 4,  "24h" , line = 8 , las = 1  , at = 1.003 , cex = 1.6)
mtext(side = 4,  "48h" , line = 14.5 , las = 1  , at = 1.003 , cex = 1.6)
for (n in c(1:length(x[,1])))
{
if ( row.names(x)[n] %in%  comp_12h )
 {
 points(1.25 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2)
 }
if ( row.names(x)[n] %in%  comp_24h )
 {
 points(1.54 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2)
 }
if ( row.names(x)[n] %in%  comp_48h )
 {
 points( 1.83 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2)
 }
}
for (i in c(1:12))
{
 axis ( 3 , at = 0.60+(i-1)/20 , labels = round(seq (min(x, na.rm = TRUE) , max(x, na.rm = TRUE) , by = (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))/11 ), digit = 2)[i] ,  las = 2 , cex.axis = 2.1 , line = 3.5 ,  col = heat.colors(12)[i] , lwd.ticks = 7)
}
 axis ( 3 , at = 0.80 , label = "z-score" ,  las = 1 , cex.axis = 2.3 , line = 0.5 ,  tick = F )
dev.off()


#make heatmap using HTseq count (exon based) for genes DE in cuffdiff
#make ID correspondance table for DE genes
gene_transcript_de <- gene_transcript[gene_transcript$V2 %in%  allDE_list,]
colnames(gene_transcript_de) <- c("g_name","transcript")
gene_transcript_de_allID <- merge(gene_transcript_de, HTseq_gtf[,c(1,2)], by.gene_transcript=g_name, by.HTseq_gtf=g_name, all = F)
#make the table for DE genes
table_TPM_de <- table_TPM[row.names(table_TPM) %in% gene_transcript_de_allID$gene,]

#get list of DE genes with the correct exon ID
comp_12h_gene <- gene_transcript_de_allID$gene[gene_transcript_de_allID$transcript %in% comp_12h]
comp_24h_gene <- gene_transcript_de_allID$gene[gene_transcript_de_allID$transcript %in% comp_24h]
comp_48h_gene <- gene_transcript_de_allID$gene[gene_transcript_de_allID$transcript %in% comp_48h]

#make plot
PlotName <- paste (dir[2], "heatmap_HTcount_cuffdiff_genes.pdf" , sep = "")
x <- zrow_funct(table_TPM_de)
pdf(PlotName  , width=25, height= 150)
par(mar=c(5,5,7,40) ,oma = c(15, 50, 7 , 5 ) , xpd=TRUE)
image(as.matrix(t(x[c(length(x[,1]):1),])) , axes = FALSE  )
axis ( 1 , at = seq(0 , 1 , by = 1/(length(x[1,])-1) ) , tick = FALSE , labels = colnames(x) ,  las = 2 , cex.axis = 1.5 )
for (n in c(1:length(x[,1])))
{
axis ( 2 , at = 1 - ((n-1)/(length(x[,1]) -1))  , tick = FALSE , labels = row.names(x)[n] ,  las = 2 , cex.axis = 1.5 )
}
mtext(side = 4,  "12h" , line = 4 , las = 1  , at = 1.003 , cex=1.6)
mtext(side = 4,  "24h" , line = 11 , las = 1  , at = 1.003 , cex = 1.6)
mtext(side = 4,  "48h" , line = 18 , las = 1  , at = 1.003 , cex = 1.6)
for (n in c(1:length(x[,1])))
{
if ( row.names(x)[n] %in%  comp_12h_gene  )
 {
 points(1.25 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2)
 }
if ( row.names(x)[n] %in%  comp_24h_gene  )
 {
 points(1.54 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2)
 }
if ( row.names(x)[n] %in%  comp_48h_gene  )
 {
 points( 1.83 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2)
 }
}
for (i in c(1:12))
{
 axis ( 3 , at = 0.60+(i-1)/20 , labels = round(seq (min(x, na.rm = TRUE) , max(x, na.rm = TRUE) , by = (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))/11 ), digit = 2)[i] ,  las = 2 , cex.axis = 2.1 , line = 3.5 ,  col = heat.colors(12)[i] , lwd.ticks = 7)
}
 axis ( 3 , at = 0.80 , label = "z-score" ,  las = 1 , cex.axis = 2.3 , line = 0.5 ,  tick = F )
dev.off()

#-----------save data in a file first checkpoint-----------
f <- paste ( dir[2] , 'savedData.RData' , sep='' )
save.image(f)


#DE experiment exon based
#define the descriptions/treatment and subset of data for each analysis
dataset_12h <- table_count[rowSums(table_count) > 0,grep("12h",  colnames(table_count) , invert = F)]
groups_12h <- c(rep(1,length(grep("OP50",  colnames(dataset_12h)))),rep(2,length(grep("HT115",  colnames(dataset_12h)))))
treatment_12h <- c(rep("OP50",length(grep("OP50",  colnames(dataset_12h)))),rep("HT115",length(grep("HT115",  colnames(dataset_12h)))))
dataset_24h <- table_count[rowSums(table_count) > 0,grep("24h",  colnames(table_count) , invert = F)]
groups_24h <- c(rep(1,length(grep("OP50",  colnames(dataset_24h)))),rep(2,length(grep("HT115",  colnames(dataset_24h)))))
treatment_24h <- c(rep("OP50",length(grep("OP50",  colnames(dataset_24h)))),rep("HT115",length(grep("HT115",  colnames(dataset_24h)))))
dataset_48h <- table_count[rowSums(table_count) > 0,grep("48h",  colnames(table_count) , invert = F)]
groups_48h <- c(rep(1,length(grep("OP50",  colnames(dataset_48h)))),rep(2,length(grep("HT115",  colnames(dataset_48h)))))
treatment_48h <- c(rep("OP50",length(grep("OP50",  colnames(dataset_48h)))),rep("HT115",length(grep("HT115",  colnames(dataset_48h)))))

#differential experiment with edgeR
#define the amount of low expression genes to discard (this is based on previous observation of the data in a plot of expression level against unadjusted P-values)
theta = 0.6
#at 12 hours
group <- factor(groups_12h)
design <- model.matrix(~group)
dge <- DGEList(counts=dataset_12h, group=group)
dge <- calcNormFactors(dge)
dge_et <- estimateCommonDisp(dge)
dge_et <- estimateTagwiseDisp(dge_et)
#prepare subset of data to be tested
rs = rowSums (dataset_12h)
use = (rs > quantile(rs, probs=theta))
dge_et_Filt = dge_et[ use, ]
#make test
DE_12h <- exactTest(dge_et_Filt)
#adjust p-values for multiple testing
nb_exons <- length(DE_12h$table$PValue)
DE_12h_adjPval <- data.frame(matrix(NA, nrow=nb_exons , ncol=3))
DE_12h_adjPval$X2 <- p.adjust(DE_12h$table$PValue , method="BH")
DE_12h_adjPval$X1 <- row.names(DE_12h)
DE_12h_adjPval$X3 <- DE_12h$table$logFC
colnames(DE_12h_adjPval) <- c("gene","PValue","logFC")

#at 24 hours
group <- factor(groups_24h)
design <- model.matrix(~group)
dge <- DGEList(counts=dataset_24h, group=group)
dge <- calcNormFactors(dge)
dge_et <- estimateCommonDisp(dge)
dge_et <- estimateTagwiseDisp(dge_et)
#prepare subset of data to be tested
rs = rowSums (dataset_24h)
use = (rs > quantile(rs, probs=theta))
dge_et_Filt = dge_et[ use, ]
#make test
DE_24h <- exactTest(dge_et_Filt)
#adjust p-values for multiple testing
nb_exons <- length(DE_24h$table$PValue)
DE_24h_adjPval <- data.frame(matrix(NA, nrow=nb_exons , ncol=3))
DE_24h_adjPval$X2 <- p.adjust(DE_24h$table$PValue , method="BH")
DE_24h_adjPval$X1 <- row.names(DE_24h)
DE_24h_adjPval$X3 <- DE_24h$table$logFC
colnames(DE_24h_adjPval) <- c("gene","PValue","logFC")

#at 48 hours
group <- factor(groups_48h)
design <- model.matrix(~group)
dge <- DGEList(counts=dataset_48h, group=group)
dge <- calcNormFactors(dge)
dge_et <- estimateCommonDisp(dge)
dge_et <- estimateTagwiseDisp(dge_et)
#prepare subset of data to be tested
rs = rowSums (dataset_48h)
use = (rs > quantile(rs, probs=theta))
dge_et_Filt = dge_et[ use, ]
#make test
DE_48h <- exactTest(dge_et_Filt)
#adjust p-values for multiple testing
nb_exons <- length(DE_48h$table$PValue)
DE_48h_adjPval <- data.frame(matrix(NA, nrow=nb_exons , ncol=3))
DE_48h_adjPval$X2 <- p.adjust(DE_48h$table$PValue , method="BH")
DE_48h_adjPval$X1 <- row.names(DE_48h)
DE_48h_adjPval$X3 <- DE_48h$table$logFC
colnames(DE_48h_adjPval) <- c("gene","PValue","logFC")


#differential experiment with DESeq

#at 12 hours
cds = newCountDataSet( dataset_12h, treatment_12h )
cds = estimateSizeFactors( cds )
cds = estimateDispersions( cds )
#prepare subset of data to be tested
rs = rowSums ( counts ( cds ))
use = (rs > quantile(rs, probs=theta))
cdsFilt = cds[ use, ]
res = data.frame(filterstat = rowMeans(counts(cds)),pvalue    = nbinomTest(cds, "OP50", "HT115" ),row.names = featureNames(cds) )
PlotName <- paste( dir[2], "deseq_dispersion_filter_12h.pdf", sep = "")
pdf(PlotName  , width=20, height= 10)
par(mfrow=c(1,2))
plotDispEsts( cds )
with(res,plot(rank(filterstat)/length(filterstat), -log10(pvalue.pval), pch=16, cex=0.45))
dev.off()
#make test
DE_12h_DESeq = nbinomTest( cdsFilt, "OP50", "HT115" )

#at 24 hours
cds = newCountDataSet( dataset_24h, treatment_24h )
cds = estimateSizeFactors( cds )
cds = estimateDispersions( cds )
#prepare subset of data to be tested
rs = rowSums ( counts ( cds ))
use = (rs > quantile(rs, probs=theta))
cdsFilt = cds[ use, ]
res = data.frame(filterstat = rowMeans(counts(cds)),pvalue    = nbinomTest(cds, "OP50", "HT115" ),row.names = featureNames(cds) )
PlotName <- paste( dir[2], "deseq_dispersion_filter_24h.pdf", sep = "")
pdf(PlotName  , width=10, height= 10)
par(mfrow=c(1,2))
plotDispEsts( cds )
with(res,plot(rank(filterstat)/length(filterstat), -log10(pvalue.pval), pch=16, cex=0.45))
dev.off()
#make test
DE_24h_DESeq = nbinomTest( cdsFilt, "OP50", "HT115" )

#at 48 hours
cds = newCountDataSet( dataset_48h, treatment_48h )
cds = estimateSizeFactors( cds )
cds = estimateDispersions( cds )
#prepare subset of data to be tested
rs = rowSums ( counts ( cds ))
use = (rs > quantile(rs, probs=theta))
cdsFilt = cds[ use, ]
res = data.frame(filterstat = rowMeans(counts(cds)),pvalue    = nbinomTest(cds, "OP50", "HT115" ),row.names = featureNames(cds) )
PlotName <- paste( dir[2], "deseq_dispersion_filter_48h.pdf", sep = "")
pdf(PlotName  , width=10, height= 10)
par(mfrow=c(1,2))
plotDispEsts( cds )
with(res,plot(rank(filterstat)/length(filterstat), -log10(pvalue.pval), pch=16, cex=0.45))
dev.off()
#make test
DE_48h_DESeq = nbinomTest( cds, "OP50", "HT115" )



#DE experiment gene based

#define the descriptions/treatment and subset of data for each analysis
dataset_12h_original <- table_count_original[rowSums(table_count_original) > 0,grep("12h",  colnames(table_count_original) , invert = F)]
groups_12h_original <- c(rep(1,length(grep("OP50",  colnames(dataset_12h_original)))),rep(2,length(grep("HT115",  colnames(dataset_12h_original)))))
treatment_12h_original <- c(rep("OP50",length(grep("OP50",  colnames(dataset_12h_original)))),rep("HT115",length(grep("HT115",  colnames(dataset_12h_original)))))
dataset_24h_original <- table_count_original[rowSums(table_count_original) > 0,grep("24h",  colnames(table_count_original) , invert = F)]
groups_24h_original <- c(rep(1,length(grep("OP50",  colnames(dataset_24h_original)))),rep(2,length(grep("HT115",  colnames(dataset_24h_original)))))
treatment_24h_original <- c(rep("OP50",length(grep("OP50",  colnames(dataset_24h_original)))),rep("HT115",length(grep("HT115",  colnames(dataset_24h_original)))))
dataset_48h_original <- table_count_original[rowSums(table_count_original) > 0,grep("48h",  colnames(table_count_original) , invert = F)]
groups_48h_original <- c(rep(1,length(grep("OP50",  colnames(dataset_48h_original)))),rep(2,length(grep("HT115",  colnames(dataset_48h_original)))))
treatment_48h_original <- c(rep("OP50",length(grep("OP50",  colnames(dataset_48h_original)))),rep("HT115",length(grep("HT115",  colnames(dataset_48h_original)))))

#differential experiment with edgeR
#define the amount of low expression genes to discard (this is based on previous observation of the data in a plot of expression level against unadjusted P-values)
theta = 0.4
#at 12 hours
group <- factor(groups_12h_original)
design <- model.matrix(~group)
dge <- DGEList(counts=dataset_12h_original, group=group)
dge <- calcNormFactors(dge)
dge_et <- estimateCommonDisp(dge)
dge_et <- estimateTagwiseDisp(dge_et)
#prepare subset of data to be tested
rs = rowSums (dataset_12h_original)
use = (rs > quantile(rs, probs=theta))
dge_et_Filt = dge_et[ use, ]
#make test
DE_12h_original <- exactTest(dge_et_Filt)
#adjust p-values for multiple testing
nb_exons <- length(DE_12h_original$table$PValue)
DE_12h_original_adjPval <- data.frame(matrix(NA, nrow=nb_exons , ncol=3))
DE_12h_original_adjPval$X2 <- p.adjust(DE_12h_original$table$PValue , method="BH")
DE_12h_original_adjPval$X1 <- row.names(DE_12h_original)
DE_12h_original_adjPval$X3 <- DE_12h_original$table$logFC
colnames(DE_12h_original_adjPval) <-  c("gene","PValue","logFC")

#at 24 hours
group <- factor(groups_24h_original)
design <- model.matrix(~group)
dge <- DGEList(counts=dataset_24h_original, group=group)
dge <- calcNormFactors(dge)
dge_et <- estimateCommonDisp(dge)
dge_et <- estimateTagwiseDisp(dge_et)
#prepare subset of data to be tested
rs = rowSums (dataset_24h_original)
use = (rs > quantile(rs, probs=theta))
dge_et_Filt = dge_et[ use, ]
#make test
DE_24h_original <- exactTest(dge_et_Filt)
#adjust p-values for multiple testing
nb_exons <- length(DE_24h_original$table$PValue)
DE_24h_original_adjPval <- data.frame(matrix(NA, nrow=nb_exons , ncol=3))
DE_24h_original_adjPval$X2 <- p.adjust(DE_24h_original$table$PValue , method="BH")
DE_24h_original_adjPval$X1 <- row.names(DE_24h_original)
DE_24h_original_adjPval$X3 <- DE_24h_original$table$logFC
colnames(DE_24h_original_adjPval) <-  c("gene","PValue","logFC")

#at 48 hours
group <- factor(groups_48h_original)
design <- model.matrix(~group)
dge <- DGEList(counts=dataset_48h_original, group=group)
dge <- calcNormFactors(dge)
dge_et <- estimateCommonDisp(dge)
dge_et <- estimateTagwiseDisp(dge_et)
#prepare subset of data to be tested
rs = rowSums (dataset_48h_original)
use = (rs > quantile(rs, probs=theta))
dge_et_Filt = dge_et[ use, ]
#make test
DE_48h_original <- exactTest(dge_et_Filt)
#adjust p-values for multiple testing
nb_exons <- length(DE_48h_original$table$PValue)
DE_48h_original_adjPval <- data.frame(matrix(NA, nrow=nb_exons , ncol=3))
DE_48h_original_adjPval$X2 <- p.adjust(DE_48h_original$table$PValue , method="BH")
DE_48h_original_adjPval$X1 <- row.names(DE_48h_original)
DE_48h_original_adjPval$X3 <- DE_48h_original$table$logFC
colnames(DE_48h_original_adjPval) <-  c("gene","PValue","logFC")


#differential experiment with DESeq

#at 12 hours
cds = newCountDataSet( dataset_12h_original, treatment_12h_original )
cds = estimateSizeFactors( cds )
cds = estimateDispersions( cds )
#prepare subset of data to be tested
rs = rowSums ( counts ( cds ))
use = (rs > quantile(rs, probs=theta))
cdsFilt = cds[ use, ]
res = data.frame(filterstat = rowMeans(counts(cds)),pvalue    = nbinomTest(cds, "OP50", "HT115" ),row.names = featureNames(cds) )
PlotName <- paste( dir[2], "deseq_dispersion_filter_12h_original.pdf", sep = "")
pdf(PlotName  , width=20, height= 10)
par(mfrow=c(1,2))
plotDispEsts( cds )
with(res,plot(rank(filterstat)/length(filterstat), -log10(pvalue.pval), pch=16, cex=0.45))
dev.off()
#make test
DE_12h_original_DESeq = nbinomTest( cdsFilt, "OP50", "HT115" )

#at 24 hours
cds = newCountDataSet( dataset_24h_original, treatment_24h_original )
cds = estimateSizeFactors( cds )
cds = estimateDispersions( cds )
#prepare subset of data to be tested
rs = rowSums ( counts ( cds ))
use = (rs > quantile(rs, probs=theta))
cdsFilt = cds[ use, ]
res = data.frame(filterstat = rowMeans(counts(cds)),pvalue    = nbinomTest(cds, "OP50", "HT115" ),row.names = featureNames(cds) )
PlotName <- paste( dir[2], "deseq_dispersion_filter_24h_original.pdf", sep = "")
pdf(PlotName  , width=10, height= 10)
par(mfrow=c(1,2))
plotDispEsts( cds )
with(res,plot(rank(filterstat)/length(filterstat), -log10(pvalue.pval), pch=16, cex=0.45))
dev.off()
#make test
DE_24h_original_DESeq = nbinomTest( cdsFilt, "OP50", "HT115" )

#at 48 hours
cds = newCountDataSet( dataset_48h_original, treatment_48h_original )
cds = estimateSizeFactors( cds )
cds = estimateDispersions( cds )
#prepare subset of data to be tested
rs = rowSums ( counts ( cds ))
use = (rs > quantile(rs, probs=theta))
cdsFilt = cds[ use, ]
res = data.frame(filterstat = rowMeans(counts(cds)),pvalue    = nbinomTest(cds, "OP50", "HT115" ),row.names = featureNames(cds) )
PlotName <- paste( dir[2], "deseq_dispersion_filter_48h_original.pdf", sep = "")
pdf(PlotName  , width=10, height= 10)
par(mfrow=c(1,2))
plotDispEsts( cds )
with(res,plot(rank(filterstat)/length(filterstat), -log10(pvalue.pval), pch=16, cex=0.45))
dev.off()
#make test
DE_48h_original_DESeq = nbinomTest( cds, "OP50", "HT115" )



#-----------save data in a file first checkpoint-----------
f <- paste ( dir[2] , 'savedData.RData' , sep='' )
save.image(f)


#make table with sample stats
path=paste(dir[2], "stat.txt", sep="")
stat_table <- read.table(path, sep="\t")
stat_table <- stat_table[,c(1,2,4,5,7,8,10,11,13,14,15)]
#fix column names
colnames(stat_table)[1] <- "sample_ID"
colnames(stat_table)[3] <- paste(unlist(strsplit(as.character(stat_table$V2[1]), " " )),collapse="_")
colnames(stat_table)[5] <- paste(unlist(strsplit(as.character(stat_table$V5[1]), " " )),collapse="_")
colnames(stat_table)[7] <- paste(unlist(strsplit(as.character(stat_table$V8[1]), " " )),collapse="_")
colnames(stat_table)[9] <- paste(unlist(strsplit(as.character(stat_table$V11[1]), " " )),collapse="_")
colnames(stat_table)[11] <- as.character(stat_table$V14[1])
stat_table <- stat_table[,c(1,3,5,7,9,11)]
stat_table$mapped <- as.character(stat_table$mapped)
for (i in c(1:length(stat_table$mapped)))
{
 stat_table$mapped[i] <- as.vector(unlist(strsplit(as.character(stat_table$mapped[i]), " " )))[4]
}
colnames(stat_table) <- gsub('-','_', colnames(stat_table))
stat_table$all_pre_trim <- (stat_table$pre_trimmed_1) + (stat_table$pre_trimmed_2)
stat_table$all_post_trim <- (stat_table$post_trimmed_1) + (stat_table$post_trimmed_2)

#take only samples used in our analysis using sample_table. This step is useful only if the R script and general script are not taking the sample list form the same source. If unaltered, both scripts use the same sample list and this step is unnecessary.
colnames(sample_table) <- c("sample_ID","sample_name")
stat_table <- merge(sample_table,stat_table, by.all = sample_ID)

#make percentage of trimmed and mapped
stat_table$percent_trimmed <- round( (stat_table$all_pre_trim-stat_table$all_post_trim)/stat_table$all_pre_trim*100 , digit = 2)
stat_table$mapped <- as.numeric(stat_table$mapped)
stat_table$percent_mapped <- (stat_table$mapped*2)/stat_table$all_post_trim*100

#make a table with counts,tpg, gene table, DEp-value, log2 fold change
#give specific column names to the exon/gene-based tables, to be able to have all of them in one same table
#First make table for exon-based analysis
DESeq_12h_exon <- DE_12h_DESeq[,c(1,6,8)]
colnames(DESeq_12h_exon) <- c("exon_name","log2FC_12h_DESeq_exon","padj_12h_DESeq_exon")
DESeq_24h_exon <- DE_24h_DESeq[,c(1,6,8)]
colnames(DESeq_24h_exon) <- c("exon_name","log2FC_24h_DESeq_exon","padj_24h_DESeq_exon")
DESeq_48h_exon <- DE_48h_DESeq[,c(1,6,8)]
colnames(DESeq_48h_exon) <- c("exon_name","log2FC_48h_DESeq_exon","padj_48h_DESeq_exon")

edgeR_12h_exon <- DE_12h_adjPval[,c(1,3,2)]
colnames(edgeR_12h_exon) <- c("exon_name","log2FC_12h_edgeR_exon","padj_12h_edgeR_exon")
edgeR_24h_exon <- DE_24h_adjPval[,c(1,3,2)]
colnames(edgeR_24h_exon) <- c("exon_name","log2FC_24h_edgeR_exon","padj_24h_edgeR_exon")
edgeR_48h_exon <- DE_48h_adjPval[,c(1,3,2)]
colnames(edgeR_48h_exon) <- c("exon_name","log2FC_48h_edgeR_exon","padj_48h_edgeR_exon")

table_TPM_name <- table_TPM
colnames(table_TPM_name) <- paste(colnames(table_TPM_name), "_TPM", sep = "")
main_table_exon <- cbind(table_count,table_TPM_name)
main_table_exon$exon_name <- row.names(main_table_exon)
exon_corr <- HTseq_gtf[,c(1,2)]
colnames(exon_corr) <- c("exon_name","gene_name")
main_table_exon <- merge(main_table_exon, exon_corr, by.all = exon_name)
colnames(main_table_exon)[c(2:33)] <- paste(colnames(main_table_exon)[c(2:33)], "_exon" , sep="")
main_table_exon <- merge(main_table_exon, DESeq_12h_exon, by.all = exon_name, all=T)
main_table_exon <- merge(main_table_exon, DESeq_24h_exon, by.all = exon_name, all =T)
main_table_exon <- merge(main_table_exon, DESeq_48h_exon, by.all = exon_name, all=T)
main_table_exon <- merge(main_table_exon, edgeR_12h_exon, by.all = exon_name, all = T)
main_table_exon <- merge(main_table_exon, edgeR_24h_exon, by.all = exon_name, all=T)
main_table_exon <- merge(main_table_exon, edgeR_48h_exon, by.all = exon_name, all=T)

#Then make table for gene-based analysis
DESeq_12h_gene <- DE_12h_original_DESeq[,c(1,6,8)]
colnames(DESeq_12h_gene) <- c("gene_id","log2FC_12h_DESeq_gene","padj_12h_DESeq_gene")
DESeq_24h_gene <- DE_24h_original_DESeq[,c(1,6,8)]
colnames(DESeq_24h_gene) <- c("gene_id","log2FC_24h_DESeq_gene","padj_24h_DESeq_gene")
DESeq_48h_gene <- DE_48h_original_DESeq[,c(1,6,8)]
colnames(DESeq_48h_gene) <- c("gene_id","log2FC_48h_DESeq_gene","padj_48h_DESeq_gene")

edgeR_12h_gene <- DE_12h_original_adjPval[,c(1,3,2)]
colnames(edgeR_12h_gene) <- c("gene_id","log2FC_12h_edgeR_gene","padj_12h_edgeR_gene")
edgeR_24h_gene <- DE_24h_original_adjPval[,c(1,3,2)]
colnames(edgeR_24h_gene) <- c("gene_id","log2FC_24h_edgeR_gene","padj_24h_edgeR_gene")
edgeR_48h_gene <- DE_48h_original_adjPval[,c(1,3,2)]
colnames(edgeR_48h_gene) <- c("gene_id","log2FC_48h_edgeR_gene","padj_48h_edgeR_gene")

table_TPM_gene_name <- table_TPM_original
colnames(table_TPM_gene_name) <- paste(colnames(table_TPM_gene_name), "_TPM", sep = "")
main_table_gene <- cbind(table_count_original,table_TPM_gene_name)
main_table_gene$gene_id <- row.names(main_table_gene)
gene_corr <- original_gtf[,c(1,2)]
colnames(gene_corr) <- c("gene_id","gene_name")
main_table_gene <- merge(main_table_gene, gene_corr, by.all = gene_id)
colnames(main_table_gene)[c(2:33)] <- paste(colnames(main_table_gene)[c(2:33)], "_gene" , sep="")
main_table_gene <- merge(main_table_gene, DESeq_12h_gene, by.all = gene_id, all=T)
main_table_gene <- merge(main_table_gene, DESeq_24h_gene, by.all = gene_id, all=T)
main_table_gene <- merge(main_table_gene, DESeq_48h_gene, by.all = gene_id, all=T)
main_table_gene <- merge(main_table_gene, edgeR_12h_gene, by.all = gene_id, all=T)
main_table_gene <- merge(main_table_gene, edgeR_24h_gene, by.all = gene_id, all=T)
main_table_gene <- merge(main_table_gene, edgeR_48h_gene, by.all = gene_id, all=T)

#merge tables from exon and gene
main_table_all <- merge(main_table_exon, main_table_gene, by.all = gene_name, all = T)

#output tables as tsv files 
TableName <- paste( dir[2], "table_exon_summary.tsv", sep="")
write.table ( main_table_exon , TableName , row.names = F , col.names = T , quote = F , sep = "\t")
TableName <- paste( dir[2], "table_gene_summary.tsv", sep="")
write.table ( main_table_gene , TableName , row.names = F , col.names = T , quote = F , sep = "\t")
TableName <- paste( dir[2], "table_all_summary.tsv", sep="")
write.table ( main_table_all , TableName , row.names = F , col.names = T , quote = F , sep = "\t")


f <- paste ( dir[2] , 'savedData.RData' , sep='' )
save.image(f)


#table for DE genes (gene-based analysis)
#make list for each analysis and if gene is up or down regulated in HT115 vs OP50
DE_12h_original_DESeq_noNA <- na.omit(DE_12h_original_DESeq)
DE_24h_original_DESeq_noNA <- na.omit(DE_24h_original_DESeq)
DE_48h_original_DESeq_noNA <- na.omit(DE_48h_original_DESeq)

DESeq_12h_original_up <- DE_12h_original_DESeq_noNA[round(DE_12h_original_DESeq_noNA$padj, digit = 2) <= pVal_threshold & DE_12h_original_DESeq_noNA$log2FoldChange > logFC_threshold,1]
DESeq_24h_original_up <- DE_24h_original_DESeq_noNA[round(DE_24h_original_DESeq_noNA$padj, digit = 2) <= pVal_threshold & DE_24h_original_DESeq_noNA$log2FoldChange > logFC_threshold,1]
DESeq_48h_original_up <- DE_48h_original_DESeq_noNA[round(DE_48h_original_DESeq_noNA$padj, digit = 2) <= pVal_threshold & DE_48h_original_DESeq_noNA$log2FoldChange > logFC_threshold,1]
DESeq_12h_original_down <- DE_12h_original_DESeq_noNA[round(DE_12h_original_DESeq_noNA$padj, digit = 2) <= pVal_threshold & DE_12h_original_DESeq_noNA$log2FoldChange < -logFC_threshold,1]
DESeq_24h_original_down <- DE_24h_original_DESeq_noNA[round(DE_24h_original_DESeq_noNA$padj, digit = 2) <= pVal_threshold & DE_24h_original_DESeq_noNA$log2FoldChange < -logFC_threshold,1]
DESeq_48h_original_down <- DE_48h_original_DESeq_noNA[round(DE_48h_original_DESeq_noNA$padj, digit = 2) <= pVal_threshold & DE_48h_original_DESeq_noNA$log2FoldChange < -logFC_threshold,1]

edgeR_12h_original_up <- DE_12h_original_adjPval[round(DE_12h_original_adjPval$PValue, digit = 2) <= pVal_threshold & DE_12h_original_adjPval$logFC > logFC_threshold,1]
edgeR_24h_original_up <- DE_24h_original_adjPval[round(DE_24h_original_adjPval$PValue, digit = 2) <= pVal_threshold & DE_24h_original_adjPval$logFC > logFC_threshold,1]
edgeR_48h_original_up <- DE_48h_original_adjPval[round(DE_48h_original_adjPval$PValue, digit = 2) <= pVal_threshold & DE_48h_original_adjPval$logFC > logFC_threshold,1]
edgeR_12h_original_down <- DE_12h_original_adjPval[round(DE_12h_original_adjPval$PValue, digit = 2) <= pVal_threshold & DE_12h_original_adjPval$logFC < -logFC_threshold,1]
edgeR_24h_original_down <- DE_24h_original_adjPval[round(DE_24h_original_adjPval$PValue, digit = 2) <= pVal_threshold & DE_24h_original_adjPval$logFC < -logFC_threshold,1]
edgeR_48h_original_down <- DE_48h_original_adjPval[round(DE_48h_original_adjPval$PValue, digit = 2) <= pVal_threshold & DE_48h_original_adjPval$logFC < -logFC_threshold,1]

genes_DE_any <- unique(c(DESeq_12h_original_up, DESeq_24h_original_up, DESeq_48h_original_up, DESeq_12h_original_down, DESeq_24h_original_down, DESeq_48h_original_down, edgeR_12h_original_up, edgeR_24h_original_up, edgeR_48h_original_up, edgeR_12h_original_down, edgeR_24h_original_down, edgeR_48h_original_down))

genes_DE_any_up <- unique(c(DESeq_12h_original_up, DESeq_24h_original_up, DESeq_48h_original_up, edgeR_12h_original_up, edgeR_24h_original_up, edgeR_48h_original_up))
genes_DE_any_down <- unique(c(DESeq_12h_original_down, DESeq_24h_original_down, DESeq_48h_original_down, DESeq_12h_original_down, DESeq_24h_original_down, DESeq_48h_original_down, edgeR_12h_original_down, edgeR_24h_original_down, edgeR_48h_original_down, edgeR_12h_original_down, edgeR_24h_original_down, edgeR_48h_original_down))

genes_DE_all_up <- unique(c(DESeq_12h_original_up,DESeq_24h_original_up,DESeq_48h_original_up))[unique(c(DESeq_12h_original_up,DESeq_24h_original_up,DESeq_48h_original_up)) %in% unique(c(edgeR_12h_original_up, edgeR_24h_original_up, edgeR_48h_original_up))]
genes_DE_all_down <- unique(c(DESeq_12h_original_down,DESeq_24h_original_down,DESeq_48h_original_down))[unique(c(DESeq_12h_original_down,DESeq_24h_original_down,DESeq_48h_original_down)) %in% unique(c(edgeR_12h_original_down, edgeR_24h_original_down, edgeR_48h_original_down))]

gene_name_any <- original_gtf$g_name[original_gtf$gene %in% genes_DE_any]
gene_name_any_up <- original_gtf$g_name[original_gtf$gene %in% genes_DE_any_up]
gene_name_any_down <- original_gtf$g_name[original_gtf$gene %in% genes_DE_any_down]
gene_name_all_up <- original_gtf$g_name[original_gtf$gene %in% genes_DE_all_up]
gene_name_all_down <- original_gtf$g_name[original_gtf$gene %in% genes_DE_all_down]

#make table of genes ordered by fold change, with the fold change value where the gene is found differentially expressed
length_table <- max(length(gene_name_any_up),length(gene_name_any_down))+2
summary_table <- data.frame(matrix(NA, nrow=length_table , ncol=8))
summary_table[1,] <- c("up in HT115","up in HT115","up in HT115","up in HT115","down in HT115","down in HT115","down in HT115","down in HT115")
summary_table[2,] <- c("gene_name", "12h", "24h", "48h", "gene_name", "12h", "24h", "48h")
summary_table$X1[c(3:(length(gene_name_any_up)+2))] <- gene_name_any_up
summary_table$X5[c(3:(length(gene_name_any_down)+2))] <- gene_name_any_down

for (i in c(3:(length(gene_name_any_up)+2)))
{
 gene <- original_gtf$gene[original_gtf$g_name %in% summary_table$X1[i]]
 if (gene %in% unique(c(edgeR_12h_original_up,DESeq_12h_original_up)) == TRUE) 
 {
  summary_table$X2[i] <- round(mean(main_table_gene$log2FC_12h_DESeq_gene[main_table_gene$gene_id == gene]), digit = 2)
 }
 if (gene %in% unique(c(edgeR_24h_original_up,DESeq_24h_original_up)) == TRUE) 
 {
  summary_table$X3[i] <- round(mean(main_table_gene$log2FC_24h_DESeq_gene[main_table_gene$gene_id == gene]), digit = 2)
 }
 if (gene %in% unique(c(edgeR_48h_original_up,DESeq_48h_original_up)) == TRUE) 
 {
  summary_table$X4[i] <- round(mean(main_table_gene$log2FC_48h_DESeq_gene[main_table_gene$gene_id == gene]), digit = 2)
 }
}

for (i in c(3:(length(gene_name_any_down)+2)))
{
 gene <- original_gtf$gene[original_gtf$g_name %in% summary_table$X5[i]]
 if (gene %in% unique(c(edgeR_12h_original_down,DESeq_12h_original_down)) == TRUE) 
 {
  summary_table$X6[i] <- round(mean(main_table_gene$log2FC_12h_DESeq_gene[main_table_gene$gene_id == gene]), digit = 2)
 }
 if (gene %in% unique(c(edgeR_24h_original_down,DESeq_24h_original_down)) == TRUE) 
 {
  summary_table$X7[i] <- round(mean(main_table_gene$log2FC_24h_DESeq_gene[main_table_gene$gene_id == gene]), digit = 2)
 }
 if (gene %in% unique(c(edgeR_48h_original_down,DESeq_48h_original_down)) == TRUE) 
 {
  summary_table$X8[i] <- round(mean(main_table_gene$log2FC_48h_DESeq_gene[main_table_gene$gene_id == gene]), digit = 2)
 }
}

summary_table_nocolname <- summary_table[c(3:length(summary_table$X1)),]
first_half <- rbind(summary_table[c(1,2),c(1:4)],summary_table_nocolname[order(summary_table_nocolname[,2] ,summary_table_nocolname[,3], summary_table_nocolname[,4], decreasing=TRUE) ,c(1:4)])
second_half <- rbind(summary_table[c(1,2),c(5:8)],summary_table_nocolname[order(summary_table_nocolname[,6] ,summary_table_nocolname[,7], summary_table_nocolname[,8], decreasing=TRUE) ,c(5:8)])
summary_table <- cbind(first_half,second_half)

summary_table[is.na(summary_table)] <- " "

TableName <- paste( dir[2], "summary_de_table.tsv", sep="")
write.table ( summary_table , TableName , row.names = F , col.names = F , quote = F , sep = "\t")

#make figures for DE exons
DE_12h_DESeq_noNA <- na.omit(DE_12h_DESeq)
DE_24h_DESeq_noNA <- na.omit(DE_24h_DESeq)
DE_48h_DESeq_noNA <- na.omit(DE_48h_DESeq)

DESeq_12h <- DE_12h_DESeq_noNA[round(DE_12h_DESeq_noNA$padj, digit = 2) <= pVal_threshold & abs(DE_12h_DESeq_noNA$log2FoldChange) > logFC_threshold,1]
DESeq_24h <- DE_24h_DESeq_noNA[round(DE_24h_DESeq_noNA$padj, digit = 2) <= pVal_threshold & abs(DE_24h_DESeq_noNA$log2FoldChange) > logFC_threshold,1]
DESeq_48h <- DE_48h_DESeq_noNA[round(DE_48h_DESeq_noNA$padj, digit = 2) <= pVal_threshold & abs(DE_48h_DESeq_noNA$log2FoldChange) > logFC_threshold,1]

edgeR_12h <- DE_12h_adjPval[round(DE_12h_adjPval$PValue, digit = 2) <= pVal_threshold & abs(DE_12h_adjPval$logFC) > logFC_threshold,1]
edgeR_24h <- DE_24h_adjPval[round(DE_24h_adjPval$PValue, digit = 2) <= pVal_threshold & abs(DE_24h_adjPval$logFC) > logFC_threshold,1]
edgeR_48h <- DE_48h_adjPval[round(DE_48h_adjPval$PValue, digit = 2) <= pVal_threshold & abs(DE_48h_adjPval$logFC) > logFC_threshold,1]


table_TPM_de <- table_TPM[row.names(table_TPM) %in% unique(c(DESeq_12h,DESeq_24h,DESeq_48h,edgeR_12h,edgeR_24h,edgeR_48h)),]
table_TPM_de <- cbind( table_TPM_de[,c(1:5)], " " = NA , table_TPM_de[,c(6:10)]," " = NA, table_TPM_de[,c(11:16)])
#labels_plot <- paste(row.names(table_TPM_de), HTseq_gtf$g_name[HTseq_gtf$gene %in% row.names(table_TPM_de)], sep = "   ")

PlotName <- paste (dir[2], "heatmap_HTcount_exons.pdf" , sep = "")
x <- zrow_funct(table_TPM_de)
pdf(PlotName  , width=20, height= 20)
par(mar=c(5,5,7,40) ,oma = c(17, 25, 7 , 5 ) , xpd=TRUE)
image(as.matrix(t(x[c(length(x[,1]):1),])) , axes = FALSE  )
axis ( 1 , at = seq(0 , 1 , by = 1/(length(x[1,])-1) ) , tick = FALSE , labels = colnames(x) ,  las = 2 , cex.axis = 1.5 )
for (n in c(1:length(x[,1])))
{
#axis ( 2 , at = 1 - ((n-1)/(length(x[,1]) -1))  , tick = FALSE , labels = row.names(x)[n] ,  las = 2 , cex.axis = 1.5 )
exon_name <- row.names(x)[n]
axis ( 2 , at = 1 - ((n-1)/(length(x[,1]) -1))  , tick = FALSE , labels = paste(exon_name, HTseq_gtf$g_name[HTseq_gtf$gene %in% exon_name], sep = "   ") ,  las = 2 , cex.axis = 1.5 )
}
mtext(side = 4,  "12h" , line = 4 , las = 1  , at = 1.04 , cex=1.6, col="blue")
mtext(side = 4,  "DESeq\n24h" , line = 11 , las = 1  , at = 1.04 , cex = 1.6, col = "blue")
mtext(side = 4,  "48h" , line = 18 , las = 1  , at = 1.04 , cex = 1.6, col="blue")
mtext(side = 4,  "12h" , line = 25 , las = 1  , at = 1.04 , cex=1.6)
mtext(side = 4,  "edgeR\n24h" , line = 32 , las = 1  , at = 1.04 , cex = 1.6)
mtext(side = 4,  "48h" , line = 39 , las = 1  , at = 1.04 , cex = 1.6)
for (n in c(1:length(x[,1])))
{
if ( row.names(x)[n] %in%  DESeq_12h  )
 {
 points(1.25 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2, col="blue")
 }
if ( row.names(x)[n] %in%  DESeq_24h  )
 {
 points(1.54 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2, col="blue")
 }
if ( row.names(x)[n] %in%  DESeq_48h  )
 {
 points( 1.83 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2, col="blue")
 }
if ( row.names(x)[n] %in%  edgeR_12h  )
 {
 points(2.12 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2)
 }
if ( row.names(x)[n] %in%  edgeR_24h  )
 {
 points(2.41 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2)
 }
if ( row.names(x)[n] %in%  edgeR_48h  )
 {
 points( 2.70 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2)
 }
}
for (i in c(1:12))
{
 axis ( 3 , at = 0.60+(i-1)/20 , labels = round(seq (min(x, na.rm = TRUE) , max(x, na.rm = TRUE) , by = (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))/11 ), digit = 2)[i] ,  las = 2 , cex.axis = 2.1 , line = 3.5 ,  col = heat.colors(12)[i] , lwd.ticks = 7)
}
 axis ( 3 , at = 0.80 , label = "z-score" ,  las = 1 , cex.axis = 2.3 , line = 0.5 ,  tick = F )
dev.off()


#make figures for DE genes

#make list for ordering genes by logFC (as seen in the table summary)

gene_list <- c(summary_table$X1[c(3:length(summary_table$X1))], summary_table$X5[c(3:length(summary_table$X5))])
gene_list <- gene_list[(gene_list == " ") == FALSE]
gene_list_corr <- original_gtf[original_gtf$g_name %in% gene_list, c(1,2)]
gene_list_corr <- gene_list_corr[match(gene_list, gene_list_corr$g_name),]


DE_12h_original_DESeq_noNA <- na.omit(DE_12h_original_DESeq)
DE_24h_original_DESeq_noNA <- na.omit(DE_24h_original_DESeq)
DE_48h_original_DESeq_noNA <- na.omit(DE_48h_original_DESeq)

DESeq_12h_original <- DE_12h_original_DESeq_noNA[round(DE_12h_original_DESeq_noNA$padj, digit =2) <= pVal_threshold & abs(DE_12h_original_DESeq_noNA$log2FoldChange) > logFC_threshold,1]
DESeq_24h_original <- DE_24h_original_DESeq_noNA[round(DE_24h_original_DESeq_noNA$padj, digit =2) <= pVal_threshold & abs(DE_24h_original_DESeq_noNA$log2FoldChange) > logFC_threshold,1]
DESeq_48h_original <- DE_48h_original_DESeq_noNA[round(DE_48h_original_DESeq_noNA$padj, digit =2) <= pVal_threshold & abs(DE_48h_original_DESeq_noNA$log2FoldChange) > logFC_threshold,1]

edgeR_12h_original <- DE_12h_original_adjPval[round(DE_12h_original_adjPval$PValue, digit =2) <= pVal_threshold & abs(DE_12h_original_adjPval$logFC) > logFC_threshold,1]
edgeR_24h_original <- DE_24h_original_adjPval[round(DE_24h_original_adjPval$PValue, digit =2) <= pVal_threshold & abs(DE_24h_original_adjPval$logFC) > logFC_threshold,1]
edgeR_48h_original <- DE_48h_original_adjPval[round(DE_48h_original_adjPval$PValue, digit =2) <= pVal_threshold & abs(DE_48h_original_adjPval$logFC) > logFC_threshold,1]


table_TPM_original_de <- table_TPM_original[row.names(table_TPM_original) %in% unique(c(DESeq_12h_original,DESeq_24h_original,DESeq_48h_original,edgeR_12h_original,edgeR_24h_original,edgeR_48h_original)),]
table_TPM_original_de <- cbind( table_TPM_original_de[,c(1:5)], " " = NA , table_TPM_original_de[,c(6:10)]," " = NA, table_TPM_original_de[,c(11:16)])
#labels_plot <- paste(row.names(table_TPM_original_de), original_gtf$g_name[original_gtf$gene %in% row.names(table_TPM_original_de)], sep = "   ")
table_TPM_original_de <- table_TPM_original_de[match(gene_list_corr$gene, row.names(table_TPM_original_de)),]

PlotName <- paste (dir[2], "heatmap_HTcount_gene.pdf" , sep = "")
x <- zrow_funct(table_TPM_original_de)
pdf(PlotName  , width=20, height= 32)
par(mar=c(5,5,7,40) ,oma = c(17, 25, 7 , 5 ) , xpd=TRUE)
image(as.matrix(t(x[c(length(x[,1]):1),])) , axes = FALSE  )
axis ( 1 , at = seq(0 , 1 , by = 1/(length(x[1,])-1) ) , tick = FALSE , labels = colnames(x) ,  las = 2 , cex.axis = 1.5 )
for (n in c(1:length(x[,1])))
{
#axis ( 2 , at = 1 - ((n-1)/(length(x[,1]) -1))  , tick = FALSE , labels = row.names(x)[n] ,  las = 2 , cex.axis = 1.5 )
gene_name <- row.names(x)[n]
axis ( 2 , at = 1 - ((n-1)/(length(x[,1]) -1))  , tick = FALSE , labels = paste(gene_name, original_gtf$g_name[original_gtf$gene %in% gene_name] , sep = "  ") ,  las = 2 , cex.axis = 1.5 )
}
mtext(side = 4,  "12h" , line = 4 , las = 1  , at = 1.04 , cex=1.6, col="blue")
mtext(side = 4,  "DESeq\n24h" , line = 11 , las = 1  , at = 1.04 , cex = 1.6, col = "blue")
mtext(side = 4,  "48h" , line = 18 , las = 1  , at = 1.04 , cex = 1.6, col="blue")
mtext(side = 4,  "12h" , line = 25 , las = 1  , at = 1.04 , cex=1.6)
mtext(side = 4,  "edgeR\n24h" , line = 32 , las = 1  , at = 1.04 , cex = 1.6)
mtext(side = 4,  "48h" , line = 39 , las = 1  , at = 1.04 , cex = 1.6)
for (n in c(1:length(x[,1])))
{
if ( row.names(x)[n] %in%  DESeq_12h_original  )
 {
 points(1.25 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2, col="blue")
 }
if ( row.names(x)[n] %in%  DESeq_24h_original  )
 {
 points(1.54 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2, col="blue")
 }
if ( row.names(x)[n] %in%  DESeq_48h_original  )
 {
 points( 1.83 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2, col="blue")
 }
if ( row.names(x)[n] %in%  edgeR_12h_original  )
 {
 points(2.12 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2)
 }
if ( row.names(x)[n] %in%  edgeR_24h_original  )
 {
 points(2.41 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2)
 }
if ( row.names(x)[n] %in%  edgeR_48h_original  )
 {
 points( 2.70 , (1 - ((n-1)/(length(x[,1]) -1))) , pch=19 , cex = 2)
 }
}
for (i in c(1:12))
{
 axis ( 3 , at = 0.60+(i-1)/20 , labels = round(seq (min(x, na.rm = TRUE) , max(x, na.rm = TRUE) , by = (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))/11 ), digit = 2)[i] ,  las = 2 , cex.axis = 2.1 , line = 3.5 ,  col = heat.colors(12)[i] , lwd.ticks = 7)
}
 axis ( 3 , at = 0.80 , label = "z-score" ,  las = 1 , cex.axis = 2.3 , line = 0.5 ,  tick = F )
dev.off()




#-----------save data in a file first checkpoint-----------
f <- paste ( dir[2] , 'savedData.RData' , sep='' )
save.image(f)




#make list for genes that show differential expression in C. elegans when changing the food source
#genes up in HT115

path <- paste ( dir[2] , 'food_change_up_HT115.csv' , sep='' )
food_up_HT115 <- read.csv(path, header=T)
colnames(food_up_HT115) <- c("gene_id","gene_name","seq_name","fold_change")

#genes down in HT115

path <- paste ( dir[2] , 'food_change_down_HT115.csv' , sep='' )
food_down_HT115 <- read.csv(path, header=T)
colnames(food_down_HT115) <- c("gene_id","gene_name","seq_name","fold_change")


#--------------Prepare table for ortholog genes------------------------

#table of ortholog genes
path <- paste ( dir[2] , 'ortholist.csv' , sep='' )
orthologs <- read.csv(path, header = T)
colnames(orthologs) <- c("sequence_ID","Human_gene_ID","identity_Ce_H","identity_H_Ce")

#table of sequence/gene Id correspondance
path <- paste ( dir[2] , 'c_elegans.PRJNA13758.WS257.geneIDs.txt' , sep='' )
correspondance <- read.csv(path, header= F)
colnames(correspondance) <- c("number","gene","gene_name","sequence_ID","lethal")
#fix sequence name
correspondance$sequence_ID <- as.character(correspondance$sequence_ID)
#correspondance[correspondance$gene_name == "unc-17", 4] <- "ZC416.8a"
correspondance[correspondance$gene_name == "cha-1", 4] <- "ZC416.8b"
correspondance[correspondance$gene_name == "unc-32", 4] <- "ZK637.8a"
correspondance[correspondance$gene_name == "dab-1", 4] <- "M110.5a"
correspondance[correspondance$gene_name == "usp-46", 4] <- "R10E11.3a"



#make table for DE genes
ortholog_de <- merge(gene_list_corr, correspondance, by.all = gene)
ortholog_de <- ortholog_de[,c(1,2,5)]
ortholog_de <- merge(ortholog_de, orthologs, by.all = sequence_ID, all.x = T)
ortholog_de <- ortholog_de[,c(1:4)]
ortholog_de  <- aggregate(ortholog_de[,c(1,3:4)], by=list(ortholog_de$gene), paste, collapse=",")
colnames(ortholog_de)[1] <- "gene"
ortholog_de$sequence_ID[ortholog_de$sequence_ID == ""] <- NA
ortholog_de$sequence_ID <- unique(unlist(strsplit(as.character(ortholog_de$sequence_ID), "," )))
ortholog_de$g_name <- unique(unlist(strsplit(as.character(ortholog_de$g_name), "," )))

#make table for gene function
path <- paste ( dir[2] , 'c_elegans.functional_descriptions.txt' , sep='' )
gene_function <- read.table(path, header= T , sep ="\t", comment.char = "", quote = "\n" )
colnames(gene_function)[1] <- "gene"
gene_function_de <- merge(ortholog_de, gene_function, by.all = gene)

#make table for writing
gene_orth_funct <- gene_function_de[,c(3,4,7)]
gene_orth_funct_up <- gene_orth_funct[gene_orth_funct$g_name %in% gene_name_any_up,]
colnames(gene_orth_funct_up) <- c("gene name","human ortholog","description")
table_name <- paste ( dir[2] , 'table_orth_funct.tex' , sep='' )
table_latex <- xtable(gene_orth_funct_up , type = "latex" , file = table_name , table.placement = "h" , caption.placement = "bottom"  )
print.xtable(table_latex , type = "latex" , file = table_name , table.placement = "h" , include.rownames = FALSE)

table_name <- paste ( dir[2] , 'table_orth_funct.tsv' , sep='' )
write.table ( gene_orth_funct_up , table_name , row.names = F , col.names = T , quote = F , sep = "\t")

#genes involved in axon regeneration
path <- paste ( dir[2] , 'regeneration_gene_nix.txt' , sep='' )
nix_list <- scan(file = path, what = "character", comment="#")
nix_part1 <- unlist(strsplit(nix_list[1], ","))
nix_part1 <- correspondance[correspondance$sequence_ID  %in% nix_part1, 2]
nix_part2 <- unlist(strsplit(nix_list[2], ","))
nix_part2 <- correspondance[correspondance$gene_name  %in% nix_part2, 2]
nix_list <- unique(c(as.character(nix_part1),as.character(nix_part2)))

path <- paste ( dir[2] , 'regeneration_gene_chen.txt' , sep='' )
chen_list <- scan(file = path, what = "character", comment="#")
chen_part1 <- unlist(strsplit(chen_list[1], ","))
chen_part1 <- original_gtf[original_gtf$g_name  %in% chen_part1, 1]
chen_part2 <- unlist(strsplit(chen_list[2], ","))
chen_part2 <- original_gtf[original_gtf$g_name  %in% chen_part2, 1]
chen_part3 <- unlist(strsplit(chen_list[3], ","))
chen_part3 <- correspondance[correspondance$sequence_ID  %in% chen_part3, 2]
chen_list <- unique(c(as.character(chen_part1),as.character(chen_part2),as.character(chen_part3)))

list_all_regeneration <- unique(c(nix_list,chen_list))


#---------Prepare list and latex variables----------------

#number of expressed genes
vect_gene <-vector()
for (i in c(1:length(table_TPM_original[1,]))) 
{
vect_gene <- c(vect_gene, length(table_TPM_original[!table_TPM_original == 0,i]))
}

#gene differentially expressed
DE_12h_original_DESeq_noNA <- na.omit(DE_12h_original_DESeq)
DE_24h_original_DESeq_noNA <- na.omit(DE_24h_original_DESeq)
DE_48h_original_DESeq_noNA <- na.omit(DE_48h_original_DESeq)

DESeq_12h_original_up <- DE_12h_original_DESeq_noNA[round(DE_12h_original_DESeq_noNA$padj, digit = 2) <= pVal_threshold & DE_12h_original_DESeq_noNA$log2FoldChange > logFC_threshold,1]
DESeq_24h_original_up <- DE_24h_original_DESeq_noNA[round(DE_24h_original_DESeq_noNA$padj, digit = 2) <= pVal_threshold & DE_24h_original_DESeq_noNA$log2FoldChange > logFC_threshold,1]
DESeq_48h_original_up <- DE_48h_original_DESeq_noNA[round(DE_48h_original_DESeq_noNA$padj, digit = 2) <= pVal_threshold & DE_48h_original_DESeq_noNA$log2FoldChange > logFC_threshold,1]
DESeq_12h_original_down <- DE_12h_original_DESeq_noNA[round(DE_12h_original_DESeq_noNA$padj, digit = 2) <= pVal_threshold & DE_12h_original_DESeq_noNA$log2FoldChange < -logFC_threshold,1]
DESeq_24h_original_down <- DE_24h_original_DESeq_noNA[round(DE_24h_original_DESeq_noNA$padj, digit = 2) <= pVal_threshold & DE_24h_original_DESeq_noNA$log2FoldChange < -logFC_threshold,1]
DESeq_48h_original_down <- DE_48h_original_DESeq_noNA[round(DE_48h_original_DESeq_noNA$padj, digit = 2) <= pVal_threshold & DE_48h_original_DESeq_noNA$log2FoldChange < -logFC_threshold,1]

edgeR_12h_original_up <- DE_12h_original_adjPval[round(DE_12h_original_adjPval$PValue, digit = 2) <= pVal_threshold & DE_12h_original_adjPval$logFC > logFC_threshold,1]
edgeR_24h_original_up <- DE_24h_original_adjPval[round(DE_24h_original_adjPval$PValue, digit = 2) <= pVal_threshold & DE_24h_original_adjPval$logFC > logFC_threshold,1]
edgeR_48h_original_up <- DE_48h_original_adjPval[round(DE_48h_original_adjPval$PValue, digit = 2) <= pVal_threshold & DE_48h_original_adjPval$logFC > logFC_threshold,1]
edgeR_12h_original_down <- DE_12h_original_adjPval[round(DE_12h_original_adjPval$PValue, digit = 2) <= pVal_threshold & DE_12h_original_adjPval$logFC < -logFC_threshold,1]
edgeR_24h_original_down <- DE_24h_original_adjPval[round(DE_24h_original_adjPval$PValue, digit = 2) <= pVal_threshold & DE_24h_original_adjPval$logFC < -logFC_threshold,1]
edgeR_48h_original_down <- DE_48h_original_adjPval[round(DE_48h_original_adjPval$PValue, digit = 2) <= pVal_threshold & DE_48h_original_adjPval$logFC < -logFC_threshold,1]

genes_DE_any <- unique(c(DESeq_12h_original_up, DESeq_24h_original_up, DESeq_48h_original_up, DESeq_12h_original_down, DESeq_24h_original_down, DESeq_48h_original_down, edgeR_12h_original_up, edgeR_24h_original_up, edgeR_48h_original_up, edgeR_12h_original_down, edgeR_24h_original_down, edgeR_48h_original_down))
genes_DE_12h_any <- unique(c(DESeq_12h_original_up, DESeq_12h_original_down, edgeR_12h_original_up, edgeR_12h_original_down))
genes_DE_24h_any <- unique(c(DESeq_24h_original_up, DESeq_24h_original_down, edgeR_24h_original_up, edgeR_24h_original_down))
genes_DE_48h_any <- unique(c(DESeq_48h_original_up, DESeq_48h_original_down, edgeR_48h_original_up, edgeR_48h_original_down))

genes_DE_any_up <- unique(c(DESeq_12h_original_up, DESeq_24h_original_up, DESeq_48h_original_up, edgeR_12h_original_up, edgeR_24h_original_up, edgeR_48h_original_up))
genes_DE_any_down <- unique(c(DESeq_12h_original_down, DESeq_24h_original_down, DESeq_48h_original_down, DESeq_12h_original_down, DESeq_24h_original_down, DESeq_48h_original_down, edgeR_12h_original_down, edgeR_24h_original_down, edgeR_48h_original_down, edgeR_12h_original_down, edgeR_24h_original_down, edgeR_48h_original_down))

genes_DE_all_up <- unique(c(DESeq_12h_original_up,DESeq_24h_original_up,DESeq_48h_original_up))[unique(c(DESeq_12h_original_up,DESeq_24h_original_up,DESeq_48h_original_up)) %in% unique(c(edgeR_12h_original_up, edgeR_24h_original_up, edgeR_48h_original_up))]
genes_DE_all_down <- unique(c(DESeq_12h_original_down,DESeq_24h_original_down,DESeq_48h_original_down))[unique(c(DESeq_12h_original_down,DESeq_24h_original_down,DESeq_48h_original_down)) %in% unique(c(edgeR_12h_original_down, edgeR_24h_original_down, edgeR_48h_original_down))]

gene_name_any <- original_gtf$g_name[original_gtf$gene %in% genes_DE_any]
gene_name_any_up <- original_gtf$g_name[original_gtf$gene %in% genes_DE_any_up]
gene_name_any_down <- original_gtf$g_name[original_gtf$gene %in% genes_DE_any_down]
gene_name_all_up <- original_gtf$g_name[original_gtf$gene %in% genes_DE_all_up]
gene_name_all_down <- original_gtf$g_name[original_gtf$gene %in% genes_DE_all_down]

#compare with previous studies finding food-related differential expression
gene_up_food_up <- gene_name_any_up[gene_name_any_up %in% food_up_HT115$gene_name]
gene_down_food_down <- gene_name_any_down[gene_name_any_down %in% food_down_HT115$gene_name]


#create a file
myLatexOutput <- paste( dir[2], "results_latex.txt", sep="")
latexOut<-file(myLatexOutput, "w")

cat("% max nb of reads mapped  \n", file = latexOut)
cat(paste("\\newcommand{\\maxReads}{", max(stat_table$mapped) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% max percent of reads mapped  \n", file = latexOut)
cat(paste("\\newcommand{\\maxPercentReads}{", max(stat_table$percent_mapped) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% min nb of reads mapped  \n", file = latexOut)
cat(paste("\\newcommand{\\minReads}{", min(stat_table$mapped) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% min percent of reads mapped  \n", file = latexOut)
cat(paste("\\newcommand{\\minPercentReads}{", min(stat_table$percent_mapped) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% mean nb of reads mapped  \n", file = latexOut)
cat(paste("\\newcommand{\\meanReads}{", mean(stat_table$mapped) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% mean percent of reads mapped  \n", file = latexOut)
cat(paste("\\newcommand{\\meanPercentReads}{", mean(stat_table$percent_mapped) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% max percent of reads trimmed  \n", file = latexOut)
cat(paste("\\newcommand{\\maxPercentTrim}{", max(stat_table$percent_trimmed) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% min percent of reads trimmed  \n", file = latexOut)
cat(paste("\\newcommand{\\minPercentTrim}{", min(stat_table$percent_trimmed) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% mean percent of reads trimmed  \n", file = latexOut)
cat(paste("\\newcommand{\\meanPercentTrim}{", mean(stat_table$percent_trimmed) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)


cat("% max number of gene expressed  \n", file = latexOut)
cat(paste("\\newcommand{\\maxGenes}{", max(vect_gene) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% sample with max number of gene expressed  \n", file = latexOut)
cat(paste("\\newcommand{\\sampleMaxGene}{", gsub("_"," ",colnames(table_TPM_original)[which.max(vect_gene)]) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% min number of gene expressed  \n", file = latexOut)
cat(paste("\\newcommand{\\minGenes}{", min(vect_gene) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% sample with min number of gene expressed  \n", file = latexOut)
cat(paste("\\newcommand{\\sampleMinGene}{", gsub("_"," ",colnames(table_TPM_original)[which.min(vect_gene)]) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% mean number of gene expressed  \n", file = latexOut)
cat(paste("\\newcommand{\\meanGenes}{", mean(vect_gene) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% nb of genes expressed across all samples  \n", file = latexOut)
cat(paste("\\newcommand{\\genesAllSamples}{", length(table_TPM_original[rowSums(table_TPM_original==0)<1,1]) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% nb of genes expressed in at least one samples  \n", file = latexOut)
cat(paste("\\newcommand{\\genesOneSample}{", length(table_TPM_original[rowSums(table_TPM_original==0)<length(table_TPM_original[1,]),1]) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)


cat("% number of genes DE in any condition  \n", file = latexOut)
cat(paste("\\newcommand{\\allDEgenes}{", length(gene_name_any) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% list of genes DE in any condition  \n", file = latexOut)
cat(paste("\\newcommand{\\allDEgeneslist}{", paste(gene_name_any, collapse=", ") ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)

cat("% number of genes DE in any condition up in HT115 \n", file = latexOut)
cat(paste("\\newcommand{\\allDEgenesUp}{", length(gene_name_any_up) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% list of genes DE in any condition up in HT115  \n", file = latexOut)
cat(paste("\\newcommand{\\allDEgenesUplist}{", paste(gene_name_any_up, collapse=", ") ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)

cat("% number of genes DE at 12h \n", file = latexOut)
cat(paste("\\newcommand{\\allDEgenesTwelve}{", genes_DE_12h_any ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% number of genes DE at 24h \n", file = latexOut)
cat(paste("\\newcommand{\\allDEgenesTwentyfour}{", genes_DE_24h_any ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% number of genes DE at 48h \n", file = latexOut)
cat(paste("\\newcommand{\\allDEgenesFourtyeight}{", genes_DE_48h_any ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)

cat("% number of genes DE in any condition down in HT115 \n", file = latexOut)
cat(paste("\\newcommand{\\allDEgenesDown}{", length(gene_name_any_down) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% list of genes DE in any condition down in HT115  \n", file = latexOut)
cat(paste("\\newcommand{\\allDEgenesDownlist}{", paste(gene_name_any_down, collapse=", ") ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)

cat("% number of genes DE in any condition up in HT115 previously found involved in regeneration\n", file = latexOut)
cat(paste("\\newcommand{\\geneFoundRegeneration}{", length(genes_DE_any_up[genes_DE_any_up %in% list_all_regeneration]) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)

cat("% number of genes DE in any condition up in HT115 previously found up in HT115 in prev studies\n", file = latexOut)
cat(paste("\\newcommand{\\geneUpFoodUp}{", length(gene_up_food_up) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% list of genes DE in any condition up in HT115 previously found up in HT115 in prev studies\n", file = latexOut)
cat(paste("\\newcommand{\\listgeneUpFoodUp}{", paste(gene_up_food_up, collapse=", ") ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)

cat("% number of genes DE in any condition down in HT115 previously found down in HT115 in prev studies\n", file = latexOut)
cat(paste("\\newcommand{\\geneDownFoodDown}{", length(gene_down_food_down) ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)
cat("% list of genes DE in any condition down in HT115 previously found down in HT115 in prev studies\n", file = latexOut)
cat(paste("\\newcommand{\\listgeneDownFoodDown}{", paste(gene_down_food_down, collapse=", ") ,"}", sep = ""), file = latexOut)
cat("\n\n", file = latexOut)


close(latexOut)

#save data in a file final
f <- paste ( dir[2] , 'savedData.RData' , sep='' )
save.image(f)
                   
