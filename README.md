# Series of Metagenomics Workshop
During this course, participants will learn about the concepts and steps involved in analyzing the metagenomic data obtained from short read high throughput sequencing.

Participants are expected to have some basic knowledge of metagenome/microbiome as well as some biological knowledge on the subject matter.

Although the material and the methods are designed to cater for the NYU Abu Dhabi High Performance Computing environment, it is possible to run this workshop on any other system provided that the neccessary software is installed, and the data is uploaded to that environment. It will also be neccessary to change the input and output directories/files to accomodate such an environment.

### Connecting to the HPC using a MAC machine ###
* Open the "Terminal" app and type ssh NetID@jubail.abudhabi.nyu.edu
* Enter your NYU password when prompted and hit enter
* Once logged in, navigate to your personal "SCRATCH" directory or the any other subdirecotry where you want to run this workshop

### Connecting to the HPC using a Windows machine ###
* Open the "Putty" app, and fill out the fields as follows Host name jubail.abudhabi.nyu.edu, Port=22, and then click on "Open".
* Enter your NetId, and your password when prompted.
* Once logged in, navigate to your personal "SCRATCH" directory or the any other subdirecotry where you want to run this workshop.

# MODULE 1 #
### The dataset ###
The shotgun metagenomic sequencing data that we will use for this hands-on training was obtained from a recently published study in the Microbiome by **Jeilu et al. 2025** [https://doi.org/10.1186/s40168-025-02276-7]. The authors deposited the data in the NCBI under the BioProject accession number PRJNA1228129.

Briefly, the study was undertaken for metagenomic profiling of airborne microbial communities from aircraft filters and face masks. This study compared the following four groups. 

- Aircraft filter (aircraft cabin filter was removed from a commercial airplane after 8039 h of usage)
- Travel (Face masks worn by air travelers)
- Hospital (Face masks worn by healthcare professionals in a hospital)
- Control (Unworn)

For the sake of time and resources, we will be using the following two groups in this workshop to analyze the microbial profiling using this dataset. 

- Aircraft filter
- Control

We will be using the NYUAD High Performance Computing (HPC) cluster for this workshop, however, you can certainly run all of the analysis on any stand alone machine (server, personal laptop/Desktop etc.) provided that you have pre-installed the necessay software packages.

```
mkdir -p /scratch/$USER/metagenomics_workshop/
mkdir -p /scratch/$USER/metagenomics_workshop/raw_fastqs/
cd /scratch/$USER/metagenomics_workshop/
cp /scratch/gencore/ma5877/metagenomic_course/raw_fastqs/* raw_fastqs/

```
Now, we will run the QC/QT of the reads. 
Basically, we will check the quality of reads first then perform quality trimming to remove low quality reads, adapter sequences and finally also decontaminate the reads meaning that we will remove host (human) sequences. At the end, this step will produce clean, high quality reads free from host contamination for the downstream analysis. There are several tools available for this purpose including metawrap, bbtools etc but we will be using **KneadData** that is part of the Biobakery suite. 

So, we need to load **KneadData** module as below. 

```
module purge
module load all
module load gencore/3
module load kneaddata/0.12.4
```
Now, if we run ``` kneaddata --help ``` command, it pulls a detailed menu of this program and we can see all the available options that we can use. Please run this help command and try to find which option should go in the provided blanks.   

```
kneaddata \
-i1 {$self->root_in_dir}/{$sample}_1.fastq \
-i2 {$self->root_in_dir}/{$sample}_2.fastq \
--reference-db /scratch/Reference_Genomes/Public/Metagenomic/kneaddata_hg39/ \
--reference-db /scratch/Reference_Genomes/Public/Metagenomic/kneaddata_phix/ \
--trimmomatic-options="SLIDINGWINDOW:4:20" --trimmomatic-options="MINLEN:50" \
--run-trim-repetitive \
--output {$self->kneaddata_dir}/ \
--output-prefix {$sample} \
--remove-intermediate-output \
--fastqc fastqc \
--threads 14 && \
fastqc {$self->kneaddata_dir}/{$sample}_paired_1.fastq -t 14 -o {$self->kneaddata_dir}/fastqc/ && \
fastqc {$self->kneaddata_dir}/{$sample}_paired_2.fastq -t 14 -o {$self->kneaddata_dir}/fastqc
```
The results of this run has been shared with you. You should compare the results that you obtain be running the QC/QT with the one ready-made provided to you. Feel free to ask question or discuss if you obtained very different results between these two sets. 

# MODULE 2 #

module load all gencore/3
module load R/4.3.1
module load biopython/1.85

mkdir -p diversity_plots/
cd diversity_plots/
cp ../data/analysis/*/kraken2/*_profile.txt .
cp ../data/analysis/*/bracken/*_S_Sh.txt .
python /scratch/Reference_Genomes/Public/Metagenomic/KrakenTools/DiversityTools/beta_diversity.py -i *_profile.txt --type kreport2 -l S > merged_beta_div.txt
(echo -e "Sample\tDiversity"; awk 'FNR==1{f=FILENAME; sub(/^.*\//,"",f); sub("_S_Sh.txt","",f)} /Shannon/{print f"\t"$NF}' *_S_Sh.txt) > merged_alpha_div.txt
rm *_profile.txt
rm *_S_Sh.txt

Rscript -e 'df<-read.table("merged_alpha_div.txt",header=TRUE,stringsAsFactors=FALSE); df$Group<-ifelse(grepl("^AF",df$Sample),"Filter",ifelse(grepl("^C",df$Sample),"Control","Other")); write.table(df,"merged_alpha_div_with_group.txt",sep="\t",quote=FALSE,row.names=FALSE)'

Rscript -e 'df<-read.table("merged_alpha_div_with_group.txt",header=TRUE,stringsAsFactors=FALSE);
ymax<-max(df$Diversity)*1.15;
png("diversity_boxplot_stats.png",width=800,height=600,res=120);
boxplot(Diversity~Group,data=df,col=c("steelblue","orange"),ylab="Shannon Diversity",main="Shannon Diversity by Group",ylim=c(min(df$Diversity),ymax));
stripchart(Diversity~Group,data=df,vertical=TRUE,method="jitter",pch=16,col="black",add=TRUE);
w<-wilcox.test(Diversity~Group,data=df);
p_txt<-paste0("Wilcoxon p = ",signif(w$p.value,3));
text(x=1.5,y=ymax,labels=p_txt);
dev.off()'

Rscript -e 'args<-commandArgs(trailingOnly=TRUE); infile<-args[1]; outfile<-args[2]; lines<-readLines(infile); samp_lines<-grep("^#", lines, value=TRUE); samps<-sub("^#\\d+\\s+","",samp_lines); samps<-sub("_profile.txt.*","",samps); start<-grep("^x\\s", lines); mat<-read.table(text=lines[start:length(lines)], header=TRUE, stringsAsFactors=FALSE, check.names=FALSE); mat<-mat[ , -1]; mat[mat=="x.xxx"]<-NA; mat<-as.matrix(sapply(mat, as.numeric)); rownames(mat)<-samps; colnames(mat)<-samps; mat[lower.tri(mat)]<-t(mat)[lower.tri(mat)]; diag(mat)<-0; write.table(mat, file=outfile, sep="\t", quote=FALSE, col.names=NA)' merged_beta_div.txt cleaned_beta_matrix.txt

Rscript -e 'library(vegan); mat<-as.matrix(read.table("cleaned_beta_matrix.txt", header=TRUE, check.names=FALSE, row.names=1)); d<-as.dist(mat); p<-cmdscale(d, k=2, eig=TRUE); scores<-as.data.frame(p$points); scores$Sample<-rownames(scores); scores$Group<-ifelse(grepl("^AF", scores$Sample), "AirFilter", ifelse(grepl("^C", scores$Sample), "Control", "Other")); ad<-adonis2(d ~ Group, data=scores); pval<-ad$`Pr(>F)`[1]; var_expl<-round(p$eig/sum(p$eig)*100,1); png("beta_pcoa.png", width=800, height=600, res=120); par(mar=c(5,5,4,2)); cols<-c(AirFilter="steelblue", Control="orange", Other="gray"); plot(scores$V1, scores$V2, xlab=paste0("PCoA1 (", var_expl[1], "%)"), ylab=paste0("PCoA2 (", var_expl[2], "%)"), pch=16, col=cols[scores$Group], main="PCoA of beta diversity"); legend("topright", legend=unique(scores$Group), col=cols[unique(scores$Group)], pch=16, bty="n"); text(x=min(scores$V1), y=max(scores$V2), labels=paste0("PERMANOVA p = ", signif(pval,3)), adj=c(0,1)); dev.off()'
```

