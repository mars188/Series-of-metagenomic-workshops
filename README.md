# Series-of-metagenomic-workshops
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
