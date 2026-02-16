# Series-of-metagenomic-workshops
```
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

