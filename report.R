#!/usr/bin/env Rscript

library('tseries')

outputpref <- Sys.getenv("outputpref")
workdir <- Sys.getenv("workdir")

kmer_info <- read.csv(file=paste(outputpref, ".stat.tab",sep=""),sep="\t", header=TRUE,colClasses=c("numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))

limi_loc<-which.max(kmer_info[,2])
siz<-length(kmer_info[,2])
df <- data.frame()
num=1
for (i in 1:limi_loc){
	for (j in limi_loc:siz){
		if (j-i>3){
			adf_res <- adf.test(kmer_info[i:j,2], k=0)
			adf_res2 <- adf.test(kmer_info[j:i,2], k=0)
			df[num,1]=kmer_info[i,1]
			df[num,2]=kmer_info[j,1]
			df[num,3]=adf_res$p.value
			df[num,4]=kmer_info[j,1]-kmer_info[i,1]
			df[num+1,1]=kmer_info[j,1]
			df[num+1,2]=kmer_info[i,1]
			df[num+1,3]=adf_res2$p.value
			df[num+1,4]=kmer_info[j,1]-kmer_info[i,1]
			num<-num+2			
		}
	}
}
sorted_indices <- order(df[,3], -df[,4])
sorted_df <- df[sorted_indices, ]

summaryFile <- paste(outputpref, ".report",sep="")
cat(paste("Limiting Value: ", max(kmer_info[,2]), sep=""), file=summaryFile, sep="\n", append=TRUE)
cat(paste("Limiting K: ", kmer_info[limi_loc,1], sep=""), file=summaryFile, sep="\n", append=TRUE)
cat(paste("Max convergent segment ranging from: ", sorted_df[1,1], "-mer", " to ", sorted_df[1,2], "-mer;", "Dickey-Fuller test: P-value = ", sorted_df[1,3], sep=""), file=summaryFile, sep="\n", append=TRUE)
cat(paste("\n\n\n","#################################################################################", sep=""), file=summaryFile, sep="\n", append=TRUE)
cat("If you use our tools in your work, we kindly ask that you also cite the following remarkable tools that have been integrated into LVgs:", file=summaryFile, sep="\n", append=TRUE)
cat("    1) FastK: https://github.com/thegenemyers/FASTK", file=summaryFile, sep="\n", append=TRUE)
cat("    2) GenomeScope2: Ranallo-Benavidez, T. R., Jaron, K. S. & Schatz, M. C. GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. Nat. Commun. 11, doi:10.1038/s41467-020-14998-3 (2020)", file=summaryFile, sep="\n", append=TRUE)
cat("    3) HiFiasm: Cheng, H., Asri, M., Lucas, J., Koren, S. & Li, H. Scalable telomere-to-telomere assembly for diploid and polyploid genomes with double graph. Nat Methods 21, 967-970, doi:10.1038/s41592-024-02269-8 (2024).", file=summaryFile, sep="\n", append=TRUE)
cat("    4) tseries: Trapletti, A. & Hornik, K. tseries: Time Series Analysis and Computational Finance. (R package version 0.10-58: https://CRAN.R-project.org/package=tseries, 2024)", file=summaryFile, sep="\n", append=TRUE)
cat("    5) KMC: Marek Kokot, Maciej Długosz, Sebastian Deorowicz, KMC 3: counting and manipulating k-mer statistics, Bioinformatics, Volume 33, Issue 17, 01 September 2017, Pages 2759–2761, https://doi.org/10.1093/bioinformatics/btx304 ", file=summaryFile, sep="\n", append=TRUE)
