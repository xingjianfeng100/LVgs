#!/usr/bin/env Rscript

outputpref <- Sys.getenv("outputpref")
workdir <- Sys.getenv("workdir")
DFtest <- Sys.getenv("DFtest")

kmer_info <- read.csv(file=paste(outputpref, ".stat.tab",sep=""),sep="\t", header=TRUE,colClasses=c("numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))

pdf(file=paste(outputpref, "_figure/", "GS_linear_plot.pdf",sep=""))
par(mar = c(5.0,5.0,4.0,5.0))
plot(kmer_info[,1],kmer_info[,2]/1000000, type = "b", pch = 23,lty = 2,
     col = "blue", col.axis = "orange1",
     col.lab = rgb(198, 156, 109, maxColorValue = 255), 
     fg = "pink", bg = "green",
     font.lab = 1, 
     cex.lab = 1.5, 
     xlab="K-mer length (Bp)", ylab="Genome size (Mb)",
     xlim = c(0,max(kmer_info[,1])*1.01), ylim = c(0,max(kmer_info[,2])*1.01/1000000))
abline(h = max(kmer_info[,2])*1.000000001/1000000,lwd = 1.5,lty = 2,col = "rosybrown")
axis(4,at = max(kmer_info[,2])*1.000000001/1000000,labels = paste(format(round(max(kmer_info[,2])*1.000000001/1000000),big.mark=","), "Mb", "Limiting Value", sep=" "),col.axis = "blue",las = 2, cex.axis = 0.5,tck = 0,hadj = 0.1)
mtext(DFtest, side = 3, line = 0, adj = 0, cex = 0.8, font = 2)
dev.off()

pdf(file=paste(outputpref, "_figure/", "Repeat_linear_plot.pdf",sep=""))
par(mar = c(5.0,5.0,4.0,5.0))
plot(kmer_info[,1],kmer_info[,3]/1000000, type = "b", pch = 23,lty = 2,
     col = "blue", col.axis = "orange1",
     col.lab = rgb(198, 156, 109, maxColorValue = 255), 
     fg = "pink", bg = "green",
     font.lab = 1,
     cex.lab = 1.5,
     xlab="K-mer length (Bp)", ylab="Repeat lenth (Mb)",
     xlim = c(0,max(kmer_info[,1])*1.01), ylim = c(0,max(kmer_info[,3])*1.01/1000000))
dev.off()


pdf(file=paste(outputpref, "_figure/", "Unique_linear_plot.pdf",sep=""))
par(mar = c(5.0,5.0,4.0,5.0))
plot(kmer_info[,1],kmer_info[,4]/1000000, type = "b", pch = 23,lty = 2,
     col = "blue", col.axis = "orange1",
     col.lab = rgb(198, 156, 109, maxColorValue = 255), 
     fg = "pink", bg = "green",
     font.lab = 1,
     cex.lab = 1.5,
     xlab="K-mer length (Bp)", ylab="Unique lenth (Mb)",
     xlim = c(0,max(kmer_info[,1])*1.01), ylim = c(0,max(kmer_info[,4])*1.01/1000000))
dev.off()


pdf(file=paste(outputpref, "_figure/", "Model_fit_plot.pdf",sep=""))
par(mar = c(5.0,5.0,4.0,5.0))
plot(kmer_info[,1],kmer_info[,5], type = "b", pch = 23,lty = 2,
     col = "blue", col.axis = "orange1",
     col.lab = rgb(198, 156, 109, maxColorValue = 255),
     fg = "pink", bg = "green",
     font.lab = 1,
     cex.lab = 1.5,
     xlab="K-mer length (Bp)", ylab="Model fit (%)",
     xlim = c(0,max(kmer_info[,1])*1.01), ylim = c(0,max(kmer_info[,5])*1.01))
dev.off()


pdf(file=paste(outputpref, "_figure/", "Read_Error_plot.pdf",sep=""))
par(mar = c(5.0,5.0,4.0,5.0))
plot(kmer_info[,1],kmer_info[,6], type = "b", pch = 23,lty = 2,
     col = "blue", col.axis = "orange1",
     col.lab = rgb(198, 156, 109, maxColorValue = 255),
     fg = "pink", bg = "green",
     font.lab = 1,
     cex.lab = 1.5,
     xlab="K-mer length (Bp)", ylab="Read error evaluation (%)",
     xlim = c(0,max(kmer_info[,1])*1.01), ylim = c(0,max(kmer_info[,6])*1.01))
dev.off()

pdf(file=paste(outputpref, "_figure/", "kmer_type_plot.pdf",sep=""))
par(mar = c(5.0,5.0,4.0,5.0))
plot(kmer_info[,1],kmer_info[,7], type = "b", pch = 23,lty = 2,
     col = "blue", col.axis = "orange1",
     col.lab = rgb(198, 156, 109, maxColorValue = 255),
     fg = "pink", bg = "green",
     font.lab = 1,
     cex.lab = 1.5,
     xlab="K-mer length (Bp)", ylab="K-mer type number",
     xlim = c(0,max(kmer_info[,1])*1.01), ylim = c(min(kmer_info[,7])*0.9,max(kmer_info[,7])*1.01))
dev.off()

pdf(file=paste(outputpref, "_figure/", "kmer_total_plot.pdf",sep=""))
par(mar = c(5.0,5.0,4.0,5.0))
plot(kmer_info[,1],kmer_info[,8], type = "b", pch = 23,lty = 2,
     col = "blue", col.axis = "orange1",
     col.lab = rgb(198, 156, 109, maxColorValue = 255),
     fg = "pink", bg = "green",
     font.lab = 1,
     cex.lab = 1.5,
     xlab="K-mer length (Bp)", ylab="K-mer total number",
     xlim = c(0,max(kmer_info[,1])*1.01), ylim = c(min(kmer_info[,8])*0.9,max(kmer_info[,8])*1.01))
dev.off()

pdf(file=paste(outputpref, "_figure/", "kcov_plot.pdf",sep=""))
par(mar = c(5.0,5.0,4.0,5.0))
plot(kmer_info[,1],kmer_info[,9], type = "b", pch = 23,lty = 2,
     col = "blue", col.axis = "orange1",
     col.lab = rgb(198, 156, 109, maxColorValue = 255),
     fg = "pink", bg = "green",
     font.lab = 1,
     cex.lab = 1.5,
     xlab="K-mer length (Bp)", ylab="Kcov (x)",
     xlim = c(0,max(kmer_info[,1])*1.01), ylim = c(0,max(kmer_info[,9])*1.01))
dev.off()
