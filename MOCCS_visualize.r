# The MIT License (MIT)
# 
# Copyright (c) <2015> <Haruka Ozaki>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# version 1.5

options(stringsAsFactors = F)

args=(commandArgs(TRUE))

if(length(args)==0){
        print("No arguments supplied.")
}else{
        for(i in 1:length(args)){
                eval(parse(text=args[[i]]))
        }
}

print("Visualizing MOCCS result...")
sprintf("file: %s", file)
sprintf("file2: %s", file2)
sprintf("label: %s", label)
sprintf("len: %s", len)
sprintf("k: %s", k)


# Universal colors
cols <- c("#000000", "#ff2800","#faf500","#35a16b","#0041ff","#66ccff","#ff99a0","#ff9900","#9a0079","#663300")


crf <- read.table(gzfile(file), header=T, stringsAsFactors=F, check.names=F, row.names=1)
crf <- as.data.frame(t(apply(crf, 1, function(x){x/sum(x)})))
d.auc <- read.table(file2, header=T, stringsAsFactors=F, row.names=1)
auc <- d.auc$auc
names(auc) <- rownames(d.auc)
auc <- auc[ order(auc, decreasing=T) ]

xlab <- "Distance from peak center [bp]"
ylab <- "Cumulative relative freqency"
end <- as.numeric(colnames(crf)[ncol(crf)])
print(sprintf("x-axis: [0, %d]", end))

# make roc-like curve plot of auc for top 10 k-mers
auc <- auc[1:10]
crf2 <- crf[names(auc),]

pdfname <- paste0(label, "_", k, "-mer_auc_plot", ".pdf")
sprintf("outfile: %s", pdfname)

pdf(pdfname)

matplot(t(crf2), type="n", col="gray", lwd=0.5, main=label, lty=1,
	sub=paste0(k, "-mers"), xlab=xlab, ylab=ylab, axes=F)
box()
# axis(1, 1:ncol(crf2), 0:end)
axis(1, seq(1, ncol(crf2), by=50), seq(0, end, by = 50)) 
axis(2)

# i.max <- ifelse(length(auc)<10, length(auc), 10)
i.max = min(length(auc), 10)
for(i in i.max:1){
	matpoints(t(crf2[grep(names(auc)[i],rownames(crf2)),]), type="l", col=cols[i], lwd=1, lty=1)
}
legend("bottomright", legend=names(auc[1:i.max]), col=cols[1:i.max], lwd=1)


# auc.th <- d.auc[d.auc$count > count.th, 1]
# if(length(auc.th) == 0){
# 	dev.off()
# 	print("No k-mers exceeds count threshold.")
# 	quit("no")
# }

# names(auc.th) <- rownames(d.auc[d.auc$count > count.th, ])

# auc.th <- auc.th[order(auc.th, decreasing=T)]
# crf <- crf[names(auc.th), ]

# auc.th <- auc.th[1:min(10, length(auc.th))]


# if(length(auc.th) > 0){
# 	matplot(t(crf[names(auc.th),]), type="n", col="gray", lwd=0.5, main=label, axes=F, lty=1,
# 		sub=paste0(k, "-mer whose count>", count.th, ": ", length(auc.th)), xlab=xlab, ylab=ylab)
# 	box()
# 	axis(1, 1:ncol(crf), 0:end)
# 	axis(2)
# 	# i.max <- ifelse(length(auc.th)<10, length(auc.th), 10)
# 	i.max = min(length(auc.th), 10)
# 	for(i in i.max:1){
# 		matpoints(t(crf[names(auc)[i],]), type="l", col=cols[i], lwd=1, lty=1,)
# 	}
# 	legend("bottomright", legend=names(auc.th[1:i.max]), col=cols[i:i.max], lwd=1)
# }

dev.off()






# # make roc-like curve plot of auc
# pdfname <- paste0(label, "_", k, "-mer_auc_plot", ".pdf")
# pdf(pdfname)

# matplot(t(crf), type="l", col="gray", lwd=0.5, main=label, lty=1,
# 	sub=paste0("all", w, "-mer"), xlab=xlab, ylab=ylab, axes=F)
# box()
# axis(1, 1:ncol(crf), 0:end)
# axis(2)

# i.max <- ifelse(length(auc)<10, length(auc), 10)
# for(i in i.max:1){
# 	matpoints(t(crf[grep(names(auc)[i],rownames(crf)),]), type="l", col=cols[i], lwd=1, lty=1)
# }
# legend("bottomright", legend=names(auc[1:i.max]), col=cols[1:i.max], lwd=1)

# auc.th <- d.auc[d.auc$count > count.th, 1]
# names(auc.th) <- rownames(d.auc[d.auc$count > count.th, ])
# auc.th <- auc.th[order(auc.th, decreasing=T)]
# if(length(auc.th) > 0){
# 	matplot(t(crf[names(auc.th),]), type="l",col="gray", lwd=0.5, main=label, axes=F, lty=1,
# 		sub=paste0(k, "-mer whose count>", count.th, ": ", length(auc.th)), xlab=xlab, ylab=ylab)
# 	box()
# 	axis(1, 1:ncol(crf), 0:end)
# 	axis(2)
# 	i.max <- ifelse(length(auc.th)<10, length(auc.th), 10)
# 	for(i in i.max:1){
# 		matpoints((crf[grep(names(auc.th[i]),rownames(crf)),]), type="l", col=cols[i], lwd=1, lty=1,)
# 	}
# 	legend("bottomright", legend=names(auc.th[1:i.max]), col=cols[i:i.max], lwd=1)
# }

# dev.off()


# pdfname <- paste0(label, "_", k, "-mer_stat", ".pdf")
# pdf(pdfname)

# # make histgram of 6-mer count
# bin <- 10
# hist(d.auc$count, breaks=i, sub=paste0("all", k, "-mer count ", label, " breaks=", bin), xlab="Count", ylab="# of k-mer")
# bin <- 100
# hist(d.auc$auc, breaks=i, sub=paste0("all", k, "-mer count ", label, " breaks=", bin), xlab="AUC", ylab="# of k-mer")

# # # make plot of auc and count
# # plot(d.auc$auc, d.auc$count, type="p", col = "#00880040", pch=19,
# # 	xlab="AUC", ylab="Count", sub=paste0("all", k, "-mer count ", label))

# dev.off() 

print("Visualized MOCCS result")

sessionInfo()
