# The MIT License (MIT)
# 
# Copyright (c) <2016> <Haruka Ozaki>
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

# version 1.7

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
crf <- as.data.frame(t(apply(crf, 1, function(x){x/x[length(x)]})))
d.auc <- read.table(file2, header=T, stringsAsFactors=F, row.names=1)
auc <- d.auc$auc
names(auc) <- rownames(d.auc)
auc <- auc[ order(auc, decreasing=T) ]

xlab <- "Distance from peak center [bp]"
ylab <- "Cumulative relative freqency"
end <- as.numeric(colnames(crf)[ncol(crf)])
print(sprintf("x-axis: [0, %d]", end))

# make roc-like curve plot of auc for top 10 k-mers
# auc <- auc[1:10]
crf2 <- crf[names(auc),]

pdfname <- paste0(label, "_", k, "-mer_auc_plot", ".pdf")
sprintf("outfile: %s", pdfname)

pdf(pdfname)

matplot(t(crf2), type="n", col="gray", lwd=0.5, main=label, lty=1,
	sub=paste0("Top ", k, "-mers based on AUC"), xlab=xlab, ylab=ylab, axes=F)
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



dev.off()

#################################
# Top 10 k-mers based on MOCCS2score

##################################################################
# Top 10 k-mers based on MOCCS2score
MOCCS2score <- d.auc$MOCCS2score
names(MOCCS2score) <- rownames(d.auc)
MOCCS2score <- MOCCS2score[ order(MOCCS2score, decreasing=T) ]

xlab <- "Distance from peak center [bp]"
ylab <- "Cumulative relative freqency"
end <- as.numeric(colnames(crf)[ncol(crf)])
print(sprintf("x-axis: [0, %d]", end))

# make roc-like curve plot of MOCCS2score for top 10 k-mers
# MOCCS2score <- MOCCS2score[1:10]
crf2 <- crf[names(MOCCS2score),]

pdfname <- paste0(label, "_", k, "-mer_MOCCS2score_plot", ".pdf")
sprintf("outfile: %s", pdfname)

pdf(pdfname)

matplot(t(crf2), type="n", col="gray", lwd=0.5, main=label, lty=1,
	sub=paste0("Top ", k, "-mers based on MOCCS2score"), xlab=xlab, ylab=ylab, axes=F)
box()
# axis(1, 1:ncol(crf2), 0:end)
axis(1, seq(1, ncol(crf2), by=50), seq(0, end, by = 50)) 
axis(2)

# i.max <- ifelse(length(MOCCS2score)<10, length(MOCCS2score), 10)
i.max = min(length(MOCCS2score), 10)
for(i in i.max:1){
	matpoints(t(crf2[grep(names(MOCCS2score)[i],rownames(crf2)),]), type="l", col=cols[i], lwd=1, lty=1)
}
legend("bottomright", legend=names(MOCCS2score[1:i.max]), col=cols[1:i.max], lwd=1)



dev.off()

#################################



print("Visualized MOCCS result")

sessionInfo()
