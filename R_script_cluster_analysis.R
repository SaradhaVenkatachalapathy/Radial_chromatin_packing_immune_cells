library(gplots)
library(RColorBrewer)
library(dendextend)
library(plot3D)  
library(factoextra)  
library(cluster)  


dird<-"/Users/saradhavenkatachalapathy/Desktop/immunecells/new_run/"
dir.create(dird)
setwd(dird)

#read in the data
protien_ratios3d<- read.csv("~/Desktop/immunecells/comb.csv")
protien_ratios3d<-protien_ratios3d[,-13]
rownames(protien_ratios3d)<-protien_ratios3d$Label

#Set color libraries
rbPal <- colorRampPalette(c('red','white','black'))#colors for the protein ratio
Colors=gray.colors(15, start = 0.9, end = 0.3)# colors for the heatmap


#Color code the foldchange
protien_ratios3d_mat<-na.omit(protien_ratios3d[,c(13,3:12)])
protien_ratios3d_mat[,1]<-log2(protien_ratios3d_mat[,1])
protien_ratios3d_mat<-protien_ratios3d_mat[order(protien_ratios3d_mat[,1]),]
protien_ratios3d_mat$cols <- rbPal(10)[as.numeric(cut(protien_ratios3d_mat[,1],breaks = 10))]

#Obtain the clusters
distmat<-as.dist(1-cor(t(protien_ratios3d_mat[,2:11]), method="spearman"))
rowclusters = hclust(distmat, method="complete")
mycl <- cutree(rowclusters, h=max(rowclusters$height/1.5))

{
  cluster_fun <- function(x, k) list(cluster = cutree(hclust(as.dist(1-cor(t(x), method="spearman")), method="complete"), k = k))
  
  png(filename="Cluster_silhouette.png", units="in",width=3, height=3 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2, font=2)
  fviz_nbclust(as.matrix(protien_ratios3d_mat[,2:11]), FUN = cluster_fun, diss=distmat,method = "silhouette")
  dev.off()
  
  png(filename="Cluster_wss.png", units="in",width=3, height=3 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2, font=2)
  fviz_nbclust(as.matrix(protien_ratios3d_mat[,2:11]), FUN = cluster_fun, method = "wss")
  dev.off()
  
  
  png(filename="Cluster_gapstatistic.png", units="in",width=3, height=3 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2, font=2)
  gap_stat <- clusGap(as.matrix(protien_ratios3d_mat[,2:11]), FUN  = cluster_fun, K.max = 10, B = 50)
  fviz_gap_stat(gap_stat,maxSE = list(method = "globalmax",maxSE=0.5))
  dev.off()
  
  
  }

clusterCols <- c("blue","forestgreen")
dend1 <- as.dendrogram(rowclusters)
dend1 <- color_branches(dend1, k = 2, col = clusterCols)
col_labels <- get_leaves_branches_col(dend1)
col_labels <- col_labels[order(order.dendrogram(dend1))]



png(filename="Heatmap_with_protein_ratios_as_rowside_colors.png", units="in",width=4, height=4 , pointsize=7, res=1200)
par(font.axis = 2,font.lab=2, font=2)
data_heat<-heatmap.2(as.matrix(protien_ratios3d_mat[,2:11]),scale="none",Rowv = dend1, Colv=FALSE, dendrogram = "row", 
                     density.info = "none",trace="none",col=Colors,labRow = NA, key.xlab="Intensity Fraction" ,keysize =1.5,key.title=" ",
                     RowSideColors = protien_ratios3d_mat$cols,srtCol = 45,margins=c(4,2),
                     labCol =c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))
colkey (col=protien_ratios3d_mat$cols, clim = range(protien_ratios3d_mat$Cor.RPL),side=3,add=T,length=0.2,cex.axis=0.8,clab = "log2(Coro1A/RPL10A)", padj=1,dist=-0.1, cex.clab =0.7)
dev.off()

png(filename="Heatmap.png", units="in",width=4, height=4 , pointsize=7, res=1200)
par(font.axis = 2,font.lab=2, font=2)
data_heat<-heatmap.2(as.matrix(protien_ratios3d_mat[,2:11]),scale="none",Rowv = dend1, Colv=FALSE, dendrogram = "row", 
                     density.info = "none",trace="none",col=Colors,labRow = NA, key.xlab="Intensity Fraction" ,keysize =1.5,key.title=" ",
                     srtCol = 45,margins=c(4,2),
                     labCol =c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))
dev.off()



#Obtain the clusters and merge the protien ratio values
data_heat_cl<-as.data.frame((t(data_heat$carpet)))
data_heat_cl$Label<-row.names(data_heat_cl)
mycl<-as.data.frame(mycl)
mycl$Label<-row.names(mycl)
comb<-merge(data_heat_cl,mycl,by="Label")
comb<-merge(comb,protien_ratios3d[,c(2,13)],by="Label")
comb$Cor.RPL<-log2(comb$Cor.RPL)

central<-subset(comb, comb$mycl==1)
perpheral<-subset(comb, comb$mycl==2)
centralmean<-colMeans(x=central[,2:11], na.rm = TRUE)
perpheralmean<-colMeans(x=perpheral[,2:11], na.rm = TRUE)
centralsd<-apply(central[,2:11],2,sd)
perpheralsd<-apply(perpheral[,2:11],2,sd)
central_CI<-1.96*(centralsd/nrow(central))
perpheral_CI<-1.96*(perpheralsd/nrow(perpheral))


#plot the clusterwise intesity change
png(filename="clusterwise_intensity_change.png", units="in",width=2, height=2 , pointsize=7, res=1200)
par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
plot(centralmean, col=clusterCols[2], type="o",pch=19,las=1,xaxt = "n",xlab="Radial Distance",main="",ylab="DNA Intensity", ylim=c(0,1))
polygon(c(1:10,10:1),c((centralmean+centralsd),rev(centralmean-centralsd)), col=alpha.col(col = clusterCols[2], alpha = 0.2), border= NA)
points(perpheralmean,col=clusterCols[1],type="o",pch=18)
polygon(c(1:10,10:1),c((perpheralmean+perpheralsd),rev(perpheralmean-perpheralsd)), col=alpha.col(col = clusterCols[1], alpha = 0.2), border= NA)
axis(1, at=1:10, labels=c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"),las=2,cex.axis=0.7)
dev.off()

#plot the log foldchange of Cor/RPL in the two clusters
png(filename="clusterwise_protein_change.png", units="in",width=2, height=2 , pointsize=7, res=1200)
par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
d1<-density(central$Cor.RPL)
d2<-density(perpheral$Cor.RPL)
plot(d1, col=clusterCols[2], type="l",las=1,xlab="Log2 (Coro1A/RPL10A)",main="",ylab="Probability Density", ylim=c(0,2.5), lwd=2)
lines(d2,col=clusterCols[1], lwd=2)
dev.off()

t.test(2^(central$Cor.RPL),2^(perpheral$Cor.RPL),alternative = "two.sided", var.equal = FALSE)



