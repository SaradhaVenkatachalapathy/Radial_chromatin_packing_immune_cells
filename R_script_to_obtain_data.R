open_and_compile <- function(dir, exp, samplen, trial,file_type){
  setwd(dir)
  file_data<-list.files()
  for (j in 1:length(file_data)){
    if (!exists("dataset")){
      if(file_type=="tsv"){
        dataset <- read.table(file_data[j], header=TRUE, sep="\t",stringsAsFactors=FALSE)
      }
      else if(file_type=="csv"){
        dataset <- read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
      }
    }
    else if (exists("dataset")){
      if(file_type=="tsv"){
        temp_dataset <- read.table(file_data[j], header=TRUE, sep="\t",stringsAsFactors=FALSE)
      }
      else if(file_type=="csv"){
        temp_dataset <- read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
      }
      if(nrow(temp_dataset)>0){
        dataset<-rbind(dataset, temp_dataset)
      }
      rm(temp_dataset)
    }
  }
  dataset$Exp<-exp
  dataset$sample<-samplen
  dataset$trial<-trial
  return(dataset)
}
open_and_compile_shells <- function(dir, exp, samplen, trial,file_type){
  setwd(dir)
  file_data<-list.files()
  for (j in 1:length(file_data)){
    if (!exists("dataset")){
      if(file_type=="tsv"){
        dataset <- read.table(file_data[j], header=TRUE, sep="\t",stringsAsFactors=FALSE)
        dataset$nucnum<-as.character(j)
      }
      else if(file_type=="csv"){
        dataset <- read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
        dataset$nucnum<-as.character(j)
      }
    }
    else if (exists("dataset")){
      if(file_type=="tsv"){
        temp_dataset <- read.table(file_data[j], header=TRUE, sep="\t",stringsAsFactors=FALSE)
        temp_dataset$nucnum<-as.character(j)
      }
      else if(file_type=="csv"){
        temp_dataset <- read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
        temp_dataset$nucnum<-as.character(j)
      }
      if(nrow(temp_dataset)>0){
        dataset<-rbind(dataset, temp_dataset)
      }
      rm(temp_dataset)
    }
  }
  dataset$Exp<-exp
  dataset$sample<-samplen
  dataset$trial<-trial
  return(dataset)
}
shell_volfraction<- function(dirg,diri, exp, sample, trial, dirres){
  
  setwd(diri)
  file_data<-list.files()
  
  geo<-open_and_compile_shells(dirg,exp, sample,trial,"tsv")
  int<-open_and_compile_shells(diri,exp, sample,trial,"tsv")
  
  intden_nom<-as.data.frame(matrix(nrow=length(file_data),ncol=23))
  colnames(intden_nom)<-c("spheroidnum","Shell0","Shell1","Shell2","Shell3","Shell4","Shell5","Shell6","Shell7","Shell8","Shell9","Shell10",
                          "Shell11","Shell12","Shell13","Shell14","Shell15","Shell16","Shell17","Shell18","Shell19","Shell20", "totint")
  intden_nom$nuc<-file_data
  vol_nom<-as.data.frame(matrix(nrow=length(file_data),ncol=24))
  colnames(vol_nom)<-c("spheroidnum","Shell0","Shell1","Shell2","Shell3","Shell4","Shell5","Shell6","Shell7","Shell8","Shell9","Shell10",
                       "Shell11","Shell12","Shell13","Shell14","Shell15","Shell16","Shell17","Shell18","Shell19","Shell20","totvol", "sphere")
  vol_nom$nuc<-file_data
  volring_nom<-as.data.frame(matrix(nrow=length(file_data),ncol=22))
  colnames(volring_nom)<-c("spheroidnum","Shell0","Shell1","Shell2","Shell3","Shell4","Shell5","Shell6","Shell7","Shell8","Shell9","Shell10",
                           "Shell11","Shell12","Shell13","Shell14","Shell15","Shell16","Shell17","Shell18","Shell19","Shell20")
  
  for( i in 1:length(file_data)){
    intsub<-subset(int,int$nucnum==i)
    intden_nom[i,23]<-intsub$IntDen[1]
    volsub<-subset(geo,geo$nucnum==i)
    vol_nom[i,23]<-volsub$Vol..unit.[1]
    vol_nom[i,24]<-volsub$Spher..unit.[1]
    
  }   
  for( i in 1:length(file_data)){
    intsub<-subset(int,int$nucnum==i)
    for(k in 1:(nrow(intsub)-1)){
      intsub$IntDen[k]<-intsub$IntDen[k]-intsub$IntDen[k+1]
    }
    intden_nom[i,1]<-i
    intden_nom[i,2:(nrow(intsub)+1)]<-intsub$IntDen
    
    volsub<-subset(geo,geo$nucnum==i)
    for(k in 1:(nrow(volsub)-1)){
      volsub$Vol..unit.[k]<-volsub$Vol..unit.[k]-volsub$Vol..unit.[k+1]
    }
    volring_nom[i,1]<-i
    volring_nom[i,2:(nrow(volsub)+1)]<-volsub$Vol..unit.
    
    volsub<-subset(geo,geo$nucnum==i)
    volsub$Vol..unit. <-volsub$Vol..unit./volsub$Vol..unit.[1]
    vol_nom[i,1]<-i
    vol_nom[i,2:(nrow(volsub)+1)]<-volsub$Vol..unit.
  }    
  
  intden_nom[,2]<-intden_nom[,2]/volring_nom[,2]
  intden_nom[,3]<-intden_nom[,3]/volring_nom[,3]
  intden_nom[,4]<-intden_nom[,4]/volring_nom[,4]
  intden_nom[,5]<-intden_nom[,5]/volring_nom[,5]
  intden_nom[,6]<-intden_nom[,6]/volring_nom[,6]
  intden_nom[,7]<-intden_nom[,7]/volring_nom[,7]
  intden_nom[,8]<-intden_nom[,8]/volring_nom[,8]
  intden_nom[,9]<-intden_nom[,9]/volring_nom[,9]
  intden_nom[,10]<-intden_nom[,10]/volring_nom[,10]
  intden_nom[,11]<-intden_nom[,11]/volring_nom[,11]
  intden_nom[,12]<-intden_nom[,12]/volring_nom[,12]
  intden_nom[,13]<-intden_nom[,13]/volring_nom[,13]
  intden_nom[,14]<-intden_nom[,14]/volring_nom[,14]
  intden_nom[,15]<-intden_nom[,15]/volring_nom[,15]
  intden_nom[,16]<-intden_nom[,16]/volring_nom[,16]
  intden_nom[,17]<-intden_nom[,17]/volring_nom[,17]
  intden_nom[,18]<-intden_nom[,18]/volring_nom[,18]
  intden_nom[,19]<-intden_nom[,19]/volring_nom[,19]
  intden_nom[,20]<-intden_nom[,20]/volring_nom[,20]
  intden_nom[,21]<-intden_nom[,21]/volring_nom[,21]
  intden_nom[,22]<-intden_nom[,22]/volring_nom[,22]
  for( i in 1:nrow(intden_nom)){
    intden_nom[i,2:22]<-intden_nom[i,2:22]/max(intden_nom[i,2:22],na.rm=T)
  }
  
  setwd(dirres)
  
  name=paste("Shell_volume_fraction",sample,trial,"lines.png",sep="_")
  png(filename=name, units="in",width=1, height=1 , pointsize=4, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
  plot(1:0,1:0,col="white",xlim=c(0,1),xlab="Volume fraction",ylab="Intensity fraction",las=1)
  cols<-rgb(runif(length(file_data)),runif(length(file_data)),runif(length(file_data))) 
  for( i in 1:length(file_data)){
    lines(vol_nom[i,2:22],intden_nom[i,2:22],col=cols[i],type="o",pch=18)
  }
  dev.off()
  
  plot(1:0,1:0,col="white",xlim=c(0,1),xlab="Volume fraction",ylab="Intensity fraction",las=1)
  
  fractions<-as.data.frame(matrix(nrow=length(file_data),ncol=11))
  fractions[,1]<-file_data
  colnames(fractions)<-c("nuc","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%")
  for( i in 1:length(file_data)){
    y<-na.omit(unlist(intden_nom[i,2:22]))
    x<-na.omit(unlist(vol_nom[i,2:22]))
    if(length(x)>=4){
      #lines(y~x, col = cols[i], pch = 19,type="o")
      points(y~x, col = cols[i], pch = 19)
      fit <- smooth.spline(x, y)
      prd<-predict(fit, x)$y
      lines(prd~x,col=cols[i])
      
      fractions[i,2]<-predict(fit, 0.1)$y
      fractions[i,3]<-predict(fit, 0.2)$y
      fractions[i,4]<-predict(fit, 0.3)$y
      fractions[i,5]<-predict(fit, 0.4)$y
      fractions[i,6]<-predict(fit, 0.5)$y
      fractions[i,7]<-predict(fit, 0.6)$y
      fractions[i,8]<-predict(fit, 0.7)$y
      fractions[i,9]<-predict(fit, 0.8)$y
      fractions[i,10]<-predict(fit, 0.9)$y
      fractions[i,11]<-predict(fit, 1)$y
    }
    
  }
  
  library(gplots)
  library(RColorBrewer)
  Colors=brewer.pal(11,"Spectral")
  
  fractions<-fractions[order(fractions[,2]),]
  row.names(fractions)<-fractions[,1]
  temp<-na.omit(as.matrix(fractions[,2:11]))
  row.names(fractions)<-fractions[,1]
  
  name=paste("Shell_volume_fraction",sample,trial,"heatmap.png",sep="_")
  png(filename=name, units="in",width=2.5, height=2 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
  data_heat<-heatmap.2(temp,scale="none",Rowv = TRUE, Colv=FALSE, dendrogram = "row", density.info = "none",trace="none",col=Colors, 
                       labRow = FALSE)
  dev.off()
  
  fractions$Exp<-exp
  fractions$sample<-sample
  fractions$trial<-trial
  return(fractions)
}

#3dcell_intensity_coro1a
{
  
  geometrical_data<-open_and_compile("E:/immune_cells/190115_11J2_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch2/", 
                                     "190115_11J2_CD4T_488Coro1A_555RPL10A", "T_cells","T1","tsv")
  temp<-open_and_compile("E:/immune_cells/190116_11J2_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch2/", 
                         "190116_11J2_CD4T_488Coro1A_555RPL10A", "T_cells","T2","tsv")
  geometrical_data<-rbind(geometrical_data,temp)
  temp<-open_and_compile("E:/immune_cells/190118_11J2_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch2/", 
                         "190118_11J2_CD4T_488Coro1A_555RPL10A", "T_cells","T3","tsv")
  geometrical_data<-rbind(geometrical_data,temp)
  temp<-open_and_compile("E:/immune_cells/190121_11J2_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch2/", 
                         "190121_11J2_CD4T_488Coro1A_555RPL10A", "T_cells","T4","tsv")
  geometrical_data<-rbind(geometrical_data,temp)
  temp<-open_and_compile("E:/immune_cells/190211_11J1_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch2/", 
                         "190211_11J1_CD4T_488Coro1A_555RPL10A", "T_cells","T5","tsv")
  geometrical_data<-rbind(geometrical_data,temp)
  temp<-open_and_compile("E:/immune_cells/190212_11J1_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch2/", 
                         "190212_11J1_CD4T_488Coro1A_555RPL10A", "T_cells","T6","tsv")
  geometrical_data<-rbind(geometrical_data,temp)
  temp<-open_and_compile("E:/immune_cells/190214_11J2_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch2/", 
                         "190214_11J2_CD4T_488Coro1A_555RPL10A", "T_cells","T7","tsv")
  geometrical_data<-rbind(geometrical_data,temp)
  temp<-open_and_compile("E:/immune_cells/190225_11J1_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch2/", 
                         "190225_11J1_CD4T_488Coro1A_555RPL10A", "T_cells","T8","tsv")
  geometrical_data<-rbind(geometrical_data,temp)
  temp<-open_and_compile("E:/immune_cells/190225_11J2_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch2/", 
                         "190225_11J2_CD4T_488Coro1A_555RPL10A", "T_cells","T9","tsv")
  geometrical_data<-rbind(geometrical_data,temp)
  
  coro1a<-geometrical_data
  rm(geometrical_data,temp)
  
}

#3dcell_intensity_rpl10
{
  
  geometrical_data<-open_and_compile("E:/immune_cells/190115_11J2_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch3/", 
                                     "190115_11J2_CD4T_488Coro1A_555RPL10A", "T_cells","T1","tsv")
  temp<-open_and_compile("E:/immune_cells/190116_11J2_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch3/", 
                         "190116_11J2_CD4T_488Coro1A_555RPL10A", "T_cells","T2","tsv")
  geometrical_data<-rbind(geometrical_data,temp)
  temp<-open_and_compile("E:/immune_cells/190118_11J2_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch3/", 
                         "190118_11J2_CD4T_488Coro1A_555RPL10A", "T_cells","T3","tsv")
  geometrical_data<-rbind(geometrical_data,temp)
  temp<-open_and_compile("E:/immune_cells/190121_11J2_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch3/", 
                         "190121_11J2_CD4T_488Coro1A_555RPL10A", "T_cells","T4","tsv")
  geometrical_data<-rbind(geometrical_data,temp)
  temp<-open_and_compile("E:/immune_cells/190211_11J1_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch3/", 
                         "190211_11J1_CD4T_488Coro1A_555RPL10A", "T_cells","T5","tsv")
  geometrical_data<-rbind(geometrical_data,temp)
  temp<-open_and_compile("E:/immune_cells/190212_11J1_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch3/", 
                         "190212_11J1_CD4T_488Coro1A_555RPL10A", "T_cells","T6","tsv")
  geometrical_data<-rbind(geometrical_data,temp)
  temp<-open_and_compile("E:/immune_cells/190214_11J2_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch3/", 
                         "190214_11J2_CD4T_488Coro1A_555RPL10A", "T_cells","T7","tsv")
  geometrical_data<-rbind(geometrical_data,temp)
  temp<-open_and_compile("E:/immune_cells/190225_11J1_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch3/", 
                         "190225_11J1_CD4T_488Coro1A_555RPL10A", "T_cells","T8","tsv")
  geometrical_data<-rbind(geometrical_data,temp)
  temp<-open_and_compile("E:/immune_cells/190225_11J2_CD4T_488Coro1A_555RPL10A/cell_2microns_measure/ch3/", 
                         "190225_11J2_CD4T_488Coro1A_555RPL10A", "T_cells","T9","tsv")
  geometrical_data<-rbind(geometrical_data,temp)
  
  rp10a<-geometrical_data
  rm(geometrical_data,temp)
  
}

#shellfracs
{
  dird<-"E:/immune_cells/R_analysis/shell_analysis"
  dir.create(dird)
  setwd(dird)
  
  shell_volfraction_D6_t1<-shell_volfraction("E:/immune_cells/190115_11J2_CD4T_488Coro1A_555RPL10A/Shells/geometric_measures/", 
                                             "E:/immune_cells/190115_11J2_CD4T_488Coro1A_555RPL10A/Shells/intensity_shells/", 
                                             "190115_11J2_CD4T_488Coro1A_555RPL10A", "T_Cells","T1",dird)
  frac<-shell_volfraction_D6_t1
  shell_volfraction_D6_t1<-shell_volfraction("E:/immune_cells/190116_11J2_CD4T_488Coro1A_555RPL10A/Shells/geometric_measures/", 
                                             "E:/immune_cells/190116_11J2_CD4T_488Coro1A_555RPL10A/Shells/intensity_shells/", 
                                             "190116_11J2_CD4T_488Coro1A_555RPL10A", "T_Cells","T2",dird)
  frac<-rbind(frac,shell_volfraction_D6_t1)
  shell_volfraction_D6_t1<-shell_volfraction("E:/immune_cells/190118_11J2_CD4T_488Coro1A_555RPL10A/Shells/geometric_measures/", 
                                             "E:/immune_cells/190118_11J2_CD4T_488Coro1A_555RPL10A/Shells/intensity_shells/", 
                                             "190118_11J2_CD4T_488Coro1A_555RPL10A", "T_Cells","T3",dird)
  frac<-rbind(frac,shell_volfraction_D6_t1)
  shell_volfraction_D6_t1<-shell_volfraction("E:/immune_cells/190121_11J2_CD4T_488Coro1A_555RPL10A/Shells/geometric_measures/", 
                                             "E:/immune_cells/190121_11J2_CD4T_488Coro1A_555RPL10A/Shells/intensity_shells/", 
                                             "190121_11J2_CD4T_488Coro1A_555RPL10A", "T_Cells","T4",dird)
  frac<-rbind(frac,shell_volfraction_D6_t1)
  shell_volfraction_D6_t1<-shell_volfraction("E:/immune_cells/190211_11J1_CD4T_488Coro1A_555RPL10A/Shells/geometric_measures/", 
                                             "E:/immune_cells/190211_11J1_CD4T_488Coro1A_555RPL10A/Shells/intensity_shells/", 
                                             "190211_11J1_CD4T_488Coro1A_555RPL10A", "T_Cells","T5",dird)
  frac<-rbind(frac,shell_volfraction_D6_t1)
  
  shell_volfraction_D6_t1<-shell_volfraction("E:/immune_cells/190212_11J1_CD4T_488Coro1A_555RPL10A/Shells/geometric_measures/", 
                                             "E:/immune_cells/190212_11J1_CD4T_488Coro1A_555RPL10A/Shells/intensity_shells/", 
                                             "190212_11J1_CD4T_488Coro1A_555RPL10A", "T_Cells","T6",dird)
  frac<-rbind(frac,shell_volfraction_D6_t1)
  
  shell_volfraction_D6_t1<-shell_volfraction("E:/immune_cells/190214_11J2_CD4T_488Coro1A_555RPL10A/Shells/geometric_measures/", 
                                             "E:/immune_cells/190214_11J2_CD4T_488Coro1A_555RPL10A/Shells/intensity_shells/", 
                                             "190214_11J2_CD4T_488Coro1A_555RPL10A", "T_Cells","T7",dird)
  frac<-rbind(frac,shell_volfraction_D6_t1)
  
  shell_volfraction_D6_t1<-shell_volfraction("E:/immune_cells/190225_11J1_CD4T_488Coro1A_555RPL10A/Shells/geometric_measures/", 
                                             "E:/immune_cells/190225_11J1_CD4T_488Coro1A_555RPL10A/Shells/intensity_shells/", 
                                             "190225_11J1_CD4T_488Coro1A_555RPL10A", "T_Cells","T8",dird)
  frac<-rbind(frac,shell_volfraction_D6_t1)
  shell_volfraction_D6_t1<-shell_volfraction("E:/immune_cells/190225_11J2_CD4T_488Coro1A_555RPL10A/Shells/geometric_measures/", 
                                             "E:/immune_cells/190225_11J2_CD4T_488Coro1A_555RPL10A/Shells/intensity_shells/", 
                                             "190225_11J2_CD4T_488Coro1A_555RPL10A", "T_Cells","T9",dird)
  frac<-rbind(frac,shell_volfraction_D6_t1)
}
rm(shell_volfraction_D6_t1)
colnames(frac)[1]<-"Label"
frac$Label<-substring(frac$Label,3, nchar(frac$Label)-8)

colnames(coro1a)[14]<-"coro1a_mean"
colnames(rp10a)[14]<-"rpl10a_mean"

fractions<-merge(coro1a[,c(3,14,19)],rp10a[,c(3,14)], by="Label")
fractions$logCor_rpl<-log(fractions$coro1a_mean/fractions$rpl10a_mean, base=2)
fractions<-merge(fractions,frac, by="Label")
rownames(fractions)<-fractions$Label

#setwd("E:/immune_cells/R_analysis/")
write.csv(fractions,file="protein_ratios4.csv")
