library("ShortRead")
library("rtracklayer")
library("Rsamtools")
library("IRanges")
library("parallel")

f.end.profiler <- function(f.filename,f.filename.rev="",f.savename,f.annotation,f.sense=TRUE,f.end="5end", f.upstream= 1000, f.downstream=1000, f.savediskspace=FALSE) {
  
  # f.filename - BigWig file 
  # f.filename.rev - if empty, f.filename will be used on both strand, if BigWig file, f.filename will be used on fw starnd, f.filename.rev will be used on the reverse strand
  # f.savename - wil be added .5profile.Rdat OR .3profile.Rdat OR .5/3profileS.Rdat or 5/3profileAS.Rdat
  # f.annotation - Rdat file from "~/Desktop/R/annotation" folder
  # f.sense - if TRUE, sense data will be generated, FALSE -AS data will be genetarted - it is used only if f.filename.rev is not empty  
  # f.end=    "5end" - fragments will be collected upstream and downstream form 5' end of the features
  #           "3end" - fragments will be collected upstream and downstream form 3' end of the features
  #           "middle" - fragments will be collected upstream and downstream from middle of the features
  # f.upstream - fragment length collected upstream in bp
  # f.downstream - fragment length collected downstream in bp
  # f.savediskspace=TRUE - will use Rle format in the result file, plotter needs more time to process it, but smaller file size
  
  
  
  # f.annotation: "~/Desktop/R/annotation/S.pombe.EF2/EF2.mRNA.Rdat"
  #               "~/Desktop/R/annotation/S.pombe.2007/sp07.mRNA.Rdat"
  #               "~/Desktop/R/annotation/S.pombe.EF2/EF2.mRNA.nucl.WT040213.Rdat"  OR 260613.Rdat
  #               "~/Desktop/R/annotation/S.pombe.2007/sp07.mRNA.nucl.WT040213.Rdat"  OR 2606613.Rdat
  #               "~/Desktop/R/annotation/N.crassa_or74a_10/Nc10.genes.Rdat"
  #               "~/Desktop/R/annotation/S.pombe.GDB090511/GDB.intron.Rdat"    this is equivalent with the EF2 annotation
  
  ###  
  #setwd("~/Documents/seq_data/seq_0008_260613_mononucleosome_Sp_WT_hrp1d_hrp3d_set2d_set1d_set9d_dCD/BigWig_files")
  #f.filename<-"N.sp_st3230_RED1RE_Rep1_polyA_Illumina2_ERCC_RNAseq.f.bw"
  #f.filename.rev<-"N.sp_st3230_RED1RE_Rep1_polyA_Illumina2_ERCC_RNAseq.r.bw"
  #f.savename<- "N.sp_st3230_RED1RE_Rep1_polyA_Illumina2_ERCC_RNAseq"
  
  #f.annotation<-"~/Desktop/R/annotation/S.pombe.EF2/EF2.mRNA.Rdat"
  #f.upstream= 500
  #f.downstream=500
  #f.savediskspace=FALSE
  #f.end="5end"
  #f.sense=TRUE
  ###  
  
  f.binnumber=0
  f.binlength=0
  f.excluded=integer(0)
  
  Rle.data <- import(f.filename,as="RleList")               # load data file
  if(f.filename.rev!=""){
    Rle.data.rev<-import(f.filename.rev,as="RleList")       # load reverse data file
    if(f.sense==TRUE){                                # Sense 
      if(f.end=="5end"){
        f.save.end<-".5profileS.Rdat"
      }
      if(f.end=="3end"){
        f.save.end<-".3profileS.Rdat"
      }
      if(f.end=="middle"){
        f.save.end<-".MprofileS.Rdat"
      }
    }else{                                            # Antisense
      if(f.end=="5end"){
        f.save.end<-".5profileAS.Rdat"
      }
      if(f.end=="3end"){
        f.save.end<-".3profileAS.Rdat"
      }
      if(f.end=="middle"){
        f.save.end<-".MprofileAS.Rdat"
      }
      Rle.data.temp<-Rle.data                         # swapping fw and rev strand information
      Rle.data<-Rle.data.rev
      Rle.data.rev<-Rle.data.temp
      rm(Rle.data.temp)
    }
  }else{
    if(f.end=="5end"){
      f.save.end<-".5profile.Rdat"
    }
    if(f.end=="3end"){
      f.save.end<-".3profile.Rdat"
    }
    if(f.end=="middle"){
      f.save.end<-".Mprofile.Rdat"
    }
  }
  
  load(f.annotation)                                      # load annotation GRanges object  ann.gr
  if(f.end=="5end"){
    suppressWarnings(temp.gr<-promoters(ann.gr,upstream=f.upstream, downstream=f.downstream))    # modified ranges, 5 end, suppress warning about trim
  }
  if(f.end=="3end"){
    suppressWarnings(tempup.gr<-flank(ann.gr,width=-f.upstream, start=FALSE, both=FALSE))    # modified ranges, 3 end suppress warning about trim
    suppressWarnings(temp.gr<-resize(tempup.gr,width=f.upstream+f.downstream, fix="start"))
    rm(tempup.gr)
  }
  if(f.end=="middle"){
    tempmid.gr<-ann.gr
    start(tempmid.gr)<-(start(ann.gr)+end(ann.gr))%/%2
    end(tempmid.gr)<-(start(ann.gr)+end(ann.gr))%/%2
    suppressWarnings(temp.gr<-promoters(tempmid.gr,upstream=f.upstream, downstream=f.downstream))    # modified ranges, middle, suppress warning about trim
    rm(tempmid.gr)
  }
  
  temp.gr<-trim(temp.gr)                                                            # trim out-of-bound ranges
  
  temp.excluded<-(1:length(temp.gr))[width(temp.gr)==0]                            # if width is 0 (start is 1 larger than end), set start to 1 and width to 1
  start(temp.gr[temp.excluded])<-1
  width(temp.gr[temp.excluded])<-1
  f.excluded<-c(f.excluded,temp.excluded)                                           # add these records to the excluded list
  
  profile.rle<-RleList(mclapply(1:length(temp.gr), function(x){
    if(as.character(strand(temp.gr[x])=="+")){                                                #fw starnd faetures from fw starnd file
      Rle.data[[as.character(seqnames(temp.gr[x]))]][start(temp.gr[x]):end(temp.gr[x])]
    }else{
      if(f.filename.rev==""){                                                                 #reverse starnd faetures also from the same file if no rev file exist
        Rle.data[[as.character(seqnames(temp.gr[x]))]][end(temp.gr[x]):start(temp.gr[x])]
      }else{
        Rle.data.rev[[as.character(seqnames(temp.gr[x]))]][end(temp.gr[x]):start(temp.gr[x])]     #if reverse file exist, reverse starnd faetures from the rev file 
      }
    }    
  }))
  
  rle.length<-sum(runLength(profile.rle))             # checking if all features are full length, non full length are excluded and exluded feature numbers are saved in f.excluded
  if(sum(rle.length!=(f.upstream+f.downstream))>0){
    temp.excluded<-(1:length(ann.gr))[rle.length!=(f.upstream+f.downstream)]
    f.excluded<-c(f.excluded,temp.excluded)
    f.excluded<-unique(f.excluded)
    print(paste("Lengths anomalie! - Excluded features:", length(f.excluded)))
    print(as.character(mcols(ann.gr[f.excluded])[["Name"]]))
    profile.rle<-profile.rle[-f.excluded]
  }
  
  if(f.savediskspace==FALSE){
    profile.rle<-matrix(unlist(mclapply(1:length(profile.rle),function(x) as.numeric(profile.rle[[x]])))
                        ,nrow=length(profile.rle), ncol=f.upstream+f.downstream, byrow=TRUE)  #converting Rle list to matrix of numerical values
  }
  
  save(profile.rle,f.annotation,f.upstream,f.binnumber,f.binlength,f.downstream,f.savediskspace,f.excluded, file=paste(f.savename, f.save.end, sep=""))
  
}


library("ShortRead")
library("rtracklayer")
library("Rsamtools")
library("IRanges")
library("parallel")

f.gene.profiler <- function(f.filename,f.filename.rev="",f.savename,f.annotation,f.sense=TRUE,f.binnumber=30, f.binlength=40, f.upstream= 1000, f.downstream=1000, f.savediskspace=FALSE, f.cores=4) {
  
  # f.filename - BigWig file 
  # f.filename.rev - if empty, f.filename will be used on both strand, if BigWig file, f.filename will be used on fw starnd, f.filename.rev will be used on the reverse strand
  # f.savename - wil be added .profile.Rdat OR .profileS.Rdat OR .profileAS.Rdat
  # f.annotation - Rdat file from "~/Desktop/R/annotation" folder
  # f.sense - if TRUE, sense data will be generated, FALSE -AS data will be genetarted - it is used only if f.filename.rev is not empty
  # f.binnumber - number of bins the feature will be divided
  # f.binlength - length of bins in bp - f.binnumber x f.binnlength = average feature length in bp
  # f.upstream - fragment length collected upstream in bp
  # f.downstream - fragment length collected downstream in bp
  # f.savediskspace=TRUE - will use Rle format in the result file, plotter needs more time to process it, but smaller file size
  # f.cores - no of processor cores
  
  
  # f.annotation: "~/Desktop/R/annotation/S.pombe.EF2/EF2.mRNA.Rdat"
  #               "~/Desktop/R/annotation/S.pombe.2007/sp07.mRNA.Rdat"
  #               "~/Desktop/R/annotation/S.pombe.EF2/EF2.mRNA.nucl.WT040213.Rdat"  OR 260613.Rdat
  #               "~/Desktop/R/annotation/S.pombe.2007/sp07.mRNA.nucl.WT040213.Rdat"  OR 2606613.Rdat
  #               "~/Desktop/R/annotation/N.crassa_or74a_10/Nc10.genes.Rdat"
  #               "~/Desktop/R/annotation/S.pombe.EF2/EF2.introns20to400.Rdat"    
  
  ###  
  #setwd("~/Documents/seq_data/seq_0008_260613_mononucleosome_Sp_WT_hrp1d_hrp3d_set2d_set1d_set9d_dCD/BigWig_files")
  #f.filename<- "N.sp_st344_WT_polyA_RNAseq_pos.bw"
  #f.filename.rev="N.sp_st344_WT_polyA_RNAseq_neg.bw"
  #f.savename<- "N.sp_st344_WT"
  #f.annotation<-"~/Desktop/R/annotation/S.pombe.EF2/EF2.mRNA.Rdat" 
  #f.sense=TRUE
  #f.binnumber=30
  #f.binlength=40
  #f.upstream= 1000
  #f.downstream=1000
  #f.savediskspace=FALSE
  #f.cores <- 4
  ###  
  
  f.excluded=integer(0)
  
  Rle.data <- import(f.filename,as="RleList")               # load data file
  if(f.filename.rev!=""){
    Rle.data.rev<-import(f.filename.rev,as="RleList")       # load reverse data file
    if(f.sense==TRUE){                                # Sense 
      f.save.end<-".profileS.Rdat"
    }else{                                            # Antisense
      f.save.end<-".profileAS.Rdat"
      Rle.data.temp<-Rle.data                         # swapping fw and rev strand information
      Rle.data<-Rle.data.rev
      Rle.data.rev<-Rle.data.temp
      rm(Rle.data.temp)
    }
  }else{
    f.save.end<-".profile.Rdat"
  }
  
  load(f.annotation)                                      # load annotation GRanges object  ann.gr
  
  ######################################## 5' probes
  
  if(f.upstream>0){
    suppressWarnings(temp.gr<-promoters(ann.gr,upstream=f.upstream, downstream=0))    # modified ranges, 5 end, suppress warning about trim
    temp.gr<-trim(temp.gr)                                                            # trim out-of-bound ranges
    
    temp.excluded<-(1:length(temp.gr))[width(temp.gr)==0]                            # if width is 0 (start is 1 larger than end), set start to 1 and width to 1
    start(temp.gr[temp.excluded])<-1
    width(temp.gr[temp.excluded])<-1
    f.excluded<-c(f.excluded,temp.excluded)                                           # add these records to the excluded list
    
    profile5.rle<-RleList(mclapply(1:length(temp.gr), function(x){
      if(as.character(strand(temp.gr[x])=="+")){                                                #fw starnd faetures from fw starnd file
        Rle.data[[as.character(seqnames(temp.gr[x]))]][start(temp.gr[x]):end(temp.gr[x])]
      }else{
        if(f.filename.rev==""){                                                                 #reverse starnd faetures also from the same file if no rev file exist
          Rle.data[[as.character(seqnames(temp.gr[x]))]][end(temp.gr[x]):start(temp.gr[x])]
        }else{
          Rle.data.rev[[as.character(seqnames(temp.gr[x]))]][end(temp.gr[x]):start(temp.gr[x])]     #if reverse file exist, reverse starnd faetures from the rev file 
        }
      }    
    }, mc.cores=f.cores))  
  }else{                                                # if f.upstream=0, we have to define profile.rle
    profile5.rle=RleList()
  }
  
  
  ######################################## feature into f.binnumber bins
  
  profilemid.rle<-RleList(mclapply(1:length(ann.gr), function(x){
    if(as.character(strand(ann.gr[x])=="+")){                                                #fw strand faetures from fw strand file
      f.genestrech(as.numeric(Rle.data[[as.character(seqnames(ann.gr[x]))]][start(ann.gr[x]):end(ann.gr[x])]),f.binnumber) 
    }else{
      if(f.filename.rev==""){                                                                 #reverse strand faetures also from the same file if no rev file exist
        f.genestrech(as.numeric(Rle.data[[as.character(seqnames(ann.gr[x]))]][end(ann.gr[x]):start(ann.gr[x])]),f.binnumber)
      }else{
        f.genestrech(as.numeric(Rle.data.rev[[as.character(seqnames(ann.gr[x]))]][end(ann.gr[x]):start(ann.gr[x])]),f.binnumber)     #if reverse file exist, reverse strand faetures from the rev file 
      }
    }    
  }, mc.cores=f.cores))  
  
  
  ######################################## 3' probes
  
  if(f.downstream>0){
    suppressWarnings(temp.gr<-flank(ann.gr,width=f.downstream, start=FALSE, both=FALSE))    # modified ranges, 3 end suppress warning about trim
    temp.gr<-trim(temp.gr)                                                            # trim out-of-bound ranges
    
    temp.excluded<-(1:length(temp.gr))[width(temp.gr)==0]                            # if width is 0 (start is 1 larger than end), set start to 1 and width to 1
    start(temp.gr[temp.excluded])<-1
    width(temp.gr[temp.excluded])<-1
    f.excluded<-c(f.excluded,temp.excluded)                                           # add these records to the excluded list
    
    profile3.rle<-RleList(mclapply(1:length(temp.gr), function(x){
      if(as.character(strand(temp.gr[x])=="+")){                                                #fw starnd faetures from fw starnd file
        Rle.data[[as.character(seqnames(temp.gr[x]))]][start(temp.gr[x]):end(temp.gr[x])]
      }else{
        if(f.filename.rev==""){                                                                 #reverse starnd faetures also from the same file if no rev file exist
          Rle.data[[as.character(seqnames(temp.gr[x]))]][end(temp.gr[x]):start(temp.gr[x])]
        }else{
          Rle.data.rev[[as.character(seqnames(temp.gr[x]))]][end(temp.gr[x]):start(temp.gr[x])]     #if reverse file exist, reverse starnd faetures from the rev file 
        }
      }    
    }, mc.cores=f.cores))  
  }else{                                                # if f.upstream=0, we have to define profile.rle
    profile3.rle=RleList()
  }
  
  
  ############## combining 5' , feature and 3' rle
  if(f.upstream>0 & f.downstream>0){
    profile.rle<-RleList(mclapply(1:length(profilemid.rle), function(x) c(profile5.rle[[x]],profilemid.rle[[x]],profile3.rle[[x]]),mc.cores=f.cores))
  }else{
    if(f.upstream>0){
      profile.rle<-RleList(mclapply(1:length(profilemid.rle), function(x) c(profile5.rle[[x]],profilemid.rle[[x]]),mc.cores=f.cores))
    }else{
      if(f.downstream>0){
        profile.rle<-RleList(mclapply(1:length(profilemid.rle), function(x) c(profilemid.rle[[x]],profile3.rle[[x]]),mc.cores=f.cores))
      }else{
        profile.rle<-profilemid.rle
      }
    }
  }
  
  
  rle.length<-sum(runLength(profile.rle))             # checking if all features are full length, non full length are excluded and exluded feature numbers are saved in f.excluded
  if(sum(rle.length!=(f.upstream+f.binnumber+f.downstream))>0){
    temp.excluded<-(1:length(ann.gr))[rle.length!=(f.upstream+f.binnumber+f.downstream)]
    f.excluded<-c(f.excluded,temp.excluded)
    f.excluded<-unique(f.excluded)
    print(paste("Lengths anomalie! - Excluded features:", length(f.excluded)))
    print(as.character(mcols(ann.gr[f.excluded])[["Name"]]))
    profile.rle<-profile.rle[-f.excluded]
  }
  
  
  if(f.savediskspace==FALSE){
    profile.rle<-matrix(unlist(mclapply(1:length(profile.rle),function(x) as.numeric(profile.rle[[x]]), mc.cores=2*f.cores))
                        ,nrow=length(profile.rle), ncol=f.upstream+f.binnumber+f.downstream, byrow=TRUE)  #converting Rle list to matrix of numerical values
  }
  
  save(profile.rle,f.annotation,f.upstream,f.binnumber,f.binlength,f.downstream,f.savediskspace,f.excluded, file=paste(f.savename, f.save.end, sep=""))
  
}


f.genestrech <- function(f.feature,f.binnumber){        # devide an f.feature to f.binnumber bins
  
  #f.binnumber<-as.numeric(f.binnumber)
  f.FElength<- length(f.feature)                        # lenght of the feature
  s.ratio<- f.binnumber/f.FElength                      # theoretical ratio of how many bins will be occupied by one bp of the feature
  a.feature <- rep(s.ratio,f.FElength)                  # s.ratio to every bp of the feature
  a.binnumber <- rep(0.000,f.binnumber)                 # how many bins did we fill up
  FEbins <- rep(0.00, f.binnumber)                      # FEbins will contain the final bins of the feature
  
  
  i <- 1										# ith bp in the feature
  j <- 1										# jth element in the bins
  
  while(a.feature[f.FElength]>0){				# repeat until we didn't divide the feature completely
    if(a.feature[i]>(1-a.binnumber[j])){ 			# if we have to divide the ith bp in the feature 
      FEbins[j]<-FEbins[j]+(1-a.binnumber[j])*f.feature[i]    #FEbins j. elembe az ORF i. elemenek egy resze
      a.feature[i] <-a.feature[i]-(1-a.binnumber[j])              #faeture array i. elemet csokkenteni a resszel amit FEbins-be beletettunk
      if(a.feature[i]<0.0001) a.feature[i]<-0		             #kerekitesi hiba miatt 0 nem mindig 0
      a.binnumber[j]<- 1						             #FEbins j. elem tele van
      j<-j+1								             #ugras az FEbins j+1. elemre - i marad mert a.feature[i] meg >0
    }else{ 									# ha az i. elem osztatlanul megy a j. elembe
      FEbins[j]<-FEbins[j]+(f.feature[i]*a.feature[i])         #FEbins j. elembe az ORF i. elemenek a maradeka
      a.binnumber[j]<- a.binnumber[j]+a.feature[i]               #FEbins array j. elemet novelni a beletett resszel (maximum 1 az eredmeny)
      a.feature[i]<-0							             #ORF i. eleme ures
      i<-i+1								             #ugras az ORF i+1. elemere
      if(a.binnumber[j]==1) j<-j+1			             #ha a FEbins j. eleme tele van, ugras a j+1. elemre
    }
  }
  
  return(FEbins)
}



library("ShortRead")
library("rtracklayer")
library(stats)

f.BW.norm <- function(f.filename1, f.filename2="",f.mode="quantile",f.normvalue=100,f.chr="",f.normfile="", f.normfilter="", 
                      f.quantile=0.5, f.qfilter=c(0.9,10),f.replace0to=0) {
  
  
  # f.filename1 - BigWig file to normalize
  # f.filename2 - if it is not empty, will be added to file1 and used these values for normalization - e.g. RNA data Fw and Rev data files 
  # f.mode -  "coverage" - calculates the total genome coverage (or estimates from a single chromosome) and normalizes to the f.norm.value
  #                         ( take care e.g. the ribosomal DNA! -do you want it in or not!)
  #           "quantile" - calculates the quantile coverage of a given chromosome and normlizes to the f.norm.value
  #                        if .Rdat file is present, the median of each row will be calculated (filtered by the f.qfilter if given) and the quantile of the rows will be the normalization base value
  #           "mean"     - calculates the mean coverage of the whole genome(f.chr="") or a given chromosome, and normlizes to the f.norm.value
  #                        similar to the coverage mode - the numbers are different but the same principle 
  #           "GmeanGmean" - only if f.normfile is provided (usually Rdat file) - calculates the geometric mean of each feature and then the geom mean of the dataset - at the moment only possible with the f.normfile .Rdat normalization
  #                      - use only if the data is NOT in log2 scale! if it is already log2 scale, normal mean will give geometric mean
  #                       - values will be given back in decimal scale!
  #           "MedGmean" - only if f.normfile is provided (usually Rdat file)- calculates the median of each feature and then the geom mean of the dataset  - at the moment only possible with the f.normfile .Rdat normalization
  #                      - use only if the data is NOT in log2 scale! if it is already log2 scale, normal mean will give geometric mean
  #                       - values will be given back in decimal scale!
  
  #                     
  # f.normvalue - this will be the normalized value (coverage/median coverage, etc.)
  # f.chr - if="", whole genome will be used as reference, otherwise program will try to use as a chromosome identification (ERROR if don't match) and use the given chr as refernce
  #               if f.mode is "quantile", use!!!!! f.chr, otherwise nonsense!!!
  # f.normfile - if it is not empty, this file will be used to normalize the dataset
  #               allowed file.extensions:
  #                   - .peak.bw OR peak.BigWig - can be used median, mean, or any quantile of the peaks as reference
  #                   - .Rdat files, generated with geneprofiler - will use every feature mean/median as reference
  #                   -             if f.normfilter is not empty, will use the filter file for .Rdat files - must use the same annotation!
  # f.normfilter - only significant if f.normfile is an .Rdat file, will use this filter for calculating mean/median values
  # f.qfilter - c(a,b) - only for .Rdat quantile normalization, it will filter every value from a row which is below the (a percentile)/b of (only the row's percentile, not the total dataset percentile)
  #             best is to use c(0.9,10)-filter out everything that is smaller than 10% of the 90 percentile - e.g filter out introns or annotation that are too long
  # f.quantile - if quantile is used, can be changed to different quantile (0.5 is median, 0.95 is 95 quantile)
  # f.replace0to - only if f.mode = "GmeanGmean" or "MedGmean" (geometric mean) - 0 will be replaced with this value and this value will be normalized with the norm factor and these values will be in the normalized file (instead of 0)
  #                                                                             - this is a dangerous part of normalization, if there are amny 0s in a data, and the norm facor is high, you can really push up the avareg of the datas.
  #                                                                             - might be better to keep this 0, in this case 0 stay 0 in the normalized data, and these 0s will be not counted in the avarege 
  
  ###  
  # f.filename1="RNAseq.Sp.wt.210214.fw.BigWig"
  # f.filename2="RNAseq.Sp.wt.210214.rev.BigWig"
  # f.mode="quantile"
  # f.qfilter=c(0.9,10)
  # f.normvalue=100
  # f.normfile="RNAseq.Sp.wt.210214.profileS.Rdat"
  # f.normfilter="~/Documents/annotation/S.pombe.EF2/filters/EF2.mRNA.NoIntrons.NoHetChr.filter"
  # f.replace0to=0.5
  # f.chr=""
  # f.quantile=0.5
  ###
  
  norm.text<-list()
  Rle.data1 <- import(f.filename1,as="RleList")
  Rle.data<-Rle.data1
  if(f.filename2!=""){
    Rle.data2 <- import(f.filename2,as="RleList")
    Rle.data<-Rle.data1+Rle.data2
    norm.text<-c(norm.text,paste("Sum of two files was used to the normalization"))
    norm.text<-c(norm.text,paste("File1:",f.filename1))
    norm.text<-c(norm.text,paste("File2:",f.filename2))
  }
  
  
  ############################################################             No Normalization file            ############################################################
  ############################################################             No Normalization file            ############################################################
  ############################################################             No Normalization file            ############################################################
  if(f.normfile==""){
    
    ############################################################             mode = "coverage"            ############################################################
    if(f.mode=="coverage"){                                     
      sample.coverage<-sum(Rle.data)
      if(f.chr==""){                                            # if all genome, sum up all chromosomes!
        sample.coverage<-sum(sample.coverage)
        print(paste("Total genome coverage:",sample.coverage))
        norm.text<-c(norm.text,paste("Original total genome coverage:",sample.coverage))
        norm.text<-c(norm.text,paste("Normalized to",f.normvalue))
      }else{                                                    # otherwise use the f.chr chromosome and estimate (extrapolte to length) the coverage of the total genome
        genome.length<-sum(sum(runLength(Rle.data)))
        chr.length<-sum(runLength(Rle.data[f.chr]))
        sample.coverage<-sample.coverage[f.chr]*(genome.length/chr.length)
        print(paste("Estimated total genome coverage:",sample.coverage,"from chromosome", f.chr))
        norm.text<-c(norm.text,paste("Estimated total genome coverage:",sample.coverage,"from chromosome", f.chr))
        norm.text<-c(norm.text,paste("Normalized to",f.normvalue))
      }
      
      norm.factor<-f.normvalue/sample.coverage
    }
    
    ############################################################             mode = "quantile"            ############################################################
    if(f.mode=="quantile"){
      if(f.chr==""){
        print("f.mode='quantile' is not compatible with f.chr='' !!!")
        return()
      }
      sample.coverage<-quantile(Rle.data[f.chr],probs=f.quantile,names=FALSE)[1,1]
      print(paste("Estimated", f.quantile, "quantile genome coverage:",sample.coverage,"from chromosome", f.chr))
      norm.factor<-f.normvalue/sample.coverage
      
      norm.text<-c(norm.text,paste("Estimated", f.quantile, "quantile genome coverage:",sample.coverage,"from chromosome", f.chr))
      norm.text<-c(norm.text,paste("Normalized to",f.normvalue))
    }
    
    ############################################################             mode = "mean"            ############################################################
    
    if(f.mode=="mean"){                                     
      sample.coverage<-mean(Rle.data)
      if(f.chr==""){                                            # if all genome, weighted average (according to chr length) of  all chromosomes
        genome.length<-sum(sum(runLength(Rle.data)))
        chr.length<-sum(runLength(Rle.data[]))
        sample.coverage<-weighted.mean(sample.coverage,chr.length/genome.length)
        
        print(paste("Mean genome coverage:",sample.coverage))
        norm.text<-c(norm.text,paste("Mean genome coverage:",sample.coverage))
        norm.text<-c(norm.text,paste("Normalized to",f.normvalue))
      }else{                                                    # otherwise use the f.chr chromosome to estimate the average coverage of the total genome
        sample.coverage<-sample.coverage[f.chr]
        print(paste("Mean genome coverage:",sample.coverage,"of chromosome", f.chr))
        norm.text<-c(norm.text,paste("Mean genome coverage:",sample.coverage,"of chromosome", f.chr))
        norm.text<-c(norm.text,paste("Normalized to",f.normvalue))
      }
      norm.factor<-f.normvalue/sample.coverage
    }
    
    
  }else{
    ############################################################             Normalization file            ############################################################
    ############################################################             Normalization file            ############################################################
    ############################################################             Normalization file            ############################################################    
    
    
    ############################################################           .peak.  normalization           ############################################################
    ############################################################           .peak.  normalization           ############################################################
    ############################################################           .peak.  normalization           ############################################################
    
    if(grepl(".peak.",f.normfile)){
      
      norm.data <- import(f.normfile,as="RleList")
      norm.data <- norm.data[norm.data>0]
      
      norm.text<-c(norm.text,paste(".peak.BigWig file normalization, using ",f.normfile))
      
      ############################################################     .peak.  quantile  normalization  ############################################################
      
      if(f.mode=="quantile"){
        if(f.chr==""){
          print("f.mode='quantile' is not compatible with f.chr='' !!!")
          return()
        }
        norm.data<-norm.data[f.chr]                                         # choose one chromosome
        peak.quantile<- quantile(norm.data, probs=f.quantile,names=FALSE)   # calculating the f.quantile quantile of peak-values on chromosome f.chr
        norm.factor<- f.normvalue/peak.quantile                             # normalization factor based on peak values
        
        print(paste(f.quantile, "quantile peak values:",peak.quantile,"from chromosome", f.chr))
        norm.text<-c(norm.text,paste(f.quantile, "quantile peak values:",peak.quantile,"from chromosome", f.chr))
        norm.text<-c(norm.text,paste("Normalized to",f.normvalue))
        
        rm(norm.data)
      }
      ############################################################     .peak.  mean  normalization  ############################################################
      
      if(f.mode=="mean"){
        mean.peak<-mean(norm.data)
        if(f.chr==""){                                            # if all genome, weighted average (according to chr length) of peak values of all chromosomes
          genome.length<-sum(sum(runLength(Rle.data)))
          chr.length<-sum(runLength(Rle.data[]))
          mean.peak<-weighted.mean(mean.peak,chr.length/genome.length)
          
          print(paste("Mean peak value:",mean.peak))
          norm.text<-c(norm.text,paste("Mean peak value:",sample.coverage))
          norm.text<-c(norm.text,paste("Normalized to",f.normvalue))
        }else{                                                    # otherwise use the f.chr chromosome to estimate the average coverage of the total genome
          mean.peak<-mean.peak[f.chr]
          print(paste("Mean peak value:",mean.peak,"of chromosome", f.chr))
          norm.text<-c(norm.text,paste("Mean peak value:",mean.peak,"of chromosome", f.chr))
          norm.text<-c(norm.text,paste("Normalized to",f.normvalue))
        }
        norm.factor<-f.normvalue/mean.peak
        rm(norm.data)
      }
      
    }    
    
    ############################################################          Rdat.  normalization           ############################################################
    ############################################################          Rdat.  normalization           ############################################################
    ############################################################          Rdat.  normalization           ############################################################
    
    if(grepl(".Rdat$",f.normfile)){
      
      load(f.normfile)
      if(f.normfilter==""){       # if no filter, make a filter with all TRUE
        gr.filter<-rep(TRUE,nrow(profile.rle))
      }else{
        load(f.normfilter)        # filter is loaded into gr.filter
        if(length(f.excluded)>0){            #f.excluded is loaded with f.normfile, excluded annotation features must be excluded also from the filter 
          gr.filter<-gr.filter[-f.excluded]   # excluded values by the profiler should be excluded from the filter 
        }
        
        if(nrow(profile.rle)!=length(gr.filter)){
          print("Annotation is not compatible with the filterfile! Annotation used in Rdat file:")
          print(f.annotation)
          return()
        }
      }
      norm.data<-profile.rle[gr.filter,]
      rm(profile.rle)
      # save(profile.rle,f.annotation,f.upstream,f.binnumber,f.binlength,f.downstream,f.savediskspace,f.excluded, file=paste(f.savename, f.save.end, sep="")
      
      norm.text<-c(norm.text,paste(".Rdat file normalization, using ",f.normfile))
      
      ############################################################     Rdat  quantile  normalization  ############################################################
      
      if(f.mode=="quantile"){
        if(f.chr!=""){
          print("f.mode='quantile' with Rdat normalization is not compatible with specified chromosomes !!!")
          return()
        }
        
        Rdat.median<-apply(norm.data,1,function(x){
          if(f.qfilter[2]!=0){                                                         # if f.qfilter is on, filter out numbers from every row 
            x<-x[!(x<(quantile(x, probs=f.qfilter[1], names=FALSE)/f.qfilter[2]))]     # that is < than (f.qfilter[1] quantile) /f.qfilter[2]
          }
          return(median(x))                          #calculate the median of every row (of a geneprofiler data)
        })                      
        
        Rdat.quantile<- quantile(Rdat.median, probs=f.quantile,names=FALSE) #calculate the f.quantile quantile of the row medians
        norm.factor<- f.normvalue/Rdat.quantile                             # normalization factor based on Rdat geneprofiler data
        
        if(f.qfilter[2]!=0){
          print(paste("rows are filtered by x> quantileA(x)/B, A=",f.qfilter[1], "B=", f.qfilter[2]))
          norm.text<-c(norm.text,paste("rows are filtered by x> quantileA(x)/B, A=",f.qfilter[1], "B=", f.qfilter[2]))
        }
        print(paste(f.quantile, "quantile median row-values:",Rdat.quantile))
        norm.text<-c(norm.text,paste(f.quantile, "quantile median row-values:",Rdat.quantile))
        norm.text<-c(norm.text,paste("Normalized to",f.normvalue))
        
        rm(norm.data)
      }
      ############################################################     Rdat  mean  normalization  ############################################################
      
      if(f.mode=="mean"){
        if(f.chr!=""){
          print("f.mode='mean' with Rdat normalization is not compatible with specified chromosomes !!!")
          return()
        }
        
        Rdat.mean<-apply(norm.data,1,mean)                            #calculate the mean of every row (of a geneprofiler data)
        Rdat.mean<- mean(Rdat.mean)                                   #calculate the mean of the row means
        norm.factor<- f.normvalue/Rdat.mean                             # normalization factor based on Rdat geneprofiler data
        
        print(paste("Mean of mean row-values:",Rdat.mean))
        norm.text<-c(norm.text,paste("mean of mean row-values:",Rdat.mean))
        norm.text<-c(norm.text,paste("Normalized to",f.normvalue))
        
        rm(norm.data)
      }
      ############################################################     Rdat geometric mean of the geom. mean normalization  ############################################################
      
      if(f.mode=="GmeanGmean"){
        if(f.chr!=""){
          print("f.mode='GmeanGmean' with Rdat normalization is not compatible with specified chromosomes !!!")
          return()
        }
        norm.data[norm.data<f.replace0to]<- f.replace0to     # also smaller than f.replace0to values will be replaced by f.replace0to
        norm.data<-log(norm.data,2)                          # changing to log2 scale
        norm.data[norm.data==-Inf]<-NA
        
        Rdat.mean<-apply(norm.data,1,mean, na.rm=TRUE)                            # calculate the geom mean of every row (of a geneprofiler data)
        Rdat.mean<- mean(Rdat.mean, na.rm=TRUE)                                   # calculate the mean of the row means
        Rdat.mean<- 2^Rdat.mean                                       # changing back to normal scale
        norm.factor<- f.normvalue/Rdat.mean                           # normalization factor based on Rdat geneprofiler data
        
        print(paste("Geometric mean of geometric mean row-values:",Rdat.mean))
        norm.text<-c(norm.text,paste("Geometric mean of geometric mean row-values:",Rdat.mean))
        norm.text<-c(norm.text,paste("Normalized to",f.normvalue))
        
        rm(norm.data)
      }
      
      ############################################################     Rdat geometric mean of the median normalization  ############################################################
      
      if(f.mode=="MedGmean"){
        if(f.chr!=""){
          print("f.mode='MedGmean' with Rdat normalization is not compatible with specified chromosomes !!!")
          return()
        }
        
        Rdat.median<-apply(norm.data,1,function(x){
          if(f.qfilter[2]!=0){                                                         # if f.qfilter is on, filter out numbers from every row 
            x<-x[!(x<(quantile(x, probs=f.qfilter[1], names=FALSE)/f.qfilter[2]))]     # that is < than (f.qfilter[1] quantile) /f.qfilter[2]
          }
          return(median(x))                          #calculate the median of every row (of a geneprofiler data)
        })    
        
        
        Rdat.median[Rdat.median<f.replace0to]<- f.replace0to     # also smaller than f.replace0to values will be replaced by f.replace0to
        
        Rdat.median<-log(Rdat.median,2)                          # changing to log2 scale
        Rdat.median[Rdat.median==-Inf]<-NA
        
        Rdat.mean<- mean(Rdat.median, na.rm=TRUE)                                   # calculate the geometric mean of the row medians
        Rdat.mean<- 2^Rdat.mean                                       # changing back to normal scale
        norm.factor<- f.normvalue/Rdat.mean                           # normalization factor based on Rdat geneprofiler data
        
        print(paste("Geometric mean of median row-values:",Rdat.mean))
        norm.text<-c(norm.text,paste("Geometric mean of median row-values:",Rdat.mean))
        norm.text<-c(norm.text,paste("Normalized to",f.normvalue))
        
        rm(norm.data)
      }
      
    }
    
  }    # end of Normalization file section!
  
  ############################################################ ############################################################ ############################################################ 
  
  rm(Rle.data)
  
  Rle.data1[Rle.data1<f.replace0to]<- f.replace0to     #  smaller than f.replace0to values will be replaced by the f.replace0to value (and later multiplied by the norm factor)
  new.data1<-Rle.data1*norm.factor
  export(new.data1,paste("N.",f.filename1,sep=""))
  write(unlist(norm.text),file=paste("N.",f.filename1,".txt",sep=""))
  rm(new.data1,Rle.data1)
  
  if(f.filename2!=""){
    Rle.data2[Rle.data2<f.replace0to]<- f.replace0to     #  smaller than f.replace0to values will be replaced by the f.replace0to value (and later multiplied by the norm factor)
    new.data2 <-Rle.data2*norm.factor
    export(new.data2,paste("N.",f.filename2,sep=""))
    write(unlist(norm.text),file=paste("N.",f.filename2,".txt",sep=""))
  }
}


library("ShortRead")
library("rtracklayer")
library(stats)
library("IRanges")


f.bw.log2 <- function(f.filename, f.savename,f.replace0to=0.5) {
  
  # f.filename - BigWig file to change to log2 scale
  # f.savename - log2.bw will be added to the end
  # f.replace0to  - 0 will be replaced by this number
  
  GR.data <- import(f.filename)
  mcols(GR.data)$score[mcols(GR.data)$score==0]<-f.replace0to
  mcols(GR.data)$score<-log(mcols(GR.data)$score,2)
  export(GR.data,paste(f.savename,".log2.bw",sep=""))
}

f.bw.ratio <- function(f.filename1,f.filename2, f.savename,f.replace0to=0.5, f.log2=TRUE) {
  
  # f.filename1 - BigWig file - the numerator (top part of the fraction) 
  # f.filename2 - BigWig file - the denominator (bottom part of the fraction) - used as a reference, e.g. WT or WCE, etc.
  # f.savename  - R will be added at the beginning, .bw at the end, if log2=TRUE .log2.bw at the end
  # f.replace0to  - anything that is smaller or equal to this number will be replaced by  this number in both files 
  # f.log2 = TRUE/FALSE  - the result will be converted to log2 scale
  
  
  Rle.data1 <- import(f.filename1,as="RleList")
  Rle.data2 <- import(f.filename2,as="RleList")
  
  suppressWarnings(Rle.data1[Rle.data1<=f.replace0to]<-f.replace0to)
  suppressWarnings(Rle.data2[Rle.data2<=f.replace0to]<-f.replace0to)
  
  result.data<- Rle.data1/Rle.data2
  rm(Rle.data1, Rle.data2)
  
  if(f.log2){
    result.data<-log(result.data,2)
    f.savename<-paste("log2.",f.savename,sep="")
  }
  
  export(result.data,paste("R.",f.savename,sep=""))
}

f.bw.average <- function(f.filename1,f.filename2, f.savename) {
  
  # f.filename1 - BigWig file - the numerator (top part of the fraction) 
  # f.filename2 - BigWig file - the denominator (bottom part of the fraction) - used as a reference, e.g. WT or WCE, etc.
  # f.savename  - Av will be added at the beginning, .bw at the end, if log2=TRUE .log2.bw at the end
  
  
  Rle.data1 <- import(f.filename1,as="RleList")
  Rle.data2 <- import(f.filename2,as="RleList")
  
  
  result.data<- (Rle.data1+Rle.data2)/2
  
  rm(Rle.data1, Rle.data2)
  
  export(result.data,paste("Av",f.savename,sep=""))
}


f.bw.subtract <- function(f.filename1,f.filename2, f.savename) {
  
  # this function makes sense to use only on log2 data!!
  # f.filename1 - BigWig file - result will be f.filename1-f.filename2
  # f.filename2 - BigWig file - result will be f.filename1-f.filename2
  # f.savename  - S will be added at the beginning, .bw at the end,
  
  
  Rle.data1 <- import(f.filename1,as="RleList")
  Rle.data2 <- import(f.filename2,as="RleList")
  
  result.data<- Rle.data1-Rle.data2
  rm(Rle.data1, Rle.data2)
  
  export(result.data,paste("S",f.savename,sep=""))
}

f.bw.add <- function(f.filename1,f.filename2, f.savename) {
  
  # f.filename1 - BigWig file - result will be f.filename1+f.filename2
  # f.filename2 - BigWig file - result will be f.filename1+f.filename2
  # f.savename  - A will be added at the beginning, .bw at the end,
  
  
  Rle.data1 <- import(f.filename1,as="RleList")
  Rle.data2 <- import(f.filename2,as="RleList")
  
  result.data<- Rle.data1+Rle.data2
  rm(Rle.data1, Rle.data2)
  
  export(result.data,paste("A",f.savename,sep=""))
}


library("lattice")
library("latticeExtra")
library("caTools")

profile.plotter<-function(f.filename,f.mode="mean",f.quantile=0.5,f.name="NA",f.printmatrix=integer(), f.printmatrix2=integer(), f.xlim=c(0,0), 
                          f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-30,0,30),f.panels=character(0), f.panel.middle=0, f.runmean.segment=TRUE,f.runmean=c(1,1,1),f.oldfilename="",
                          f.finalplot=FALSE, f.log2=FALSE,f.replace0to=0.5,f.convertbackto10=FALSE, f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="",f.mod.annotation="",...){
  
  # f.filename - Rdat filename, created by geneprofiler; if f.filename="" - it will only print the f.printmatrix
  # f.mode="quantile" or "mean",f.quantile=0.5, if f.mode = Quantile, the f.quantile will be generated of every column, otherwise the mean will be calculated (default)
  # if f.log2 = TRUE, data will be converted to log2 scale and geometric averge will be calculated; 0 reads (and < than f.replace0to) will be adjusted to f.replace0to
  # if f.convertbackto10 == TRUE - if f.log2 was TRUE, it willconvert back to decimal before printing (geom average but decimal printing)
  # f.name - name of the series on the graph and in data frame returned from the function
  # f.printmatrix - should be the same length as profile from f.filename, except if adjusted with f.xtrim, new data will be added as last column to f.printmatrix
  #               - if f.filename=="", only f.printmatrix will be printed
  # f.printmatrix2 - if not empty, it should have x values and as many other columns as panels, this will be colored by filling until 0 with grey color, f.xticks should also be given
  #                 - if less colums than panels, only those panels will have this shaded print that has data 
  # f.xtrim(min, max) will trim the data before printing, also returned dataframe will be trimed (x values outside of the range will be lost), if f.xtrim==(0,0), no triming
  # f.xlim(min, max) will set the x limits on the graph, no triming in the data frame
  # f.ylim(min, max) will set the y limits on the graph, no triming in the data frame
  # f.panels = c("Name1", "Name2", ...) - should give name to all columns (including new data), and it will be printed on the panel according the name
  #                                     - if 4 different names, 4 panel will be generated; 1st column(X not included) is the 1st name,  new data is the last name
  #                                     - if more columns than panels, panels will be recycled -e.g. if 2 panels, 1st col is 1st panel, 2nd is 2nd, 3rd col is 1st panel, etc.
  # f.panel.middle - if one number, all middle lines will be at this, if a list, then each panel should have a number where the middle line will be drawn
  # f.xlab; f.ylab  - x and y label
  # f.runmean.segment=TRUE - if TRUE, the segments (such as 5', gene, 3') will be separated and runmean will stop at segment boundaries
  # f.runmean - moving window average, the size of the window in bp - only for printing, the returned data will be the unaveraged, if segments=TRUE, vector with 3 numbers: u1e9pstream, bins, downstream
  # f.oldfilename="" - this has to be provided if f.filename is empty, but we need some of the variables from this file - any of them is ok
  # f.finalplot - TRUE/FALSE - if TRUE, finalpot is given back to the caller, otherwise the f.matrix
  # f.replace0to - used only of f.log2=TRUE, 0 values will be replaced by this number (default is 0.5 - -1 in log2 scale) 
  # f.legendcols - if 0, as many legend columns as the number of columns, otherwise the given number
  # f.middle - only important if no file is loaded, determines the position of the 2nd dotted line, if 0, no line
  # f.mod.annotation - if the Rdat files were generated on the server, you have to add here the same annotation filename/path on your computer!
  # f.xticks - c(a,b,c,d..) - the ticks on x axis 
  # f.filter - a filter file (variable name should be gr.filter) made for the annotation used in the profiler to generate the Rdat file
  
  # ... will be in the par.settings=simpleTheme(); e.g. col=c("green, "red") or other parameters(col, alpha, cex, pch, lty, lwd, font, fill, 
  #                                                                                              border,col.points, col.line, alpha.points, alpha.line)
  # ... e.g. col=c("green","red"); lwd=5
  
  ####
  #f.filter="~/Documents/annotation/S.pombe.EF2/filters/EF2.mRNA.NoHetChr.filter"
  #f.filename<-"exp.hrp1dhrp3d.c043.profileAS.Rdat"
  #f.mode="mean"
  #f.quantile=0.5
  #f.name="NA"
  #f.printmatrix=integer()
  #f.printmatrix2=integer()
  #f.xlim=c(0,0)
  #f.ylim=c(0,0)
  #f.xtrim=c(0,0)
  #f.xticks=c(-30,0,30)
  #f.panels=character(0)
  #f.runmean.segment=TRUE
  #f.runmean=c(1,1,1)
  #f.log2=FALSE
  #f.replace0to=0.5
  #f.xlab="bp"
  #f.ylab="Relative occupancy (log2)"
  #f.legendcols=0
  #f.middle=0
  #f.filter=""
  #f.mod.annotation=""
  
  ####  
  
  
  if(f.filename !=""){
    load(f.filename)    					# data loaded into profile.rle; annotation file name is in f.annotation
    # profile.rle   - profiles for all features in ann.ger (except f.excluded) in rle format or if f.savediskspace= FALSE in matrix format
    # f.annotation  - name and path of annotation file used in the profile generation; will load ann.gr
    # f.upstream, f.middle, f.downstream - upstream, middle and downstream sequence length in the profile in bp.
    # f.excluded - integer array  of the excluded features in f.annotation file (ann.gr)
    # f.savediskspace - if TRUE, profiles in Rle format, if FALSE, profiles in matrix 
    
    f.middle<-f.binnumber*f.binlength                     # the second dotted line should be f.binnumber*f.binlength
    
    #if(f.mod.annotation==""){
    #load(f.annotation)            # annotation GRanges in ann.gr
    # ann.gr - annotation file used in the profile generation
    #}else{                                                        # if Rdat files were done on an other computer, you have to load the annotation
    #load(f.mod.annotation)
    #}
    
    ################
    # filtering here
    ################
    
    if(f.filter!=""){
      load(f.filter)     # gr.filter is loaded length is the original length of the annotation file
      if(length(f.excluded)>0){
        gr.filter<-gr.filter[-f.excluded]         # excluded values by the profiler should be excluded from the filter 
      }
      profile.rle<-profile.rle[gr.filter,]         # profile.rle is filtered 
    }
    
    print(nrow(profile.rle))
    
    ################
    
    ################
    
    if(f.savediskspace==TRUE){
      profile.num<-matrix(unlist(mclapply(1:length(profile.rle),function(x) as.numeric(profile.rle[[x]]), mc.cores=2*f.cores))
                          ,nrow=length(profile.rle), ncol=f.upstream+f.downstream, byrow=TRUE)  #converting Rle list to matrix of numerical values
      rm(profile.rle)
    }else{                                   # if f.savediskspace = FALSE, profile.rle is already converted to integer matrix
      profile.num<-profile.rle
      rm(profile.rle)
    }
    
    if(f.log2==TRUE){                       # if f.log2 = TRUE, data converted to log2 scale, geometric averge is calculated; 0 reads adjusted to f.replace0to reads
      profile.num[profile.num<f.replace0to]<- f.replace0to     # also smaller than f.replace0to values will be replaced by f.replace0to
      profile.num<-log(profile.num,2)
      profile.num[profile.num==-Inf]<-NA
    }
    
    if(f.mode=="quantile"){
      profile=apply(profile.num,2,function(x) quantile(x,probs=f.quantile, names=FALSE))    #calculate the median of every column (of a geneprofiler data)    
    }else{                         # mean
      profile<-colMeans(profile.num, na.rm=TRUE)
    }
    
    if(f.convertbackto10==TRUE&f.log2==TRUE){   #if we changed it to log2 scale but want to print normal scale (e,g, geom average but print linear scale)
      profile<-2^profile
    }  
    
    if(f.middle==0|f.binlength==1){                        # if no "average gene" feature (only 5 or 3 prime endprofiler) OR binlength=1
      X<-(-f.upstream:(f.middle+f.downstream-1))
    }else{                                  # if bins are used and length !=1
      if(f.upstream>0){
        X<-c(-f.upstream:-1,seq(from=f.binlength/2, to=f.middle-f.binlength/2, by=f.binlength), (f.middle):(f.middle+f.downstream-1))
      }else{
        X<-c(seq(from=f.binlength/2, to=f.middle-f.binlength/2, by=f.binlength), (f.middle):(f.middle+f.downstream-1))
      }
    }
    
    
    avg.data<-as.data.frame(cbind(X,profile))
    colnames(avg.data)<-c("X",f.name)
    if(sum(abs(f.xtrim))>0){               # if f.xtrim()is not (0,0) --> trim the data
      avg.data<-avg.data[avg.data[,1]>=f.xtrim[1]&avg.data[,1]<=f.xtrim[2],]
    } 
  }
  
  ###########################################################################################################################################
  ###############                                if f.printmatrix brings other data series                                    ###############
  ###########################################################################################################################################
  
  if(length(f.printmatrix)>0){                   # if f.printmatrix brings other data series  
    if(sum(abs(f.xtrim))>0){                     # if f.xtrim()is not (0,0) --> trim f.printmatrix
      f.printmatrix<-f.printmatrix[f.printmatrix[,1]>=f.xtrim[1]&f.printmatrix[,1]<=f.xtrim[2],]
    }
    if(f.filename==""){
      avg.data<-f.printmatrix                   # if f.filename=="" - no new data, print only f.printmatrix
    }else{                                                      # otherwise 
      if(nrow(avg.data)!=nrow(f.printmatrix)){                  # check if f.printmatrix length is the same as new profile
        print("X coordinates in new profile don't match X coordintaes in f.printmatrix !!!")
        return(f.printmatrix)
      }else{                                                    # append new profile to f.printmatrix as last column
        while(sum(colnames(f.printmatrix)==f.name)>0){          # if f.printmatrix has already column with f.name, append ".new"
          f.name<-paste(f.name,".new",sep="" )
        }
        avg.data<-cbind(f.printmatrix,avg.data[,2])
        colnames(avg.data)<-c(colnames(f.printmatrix),f.name)
      }
    } 
  }
  # new data plus print.matrix in avg.data dataframe, 1st column is X
  f.printmatrix<-avg.data                                     # f.printmatrix will be handed back to caller
  
  if(f.filename ==""){               # if f.filename is empty, we neeed to load one compatible rle file to get f.upsteam, and other variables
    load(f.oldfilename)    					# data loaded into profile.rle; annotation file name is in f.annotation
  }
  
  if(f.runmean.segment){              # if runmean should stop at segment boundaries
    if(sum(f.runmean>1)>0){
      if(f.upstream>0){
        avg.data[1:f.upstream,]<-as.data.frame(cbind(avg.data[1:f.upstream,1],sapply(2:ncol(avg.data), function(x) runmean(avg.data[1:f.upstream,x],f.runmean[1]))))     #moving window average
        if(length(f.printmatrix2)>0){
          f.printmatrix2[1:f.upstream,2]<- runmean(f.printmatrix2[1:f.upstream,2],f.runmean[1])     #moving window average
        }
      }
      if(f.binnumber>0){
        avg.data[(f.upstream+1):(f.upstream+f.binnumber),]<-as.data.frame(cbind(avg.data[(f.upstream+1):(f.upstream+f.binnumber),1],
                                                                                sapply(2:ncol(avg.data), function(x) runmean(avg.data[(f.upstream+1):(f.upstream+f.binnumber),x],f.runmean[2]))))     #moving window average
        if(length(f.printmatrix2)>0){
          f.printmatrix2[(f.upstream+1):(f.upstream+f.binnumber),2]<- runmean(f.printmatrix2[(f.upstream+1):(f.upstream+f.binnumber),2],f.runmean[2])     #moving window average
        }
      }
      if(f.downstream>0){
        avg.data[(f.upstream+f.binnumber+1):(f.upstream+f.binnumber+f.downstream),]<-as.data.frame(cbind(avg.data[(f.upstream+f.binnumber+1):(f.upstream+f.binnumber+f.downstream),1],
                                                                                                         sapply(2:ncol(avg.data), function(x) runmean(avg.data[(f.upstream+f.binnumber+1):(f.upstream+f.binnumber+f.downstream),x],f.runmean[3]))))     #moving window average
        if(length(f.printmatrix2)>0){
          f.printmatrix2[(f.upstream+f.binnumber+1):(f.upstream+f.binnumber+f.downstream),2]<- runmean(f.printmatrix2[(f.upstream+f.binnumber+1):(f.upstream+f.binnumber+f.downstream),2],f.runmean[3])     #moving window average
        }
      }
      
      colnames(avg.data)<-colnames(f.printmatrix)
    }
  }else{                            # if runmean should stop go over segment boundaries
    if(f.runmean[1]>1){
      avg.data<-as.data.frame(cbind(avg.data[,1],sapply(2:ncol(avg.data), function(x) runmean(avg.data[,x],f.runmean[1]))))     #moving window average
      colnames(avg.data)<-colnames(f.printmatrix)
      if(length(f.printmatrix2)>0){
        i<-2
        while(i<=length(f.printmatrix2)){
          f.printmatrix2[,i]<- runmean(f.printmatrix2[,i],f.runmean[1])     #moving window average
          i<-i+1
        }
        
      }
    }
  }
  
  
  ###########################################################################################################################################
  ###############                                            Preparation for printing                                         ###############
  ###########################################################################################################################################
  
  i<-2
  iPanel<-1
  X<-integer(0)
  Reads<-integer(0)
  Panel<-integer(0)
  Experiment<-integer(0)
  
  while(i<=ncol(avg.data)){                                           # converting horizontal data_frame to vertical df 
    X<-c(X,avg.data$X)                                                # Experiment - Panel
    Reads<-c(Reads,avg.data[,i])
    if(length(f.panels)==0){                                           # Panels
      Panel<-c(Panel,rep("Composite plot",nrow(avg.data)))
    }else{
      Panel<-c(Panel,rep(f.panels[iPanel],nrow(avg.data)))
    }
    Experiment<-c(Experiment,rep(colnames(avg.data)[i],nrow(avg.data)))
    i<-i+1
    iPanel<-iPanel+1
    if(iPanel>length(f.panels)){                                   # dataset is distributed between the panels, every new column is to a new panel
      iPanel<-1
    }
  }
  
  print.data<- as.data.frame(cbind(X,Reads,Panel,Experiment))
  print.data$X<-as.numeric(as.character(print.data$X))
  print.data$Reads<-as.numeric(as.character(print.data$Reads))
  
  
  if(f.legendcols==0){            # if f.legendcols not specified, the number of data columns
    f.legendcols<-i-2
  }
  
  ###########################################################################################################################################
  ###############                                            P R I N T I N G                                                  ###############
  ###########################################################################################################################################
  
  if(length(f.printmatrix2)==0){
    final.plot<-xyplot(Reads~X | Panel, data=print.data, groups=Experiment[,drop=TRUE], type="l",par.settings=simpleTheme(...),
                       auto.key=list(columns=f.legendcols, points=FALSE,lines=TRUE),xlab=f.xlab,ylab=f.ylab,
                       scales=list(x=list(relation="same",at=f.xticks),y="same"),
                       panel=function(x,y,...){
                         #panel.xyarea(genemed.matrix2[,1],genemed.matrix2[,2], origin=0,col="lightgrey",col.line="lightgrey",border="lightgrey",lwd=8)
                         panel.xyplot(x,y,...)
                         
                         if(length(f.panel.middle)>1){                                          # if multiple values for f.panel.middle, then use them one by one for the panels
                           if(panel.number()<=length(f.panel.middle)){
                             panel.lines(c(f.panel.middle[panel.number()],f.panel.middle[panel.number()]),c(-100,100), col="black", lwd=2, lty=2)
                           }
                         }else{
                           panel.lines(c(f.panel.middle,f.panel.middle),c(-100,100), col="black", lwd=2, lty=2)    # if only 1 value for f.panel.middle, then use  that for all panels
                         }
                         
                         if(f.middle!=0){
                           panel.lines(c(f.middle,f.middle),c(-100,100), col="black", lwd=2, lty=3)
                         } 
                         
                         
                         panel.lines(c(-10000,10000),c(0,0), col="black", lwd=2, lty=2)        # horizontal line at 0
                       })
  }else{
    final.plot<-xyplot(Reads~X | Panel, data=print.data, groups=Experiment[,drop=TRUE], type="l",par.settings=simpleTheme(...),
                       auto.key=list(columns=f.legendcols, points=FALSE,lines=TRUE),xlab=f.xlab,ylab=f.ylab,
                       scales=list(x=list(relation="same",at=f.xticks),y="same"),
                       
                       panel=function(x,y,...){
                         
                         if(panel.number()<=(ncol(f.printmatrix2)-1)){
                           panel.xyarea(f.printmatrix2[,1],f.printmatrix2[,(panel.number()+1)], origin=0,col="lightgrey",col.line="lightgrey",border="lightgrey",lwd=8)
                         }
                         
                         panel.xyplot(x,y,...)
                         
                         if(length(f.panel.middle)>1){                                          # if multiple values for f.panel.middle, then use them one by one for the panels
                           if(panel.number()<=length(f.panel.middle)){
                             panel.lines(c(f.panel.middle[panel.number()],f.panel.middle[panel.number()]),c(-100,100), col="black", lwd=2, lty=2)
                           }
                         }else{
                           panel.lines(c(f.panel.middle,f.panel.middle),c(-100,100), col="black", lwd=2, lty=2)    # if only 1 value for f.panel.middle, then use  that for all panels
                         }
                         
                         if(length(f.panel.middle)>1){                                          # if multiple values for f.panel.middle, then use them one by one for the panels
                           if(panel.number()<=length(f.panel.middle)){
                             panel.lines(c(f.panel.middle[panel.number()],f.panel.middle[panel.number()]),c(-100,100), col="black", lwd=2, lty=2)
                           }
                         }else{
                           panel.lines(c(f.panel.middle,f.panel.middle),c(-100,100), col="black", lwd=2, lty=2)    # if only 1 value for f.panel.middle, then use  that for all panels
                         }
                         
                         if(f.middle!=0){
                           panel.lines(c(f.middle,f.middle),c(-100,100), col="black", lwd=2, lty=3)
                         } 
                         
                       })
  }
  
  
  
  if(sum(abs(f.ylim))>0){                               # if ylim was specified
    final.plot<-update(final.plot, ylim=f.ylim, evaluate=FALSE)
  }
  if(sum(abs(f.xlim))>0){                               # if xlim was specified
    if(sum(abs(f.ylim))>0){                               # if ylim was specified
      final.plot<-update(final.plot, ylim=f.ylim, xlim=f.xlim, evaluate=FALSE)
    }else{
      final.plot<-update(final.plot, xlim=f.xlim, evaluate=FALSE)
    }
  }
  

  print(final.plot)
  
 
  
  ###########################################################################################################################################
  ###############                                              R E T U R N                                                    ###############
  ########################################################################################################################################### 
  
  
  if(f.finalplot==TRUE){
    return(final.plot)
  }else{
    return(f.printmatrix)                               # f.printmatrix has previous f.printmatrix with new trim and new profile in the last column 
  }
}  


peaks<-function(series,span=3, ties.method = "first")
{
  if((span <- as.integer(span)) %% 2 != 1) stop("'span' must be odd")
  z <- embed(series, span)
  s <- span%/%2
  v <- max.col(z, ties.method=ties.method) == 1 + s
  pad <- rep(FALSE, s)
  result <- c(pad, v, pad)
  result
}

