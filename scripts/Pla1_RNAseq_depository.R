# processing RNAseq data

# Generating profile.Rdat files from ORFs for normalization - Gene profiler

file.list1<-dir(pattern="*.f.bw$")
file.list2<-dir(pattern="*.r.bw$")
file.list3<-substr(file.list1, start=1, stop=nchar(file.list1)-5)

f.annotation="~/S.pombe.EF2/EF2.mRNA.Rdat"

lapply(1:length(file.list1), function(x) f.gene.profiler(file.list1[[x]], file.list2[[x]], file.list3[[x]],f.annotation,f.sense=TRUE,f.binnumber=25,
                                                         f.binlength=1,f.upstream=0, f.downstream=0,f.cores=4))

# Normalization using mRNA profileS.Rdat files 

file.list1<-dir(pattern="*.f.bw$")
file.list2<-dir(pattern="*.r.bw$")
file.list3<-dir(pattern="*.profileS.Rdat$")

#GmeanGmean500
mclapply(1:length(file.list1), function(x) f.BW.norm(file.list1[[x]], f.filename2=file.list2[[x]],f.mode="GmeanGmean",f.normvalue=500,
                                       f.chr="",f.normfile=file.list3[[x]],f.quantile=0.5,f.qfilter=c(0.5,10),f.replace0to=0.5),mc.cores=8)


setwd("~/GmeanGmean500")
# creating avareage bw files from Rep1 and Rep2
file.list1<-dir(pattern="*.f.bw$")

f.bw.average("N.sp_st3230_RED1RE_Rep1_polyA_Illumina2_ERCC_RNAseq.f.bw","N.sp_st3230_RED1RE_Rep2_polyA_Illumina2_ERCC_RNAseq.f.bw","N.sp_st3230_RED1RE_Rep1_Rep2_average.f.bw")
f.bw.average("N.sp_st3230_RED1RE_Rep1_polyA_Illumina2_ERCC_RNAseq.r.bw","N.sp_st3230_RED1RE_Rep2_polyA_Illumina2_ERCC_RNAseq.r.bw","N.sp_st3230_RED1RE_Rep1_Rep2_average.r.bw")

f.bw.average("N.sp_st3232_RED1DI_Rep1_polyA_Illumina2_ERCC_RNAseq.f.bw","N.sp_st3232_RED1DI_Rep2_polyA_Illumina2_ERCC_RNAseq.f.bw","N.sp_st3232_RED1DI_Rep1_Rep2_average.f.bw")
f.bw.average("N.sp_st3232_RED1DI_Rep1_polyA_Illumina2_ERCC_RNAseq.r.bw","N.sp_st3232_RED1DI_Rep2_polyA_Illumina2_ERCC_RNAseq.r.bw","N.sp_st3232_RED1DI_Rep1_Rep2_average.r.bw")

f.bw.average("N.sp_st3436_RED1DEL_Rep1_polyA_Illumina2_ERCC_RNAseq.f.bw","N.sp_st3436_RED1DEL_Rep2_polyA_Illumina2_ERCC_RNAseq.f.bw","N.sp_st3436_RED1DEL_Rep1_Rep2_average.f.bw")
f.bw.average("N.sp_st3436_RED1DEL_Rep1_polyA_Illumina2_ERCC_RNAseq.r.bw","N.sp_st3436_RED1DEL_Rep2_polyA_Illumina2_ERCC_RNAseq.r.bw","N.sp_st3436_RED1DEL_Rep1_Rep2_average.r.bw")



# GeneProfiler
file.list1<-dir(pattern="*.f.bw$")
file.list2<-dir(pattern="*.r.bw$")
file.list3<-substr(file.list1, start=1, stop=nchar(file.list1)-5)

file.list1<-file.list1[1:4]   #only average files
file.list2<-file.list2[1:4]
file.list3<-file.list3[1:4]

f.annotation="~/S.pombe.EF2/EF2.mRNA.Rdat"

lapply(1:length(file.list1), function(x) f.gene.profiler(file.list1[[x]], file.list2[[x]], file.list3[[x]],f.annotation,f.sense=TRUE,f.binnumber=60,
                                                         f.binlength=20,f.upstream=500, f.downstream=500,f.cores=8))

lapply(1:length(file.list1), function(x) f.gene.profiler(file.list1[[x]], file.list2[[x]], file.list3[[x]],f.annotation,f.sense=FALSE,f.binnumber=60,
                                                         f.binlength=20,f.upstream=500, f.downstream=500,f.cores=8))

#[1] "Lengths anomalie! - Excluded features: 3"
#[1] "SPAC750.08c" "tlh1"        "SPBC1348.15"

# Profile plotter - Average genes

file.list1<-dir(pattern="*average.profileS.Rdat$")
file.list2<-dir(pattern="*average.profileAS.Rdat$")
print.S.WT=integer()
print.AS.WT=integer()
WT.order=1

print.S.WT<-profile.plotter(file.list1[WT.order],f.mode="mean",f.name="WT",f.printmatrix=print.S.WT, f.printmatrix2=integer(), f.xlim=c(0,0), 
                            f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-500,0,1200,1700),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.1, 
                            f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_ASs.filter",f.mod.annotation="")

print.AS.WT<-profile.plotter(file.list2[WT.order],f.mode="mean",f.name="WT",f.printmatrix=print.AS.WT, f.printmatrix2=integer(), f.xlim=c(0,0), 
                             f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-500,0,1200,1700),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.1, 
                             f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_ASs.filter",f.mod.annotation="")

print.S=integer()
ser.order=c(1,2,4,3)

ser.name<-c("1_WT","2_DI","4_KO","3_DEL")

for (x in ser.order){
  print.S<-profile.plotter(file.list1[x],f.mode="mean",f.name=ser.name[x],f.printmatrix=print.S, f.printmatrix2=print.S.WT, f.xlim=c(0,0), 
                           f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-500,0,1200,1700),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.1, 
                           f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_ASs.filter",f.mod.annotation="")
}


aaa<-profile.plotter("",f.mode="mean",f.name="",f.printmatrix=print.S, f.printmatrix2=print.S.WT, f.xlim=c(0,0), 
                     f.ylim=c(0,0),f.xtrim=c(0,0),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.oldfilename=file.list1[1],f.log2=TRUE,f.replace0to=0.1, 
                     f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=4, f.xticks=c(-500,0,1200,1700),f.middle=1200,f.filter="EF2.mRNA.red1KO_affected_ASs.filter",f.mod.annotation=""
                     ,col=c("black","orange","blue","red"),lwd=c(6,6,6,6))

print.AS=integer()

for (x in ser.order){
  print.AS<-profile.plotter(file.list2[x],f.mode="mean",f.name=ser.name[x],f.printmatrix=print.AS, f.printmatrix2=print.AS.WT, f.xlim=c(0,0), 
                            f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-500,0,1200,1700),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.1, 
                            f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_ASs.filter",f.mod.annotation="")
}

aaa<-profile.plotter("",f.mode="mean",f.name="",f.printmatrix=print.AS, f.printmatrix2=print.AS.WT, f.xlim=c(0,0), 
                     f.ylim=c(0,0),f.xtrim=c(0,0),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.oldfilename=file.list2[1],f.log2=TRUE,f.replace0to=0.1, 
                     f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=4, f.xticks=c(-500,0,1200,1700),f.middle=1200,f.filter="EF2.mRNA.red1KO_affected_ASs.filter",f.mod.annotation=""
                     ,col=c("black","orange","blue","red"),lwd=c(6,6,6,6))

############################################### +/- 500 gene start and gene end #######################################
############################################### +/- 500 gene start and gene end #######################################
############################################### +/- 500 gene start and gene end #######################################
############################################### +/- 500 gene start and gene end #######################################
############################################### +/- 500 gene start and gene end #######################################
############################################### +/- 500 gene start and gene end #######################################


setwd("~/MedGmean500")

# EndProfiler
file.list1<-dir(pattern="*.f.bw$")
file.list2<-dir(pattern="*.r.bw$")
file.list3<-substr(file.list1, start=1, stop=nchar(file.list1)-5)
f.annotation="~/S.pombe.EF2/EF2.mRNA.Rdat"

lapply(1:length(file.list1), function(x) f.end.profiler(file.list1[[x]], file.list2[[x]], file.list3[[x]],f.annotation,f.sense=TRUE,f.end="5end",
                                                        f.upstream=500, f.downstream=750))

lapply(1:length(file.list1), function(x) f.end.profiler(file.list1[[x]], file.list2[[x]], file.list3[[x]],f.annotation,f.sense=TRUE,f.end="3end",
                                                        f.upstream=750, f.downstream=500))

lapply(1:length(file.list1), function(x) f.end.profiler(file.list1[[x]], file.list2[[x]], file.list3[[x]],f.annotation,f.sense=FALSE,f.end="5end",
                                                        f.upstream=500, f.downstream=750))

lapply(1:length(file.list1), function(x) f.end.profiler(file.list1[[x]], file.list2[[x]], file.list3[[x]],f.annotation,f.sense=FALSE,f.end="3end",
                                                        f.upstream=750, f.downstream=500))
########## Sense end profiles #############

file.list1<-dir(pattern="*.5profileS.Rdat$")
file.list2<-dir(pattern="*.3profileS.Rdat$")

WT.order=c(1,2)                     #Re1 and Re2 - we need the average of these 2

print.S.WTP1=integer()               #WT profile for panel 1 (5' profiles)
print.S.WTP2=integer()               #WT profile for panel 2 (3' profiles)

for (x in WT.order){
  print.S.WTP1<-profile.plotter(file.list1[x],f.mode="mean",f.name="WT",f.printmatrix=print.S.WTP1, f.printmatrix2=integer(), f.xlim=c(0,0), 
                            f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-500,0,750),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.5, 
                            f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation="")


  print.S.WTP2<-profile.plotter(file.list2[x],f.mode="mean",f.name="WT",f.printmatrix=print.S.WTP2, f.printmatrix2=integer(), f.xlim=c(0,0), 
                            f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-750,0,500),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.5, 
                            f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation="")
}

print.S.WTP1[,2]<-(print.S.WTP1[,2]+print.S.WTP1[,3])/2       # making avarage 5'WT profiles and dropping 2nd column
print.S.WTP1<-print.S.WTP1[,-3]
  
print.S.WTP2[,2]<-(print.S.WTP2[,2]+print.S.WTP2[,3])/2       # making avarage 3' WT profiles and dropping 2nd column
print.S.WTP2<-print.S.WTP2[,-3]
 
print.SP1=integer()
print.SP2=integer()

ser.order=c(5,1,2,3,4,6,7)

ser.name<-c("5_RE1","6_RE2","7_DI1","8_DI2","4_KO","9_DEL1", "A_DEL2")

for (x in ser.order){
  print.SP1<-profile.plotter(file.list1[x],f.mode="mean",f.name=ser.name[x],f.printmatrix=print.SP1, f.printmatrix2=print.S.WTP1, f.xlim=c(0,0), 
                           f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-500,0,750),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.5, 
                           f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation="")
}

for (x in ser.order){
  print.SP2<-profile.plotter(file.list2[x],f.mode="mean",f.name=ser.name[x],f.printmatrix=print.SP2, f.printmatrix2=print.S.WTP2, f.xlim=c(0,0), 
                           f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-750,0,500),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.5, 
                           f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation="")
}


print.SP1$RE<-(print.SP1[,3]+print.SP1[,4])/2
print.SP1$DI<-(print.SP1[,5]+print.SP1[,6])/2
print.SP1$DEL<-(print.SP1[,7]+print.SP1[,8])/2
names(print.SP1)[9]<-"1_RE"
names(print.SP1)[10]<-"2_DI"
names(print.SP1)[11]<-"3_DEL"

print.SP2$RE<-(print.SP2[,3]+print.SP2[,4])/2
print.SP2$DI<-(print.SP2[,5]+print.SP2[,6])/2
print.SP2$DEL<-(print.SP2[,7]+print.SP2[,8])/2
names(print.SP2)[9]<-"1_RE"
names(print.SP2)[10]<-"2_DI"
names(print.SP2)[11]<-"3_DEL"

print.S=print.SP1[1]
for(x in 2:length(print.SP1)){
  print.S<-cbind(print.S,print.SP1[x],print.SP2[x])
}

print.S.WT<-cbind(print.S.WTP1[,1:2], print.S.WTP2[,2])  
names(print.S.WT)[2]<-"WTP1"
names(print.S.WT)[3]<-"WTP2"

aaa<-profile.plotter("",f.mode="mean",f.name="",f.printmatrix=print.S, f.printmatrix2=print.S.WT, f.xlim=c(0,0), 
                     f.ylim=c(0,0),f.xtrim=c(0,0),f.panels=c("A_5_profile","B_3_profile"),f.panel.middle=c(0,250), f.runmean.segment=FALSE,f.runmean=c(30,30,30),f.oldfilename=file.list1[1],f.log2=TRUE,
                     f.replace0to=0.5, f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=5, f.xticks=c(-500,0,250,750),f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation=""
                     ,col=c("black","orange","blue","red","darkgrey","darkgrey","yellow","yellow","lightblue","lightblue"), lty=c(1,1,1,1,3,3,3,3,3,3),lwd=c(6,6,6,6,2,2,2,2,2,2))



########## Antisense end profiles #############

file.list1<-dir(pattern="*.5profileAS.Rdat$")
file.list2<-dir(pattern="*.3profileAS.Rdat$")

WT.order=c(1,2)                     #Re1 and Re2 - we need the average of these 2

print.AS.WTP1=integer()               #WT profile for panel 1 (5' profiles)
print.AS.WTP2=integer()               #WT profile for panel 2 (3' profiles)

for (x in WT.order){
  print.AS.WTP1<-profile.plotter(file.list1[x],f.mode="mean",f.name="WT",f.printmatrix=print.AS.WTP1, f.printmatrix2=integer(), f.xlim=c(0,0), 
                                f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-500,0,750),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.5, 
                                f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation="")
  
  
  print.AS.WTP2<-profile.plotter(file.list2[x],f.mode="mean",f.name="WT",f.printmatrix=print.AS.WTP2, f.printmatrix2=integer(), f.xlim=c(0,0), 
                                f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-750,0,500),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.5, 
                                f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation="")
}

print.AS.WTP1[,2]<-(print.AS.WTP1[,2]+print.AS.WTP1[,3])/2       # making avarage 5'WT profiles and dropping 2nd column
print.AS.WTP1<-print.AS.WTP1[,-3]

print.AS.WTP2[,2]<-(print.AS.WTP2[,2]+print.AS.WTP2[,3])/2       # making avarage 3' WT profiles and dropping 2nd column
print.AS.WTP2<-print.AS.WTP2[,-3]

print.ASP1=integer()
print.ASP2=integer()

ser.order=c(5,1,2,3,4,6,7)

ser.name<-c("5_RE1","6_RE2","7_DI1","8_DI2","4_KO","9_DEL1", "A_DEL2")

for (x in ser.order){
  print.ASP1<-profile.plotter(file.list1[x],f.mode="mean",f.name=ser.name[x],f.printmatrix=print.ASP1, f.printmatrix2=print.AS.WTP1, f.xlim=c(0,0), 
                             f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-500,0,750),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.5, 
                             f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation="")
}

for (x in ser.order){
  print.ASP2<-profile.plotter(file.list2[x],f.mode="mean",f.name=ser.name[x],f.printmatrix=print.ASP2, f.printmatrix2=print.AS.WTP2, f.xlim=c(0,0), 
                             f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-750,0,500),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.5, 
                             f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation="")
}


print.ASP1$RE<-(print.ASP1[,3]+print.ASP1[,4])/2
print.ASP1$DI<-(print.ASP1[,5]+print.ASP1[,6])/2
print.ASP1$DEL<-(print.ASP1[,7]+print.ASP1[,8])/2
names(print.ASP1)[9]<-"1_RE"
names(print.ASP1)[10]<-"2_DI"
names(print.ASP1)[11]<-"3_DEL"

print.ASP2$RE<-(print.ASP2[,3]+print.ASP2[,4])/2
print.ASP2$DI<-(print.ASP2[,5]+print.ASP2[,6])/2
print.ASP2$DEL<-(print.ASP2[,7]+print.ASP2[,8])/2
names(print.ASP2)[9]<-"1_RE"
names(print.ASP2)[10]<-"2_DI"
names(print.ASP2)[11]<-"3_DEL"

print.AS=print.ASP1[1]
for(x in 2:length(print.ASP1)){
  print.AS<-cbind(print.AS,print.ASP1[x],print.ASP2[x])
}

print.AS.WT<-cbind(print.AS.WTP1[,1:2], print.AS.WTP2[,2])  
names(print.AS.WT)[2]<-"WTP1"
names(print.AS.WT)[3]<-"WTP2"

aaa<-profile.plotter("",f.mode="mean",f.name="",f.printmatrix=print.AS, f.printmatrix2=print.AS.WT, f.xlim=c(0,0), 
                     f.ylim=c(0,0),f.xtrim=c(0,0),f.panels=c("A_5_profile","B_3_profile"),f.panel.middle=c(0,250), f.runmean.segment=FALSE,f.runmean=c(30,30,30),f.oldfilename=file.list1[1],f.log2=TRUE,
                     f.replace0to=0.5, f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=5, f.xticks=c(-500,0,250,750),f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation=""
                     ,col=c("black","orange","blue","red","darkgrey","darkgrey","yellow","yellow","lightblue","lightblue"), lty=c(1,1,1,1,3,3,3,3,3,3),lwd=c(6,6,6,6,2,2,2,2,2,2))

############################################### PROMT filter #######################################
############################################### PROMT filter #######################################
############################################### PROMT filter #######################################
############################################### PROMT filter #######################################
############################################### PROMT filter #######################################
############################################### PROMT filter #######################################
############################################### PROMT filter #######################################
############################################### PROMT filter #######################################

setwd("~/GmeanGmean500")

file.list1<-dir(pattern="*.5profileAS.Rdat$")

load(file.list1[1])
profile5.RE1<-profile.rle[,1:500]
load(file.list1[2])
profile5.RE2<-profile.rle[,1:500]

profile5.RE<-(profile5.RE1+profile5.RE2)/2          # average of RE1 and RE2

load(file.list1[5])
profile5.KO<-profile.rle[,1:500]

profile5.delta<-profile5.KO-profile5.RE           # delta between red1ko and RE average


delta.1<-rowMedians(profile5.delta[,1:250])       #rowMedians of delta 1-250 (upstream of gene on antisense direction) 
delta.2<-rowMedians(profile5.delta[,150:400])       #rowMedians of delta 150-400 (upstream of gene on antisense direction) 
delta.3<-rowMedians(profile5.delta[,250:500])       #rowMedians of delta 250-500 (upstream of gene on antisense direction) 

RE.1<-rowMedians(profile5.RE[,1:250])       #rowMedians of RE 1-250 (upstream of gene on antisense direction) 
RE.2<-rowMedians(profile5.RE[,150:400])       #rowMedians of RE 150-400 (upstream of gene on antisense direction) 
RE.3<-rowMedians(profile5.RE[,250:500])       #rowMedians of RE 250-500 (upstream of gene on antisense direction)

delta.1[delta.1<10]<-0                      #anything that is smaller than 10 differnece in median delta is set to 0
delta.2[delta.2<10]<-0
delta.3[delta.3<10]<-0

RE.1[RE.1<1]<-1                           #if RE median is smaller thna 1 then set to 1 (so that no division with 0)
RE.2[RE.2<1]<-1  
RE.3[RE.3<1]<-1  

ratio.1<-delta.1/RE.1                     # ratio between delat and WT expression level
ratio.2<-delta.1/RE.2
ratio.3<-delta.1/RE.3

ratio.1[ratio.1<0.5]<-0                   # minimum 50% more in KO than in WT (RE) delta is at least 0.5 ratio of WT expression AND more than 10
ratio.2[ratio.2<0.5]<-0 
ratio.3[ratio.3<0.5]<-0 

ratio.1.ann<-c(ratio.1[1:3261],0,ratio.1[3262:3524],0,ratio.1[3525:5116])       # we have to put back the 2 excluded gene (no 3262 and 3526 - these numbers are stored in f.excluded) as 0
ratio.2.ann<-c(ratio.2[1:3261],0,ratio.2[3262:3524],0,ratio.2[3525:5116])
ratio.3.ann<-c(ratio.3[1:3261],0,ratio.3[3262:3524],0,ratio.3[3525:5116])

gr.filter<-(ratio.1.ann>0|ratio.2.ann>0|ratio.3.ann>0)

save(gr.filter, file="EF2.mRNA.red1KO_affected_PROMTs.filter")  # 2091 genes



############################################### Introns #######################################
############################################### Introns #######################################
############################################### Introns #######################################
############################################### Introns #######################################
############################################### Introns #######################################
############################################### Introns #######################################
############################################### Introns #######################################
setwd("~/unnormalized")

# Generating profile.Rdat files from ORFs for normalization - Gene profiler

file.list1<-dir(pattern="*.f.bw$")
file.list2<-dir(pattern="*.r.bw$")
file.list3<-substr(file.list1, start=1, stop=nchar(file.list1)-5)
f.annotation="~/S.pombe.EF2/EF2.ORF.Rdat"

lapply(1:length(file.list1), function(x) f.gene.profiler(file.list1[[x]], file.list2[[x]], file.list3[[x]],f.annotation,f.sense=TRUE,f.binnumber=25,
                                                         f.binlength=1,f.upstream=0, f.downstream=0,f.cores=4))
# Normalization using mRNA profileS.Rdat files 

file.list1<-dir(pattern="*.f.bw$")
file.list2<-dir(pattern="*.r.bw$")
file.list3<-dir(pattern="*.profileS.Rdat$")

#MedGmean500
mclapply(1:length(file.list1), function(x) f.BW.norm(file.list1[[x]], f.filename2=file.list2[[x]],f.mode="MedGmean",f.normvalue=500,
                                       f.chr="",f.normfile=file.list3[[x]],f.quantile=0.5,f.qfilter=c(0.5,10),f.replace0to=0.5),mc.cores=8)


setwd("~/ORFMedGmean500")


# GeneProfiler
file.list1<-dir(pattern="*.f.bw$")
file.list2<-dir(pattern="*.r.bw$")
file.list3<-substr(file.list1, start=1, stop=nchar(file.list1)-5)
file.list3<-paste(file.list3,".intron",sep="")
f.annotation="~/Desktop/R/annotation/S.pombe.EF2/EF2.ORFintrons20to400.Rdat"

lapply(1:length(file.list1), function(x) f.gene.profiler(file.list1[[x]], file.list2[[x]], file.list3[[x]],f.annotation,f.sense=TRUE,f.binnumber=30,
                                                         f.binlength=1,f.upstream=20, f.downstream=20,f.cores=4))

# Profile plotter
setwd("~/ORFMedGmean500/Rep2")

file.list1<-dir(pattern="*.intron.profileS.Rdat$")
print.S.WT=integer()
ser.order=c(1)
ser.name<-c("WT")

print.S.WT<-profile.plotter(file.list1[ser.order[1]],f.mode="mean",f.name=ser.name[1],f.printmatrix=print.S.WT, f.printmatrix2=integer(), f.xlim=c(0,0), 
                            f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-20,0,30,50),f.panels=character(0), f.runmean.segment=TRUE,f.runmean=c(1,1,1),f.log2=TRUE,f.replace0to=0.5, 
                            f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="",f.mod.annotation="")

print.S=integer()
ser.order=c(1,2,4,3)
ser.name<-c("1. WT2","2. DI2", "4. KO2","3. DEL2")

for (x in ser.order){
  print.S<-profile.plotter(file.list1[x],f.mode="mean",f.name=ser.name[x],f.printmatrix=print.S, f.printmatrix2=print.S.WT, f.xlim=c(0,0), 
                           f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-20,0,30,50),f.panels=character(0), f.runmean.segment=TRUE,f.runmean=c(1,1,1),f.log2=TRUE,f.replace0to=0.5, 
                           f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="",f.mod.annotation="")
}


aaa<-profile.plotter("",f.mode="mean",f.name="",f.printmatrix=print.S, f.printmatrix2=print.S.WT, f.xlim=c(0,0), 
                     f.ylim=c(0,0),f.xtrim=c(0,0),f.panels=character(0), f.runmean.segment=TRUE,f.runmean=c(1,1,1),f.oldfilename=file.list1[1],f.log2=TRUE,f.replace0to=0.5,
                     f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.xticks=c(-20,0,29,49),f.middle=29,f.filter="",f.mod.annotation=""
                     ,col=c("black","orange","blue","red"), lwd=6)



###############################################  total RNAseq profiles #######################################
###############################################  total RNAseq profiles #######################################
###############################################  total RNAseq profiles #######################################
###############################################  total RNAseq profiles #######################################
###############################################  total RNAseq profiles #######################################
###############################################  total RNAseq profiles #######################################
###############################################  total RNAseq profiles #######################################
###############################################  total RNAseq profiles #######################################

setwd("~/totalRNA")

# Generating profile.Rdat files from ORFs for normalization - Gene profiler

file.list1<-dir(pattern="*.f.bw$")
file.list2<-dir(pattern="*.r.bw$")
file.list3<-substr(file.list1, start=1, stop=nchar(file.list1)-5)

f.annotation="~/Desktop/R/annotation/S.pombe.EF2/EF2.mRNA.Rdat"

lapply(1:length(file.list1), function(x) f.gene.profiler(file.list1[[x]], file.list2[[x]], file.list3[[x]],f.annotation,f.sense=TRUE,f.binnumber=25,
                                                         f.binlength=1,f.upstream=0, f.downstream=0,f.cores=4))

# Normalization using mRNA profileS.Rdat files 

file.list1<-dir(pattern="*.f.bw$")
file.list2<-dir(pattern="*.r.bw$")
file.list3<-dir(pattern="*.profileS.Rdat$")

#GmeanGmean500
mclapply(1:length(file.list1), function(x) f.BW.norm(file.list1[[x]], f.filename2=file.list2[[x]],f.mode="GmeanGmean",f.normvalue=500,
                                                     f.chr="",f.normfile=file.list3[[x]],f.quantile=0.5,f.qfilter=c(0.5,10),f.replace0to=0.5),mc.cores=8)


setwd("~/totalRNA/GmeanGmean500")


file.list1<-dir(pattern="*.f.bw$")
file.list2<-dir(pattern="*.r.bw$")
file.list3<-substr(file.list1, start=1, stop=nchar(file.list1)-5)
f.annotation="~/S.pombe.EF2/EF2.mRNA.Rdat"

lapply(1:length(file.list1), function(x) f.end.profiler(file.list1[[x]], file.list2[[x]], file.list3[[x]],f.annotation,f.sense=TRUE,f.end="5end",
                                                        f.upstream=500, f.downstream=750))

lapply(1:length(file.list1), function(x) f.end.profiler(file.list1[[x]], file.list2[[x]], file.list3[[x]],f.annotation,f.sense=TRUE,f.end="3end",
                                                        f.upstream=750, f.downstream=500))

lapply(1:length(file.list1), function(x) f.end.profiler(file.list1[[x]], file.list2[[x]], file.list3[[x]],f.annotation,f.sense=FALSE,f.end="5end",
                                                        f.upstream=500, f.downstream=750))

lapply(1:length(file.list1), function(x) f.end.profiler(file.list1[[x]], file.list2[[x]], file.list3[[x]],f.annotation,f.sense=FALSE,f.end="3end",
                                                        f.upstream=750, f.downstream=500))

########## Sense end profiles #############

file.list1<-dir(pattern="*.5profileS.Rdat$")
file.list2<-dir(pattern="*.3profileS.Rdat$")

WT.order=c(1,2)                     #Re1 and Re2 - we need the average of these 2

print.S.WTP1=integer()               #WT profile for panel 1 (5' profiles)
print.S.WTP2=integer()               #WT profile for panel 2 (3' profiles)

for (x in WT.order){
  print.S.WTP1<-profile.plotter(file.list1[x],f.mode="mean",f.name="WT",f.printmatrix=print.S.WTP1, f.printmatrix2=integer(), f.xlim=c(0,0), 
                                f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-500,0,750),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.5, 
                                f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation="")
  
  
  print.S.WTP2<-profile.plotter(file.list2[x],f.mode="mean",f.name="WT",f.printmatrix=print.S.WTP2, f.printmatrix2=integer(), f.xlim=c(0,0), 
                                f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-750,0,500),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.5, 
                                f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation="")
}

print.S.WTP1[,2]<-(print.S.WTP1[,2]+print.S.WTP1[,3])/2       # making avarage 5'WT profiles and dropping 2nd column
print.S.WTP1<-print.S.WTP1[,-3]

print.S.WTP2[,2]<-(print.S.WTP2[,2]+print.S.WTP2[,3])/2       # making avarage 3' WT profiles and dropping 2nd column
print.S.WTP2<-print.S.WTP2[,-3]

print.SP1=integer()
print.SP2=integer()

ser.order=c(5,6,1,2,3,4,7,8)

ser.name<-c("5_RE1","6_RE2","7_DI1","8_DI2","B_KO1","C_KO2","9_DEL1", "A_DEL2")

for (x in ser.order){
  print.SP1<-profile.plotter(file.list1[x],f.mode="mean",f.name=ser.name[x],f.printmatrix=print.SP1, f.printmatrix2=print.S.WTP1, f.xlim=c(0,0), 
                             f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-500,0,750),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.5, 
                             f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation="")
}

for (x in ser.order){
  print.SP2<-profile.plotter(file.list2[x],f.mode="mean",f.name=ser.name[x],f.printmatrix=print.SP2, f.printmatrix2=print.S.WTP2, f.xlim=c(0,0), 
                             f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-750,0,500),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.5, 
                             f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation="")
}

print.SP1$KO<-(print.SP1[,2]+print.SP1[,3])/2
print.SP1$RE<-(print.SP1[,4]+print.SP1[,5])/2
print.SP1$DI<-(print.SP1[,6]+print.SP1[,7])/2
print.SP1$DEL<-(print.SP1[,8]+print.SP1[,9])/2
names(print.SP1)[10]<-"4_KO"
names(print.SP1)[11]<-"1_RE"
names(print.SP1)[12]<-"2_DI"
names(print.SP1)[13]<-"3_DEL"

print.SP2$KO<-(print.SP2[,2]+print.SP2[,3])/2
print.SP2$RE<-(print.SP2[,4]+print.SP2[,5])/2
print.SP2$DI<-(print.SP2[,6]+print.SP2[,7])/2
print.SP2$DEL<-(print.SP2[,8]+print.SP2[,9])/2
names(print.SP2)[10]<-"4_KO"
names(print.SP2)[11]<-"1_RE"
names(print.SP2)[12]<-"2_DI"
names(print.SP2)[13]<-"3_DEL"

print.S=print.SP1[1]
for(x in 2:length(print.SP1)){
  print.S<-cbind(print.S,print.SP1[x],print.SP2[x])
}

print.S.WT<-cbind(print.S.WTP1[,1:2], print.S.WTP2[,2])  
names(print.S.WT)[2]<-"WTP1"
names(print.S.WT)[3]<-"WTP2"

aaa<-profile.plotter("",f.mode="mean",f.name="",f.printmatrix=print.S, f.printmatrix2=print.S.WT, f.xlim=c(0,0), 
                     f.ylim=c(0,0),f.xtrim=c(0,0),f.panels=c("A_5_profile","B_3_profile"),f.panel.middle=c(0,250), f.runmean.segment=FALSE,f.runmean=c(30,30,30),f.oldfilename=file.list1[1],f.log2=TRUE,
                     f.replace0to=0.5, f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=5, f.xticks=c(-500,0,250,750),f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation=""
                     ,col=c("black","orange","blue","red","darkgrey","darkgrey","yellow","yellow","lightblue","lightblue","pink","pink"), lty=c(1,1,1,1,3,3,3,3,3,3,3,3),lwd=c(6,6,6,6,2,2,2,2,2,2,2,2))



########## Antisense end profiles #############

file.list1<-dir(pattern="*.5profileAS.Rdat$")
file.list2<-dir(pattern="*.3profileAS.Rdat$")

WT.order=c(1,2)                     #Re1 and Re2 - we need the average of these 2

print.AS.WTP1=integer()               #WT profile for panel 1 (5' profiles)
print.AS.WTP2=integer()               #WT profile for panel 2 (3' profiles)

for (x in WT.order){
  print.AS.WTP1<-profile.plotter(file.list1[x],f.mode="mean",f.name="WT",f.printmatrix=print.AS.WTP1, f.printmatrix2=integer(), f.xlim=c(0,0), 
                                f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-500,0,750),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.5, 
                                f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation="")
  
  
  print.AS.WTP2<-profile.plotter(file.list2[x],f.mode="mean",f.name="WT",f.printmatrix=print.AS.WTP2, f.printmatrix2=integer(), f.xlim=c(0,0), 
                                f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-750,0,500),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.5, 
                                f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation="")
}

print.AS.WTP1[,2]<-(print.AS.WTP1[,2]+print.AS.WTP1[,3])/2       # making avarage 5'WT profiles and dropping 2nd column
print.AS.WTP1<-print.AS.WTP1[,-3]

print.AS.WTP2[,2]<-(print.AS.WTP2[,2]+print.AS.WTP2[,3])/2       # making avarage 3' WT profiles and dropping 2nd column
print.AS.WTP2<-print.AS.WTP2[,-3]

print.ASP1=integer()
print.ASP2=integer()

ser.order=c(5,6,1,2,3,4,7,8)

ser.name<-c("5_RE1","6_RE2","7_DI1","8_DI2","B_KO1","C_KO2","9_DEL1", "A_DEL2")

for (x in ser.order){
  print.ASP1<-profile.plotter(file.list1[x],f.mode="mean",f.name=ser.name[x],f.printmatrix=print.ASP1, f.printmatrix2=print.AS.WTP1, f.xlim=c(0,0), 
                             f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-500,0,750),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.5, 
                             f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation="")
}

for (x in ser.order){
  print.ASP2<-profile.plotter(file.list2[x],f.mode="mean",f.name=ser.name[x],f.printmatrix=print.ASP2, f.printmatrix2=print.AS.WTP2, f.xlim=c(0,0), 
                             f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-750,0,500),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=0.5, 
                             f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation="")
}

print.ASP1$KO<-(print.ASP1[,2]+print.ASP1[,3])/2
print.ASP1$RE<-(print.ASP1[,4]+print.ASP1[,5])/2
print.ASP1$DI<-(print.ASP1[,6]+print.ASP1[,7])/2
print.ASP1$DEL<-(print.ASP1[,8]+print.ASP1[,9])/2
names(print.ASP1)[10]<-"4_KO"
names(print.ASP1)[11]<-"1_RE"
names(print.ASP1)[12]<-"2_DI"
names(print.ASP1)[13]<-"3_DEL"

print.ASP2$KO<-(print.ASP2[,2]+print.ASP2[,3])/2
print.ASP2$RE<-(print.ASP2[,4]+print.ASP2[,5])/2
print.ASP2$DI<-(print.ASP2[,6]+print.ASP2[,7])/2
print.ASP2$DEL<-(print.ASP2[,8]+print.ASP2[,9])/2
names(print.ASP2)[10]<-"4_KO"
names(print.ASP2)[11]<-"1_RE"
names(print.ASP2)[12]<-"2_DI"
names(print.ASP2)[13]<-"3_DEL"

print.AS=print.ASP1[1]
for(x in 2:length(print.ASP1)){
  print.AS<-cbind(print.AS,print.ASP1[x],print.ASP2[x])
}

print.AS.WT<-cbind(print.AS.WTP1[,1:2], print.AS.WTP2[,2])  
names(print.AS.WT)[2]<-"WTP1"
names(print.AS.WT)[3]<-"WTP2"

aaa<-profile.plotter("",f.mode="mean",f.name="",f.printmatrix=print.AS, f.printmatrix2=print.AS.WT, f.xlim=c(0,0), 
                     f.ylim=c(0,0),f.xtrim=c(0,0),f.panels=c("A_5_profile","B_3_profile"),f.panel.middle=c(0,250), f.runmean.segment=FALSE,f.runmean=c(30,30,30),f.oldfilename=file.list1[1],f.log2=TRUE,
                     f.replace0to=0.5, f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=5, f.xticks=c(-500,0,250,750),f.middle=0,f.filter="EF2.mRNA.red1KO_affected_PROMTs.filter",f.mod.annotation=""
                     ,col=c("black","orange","blue","red","darkgrey","darkgrey","yellow","yellow","lightblue","lightblue","pink","pink"), lty=c(1,1,1,1,3,3,3,3,3,3,3,3),lwd=c(6,6,6,6,2,2,2,2,2,2,2,2))




################################################### ChIPseq heterochr islands composite plot ##################################################################
################################################### ChIPseq heterochr islands composite plot ##################################################################
################################################### ChIPseq heterochr islands composite plot ##################################################################
################################################### ChIPseq heterochr islands composite plot ##################################################################
################################################### ChIPseq heterochr islands composite plot ##################################################################
################################################### ChIPseq heterochr islands composite plot ##################################################################
################################################### ChIPseq heterochr islands composite plot ##################################################################
################################################### ChIPseq heterochr islands composite plot ##################################################################
################################################### ChIPseq heterochr islands composite plot ##################################################################
################################################### ChIPseq heterochr islands composite plot ##################################################################

setwd("~/ChIP-seq/bigwigs/covnorm_dusty_lt_4")

# GeneProfiler
file.list1<-dir(pattern="*.bw$")
file.list3<-substr(file.list1, start=1, stop=nchar(file.list1)-3)

f.annotation="~/S.pombe.EF2/EF2.heterochromaticIslands.Rdat"

lapply(1:length(file.list1), function(x) f.gene.profiler(file.list1[[x]],f.filename.rev="", f.savename=file.list3[[x]],f.annotation,f.sense=TRUE,f.binnumber=100,
                                                         f.binlength=20,f.upstream=500, f.downstream=500,f.cores=8))

# Profile plotter - Average genes

file.list1<-dir(pattern="*profile.Rdat$")
print.S.WT=integer()
WT.order=3

print.S.WT<-profile.plotter(file.list1[WT.order],f.mode="mean",f.name="WT",f.printmatrix=print.S.WT, f.printmatrix2=integer(), f.xlim=c(0,0), 
                            f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-500,0,2000,2500),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=1,f.convertbackto10=TRUE,
                            f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="~/Desktop/R/annotation/S.pombe.EF2/filters/EF2.heterochromaticIslands.red1depenndent.filter",f.mod.annotation="")

print.S=integer()
ser.order=c(3,4,5,1)

ser.name<-c("6_DEL","2_WT2","1_WT1","4_KO","5_DI","3_WT3")

for (x in ser.order){
  print.S<-profile.plotter(file.list1[x],f.mode="mean",f.name=ser.name[x],f.printmatrix=print.S, f.printmatrix2=print.S.WT, f.xlim=c(0,0), 
                           f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-500,0,2000,2500),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(5,5,5),f.log2=TRUE,f.replace0to=1,f.convertbackto10=TRUE, 
                           f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="~/Desktop/R/annotation/S.pombe.EF2/filters/EF2.heterochromaticIslands.red1depenndent.filter",f.mod.annotation="")
}


aaa<-profile.plotter("",f.mode="mean",f.name="",f.printmatrix=print.S, f.printmatrix2=print.S.WT, f.xlim=c(0,0), 
                     f.ylim=c(0,0),f.xtrim=c(0,0),f.panels=character(0), f.runmean.segment=FALSE,f.runmean=c(15,5,15),f.oldfilename=file.list1[1],f.log2=TRUE,f.replace0to=1,f.convertbackto10=TRUE,
                     f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=4, f.xticks=c(-500,0,2000,2500),f.middle=2000,f.filter="~/Desktop/R/annotation/S.pombe.EF2/filters/EF2.heterochromaticIslands.red1depenndent.filter",f.mod.annotation=""
                     ,col=c("black","red","orange","blue"),lwd=c(6,6,6,6))


