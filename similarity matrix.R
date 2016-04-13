setwd("~/Population Structure/Pasteuria")

library("proxy", lib.loc="~/R/win-library/3.2")
library("nomclust", lib.loc="~/R/win-library/3.2")
library("dplyr", lib.loc="~/R/win-library/3.2")

dat<- read.table("pastmatrix4.11.txt", header=TRUE, sep="\t", dec=".", strip.white=TRUE) # loads the data
info<- read.table("info.txt", header=TRUE, sep="\t", dec=".", strip.white=TRUE) # loads the data



SimMatch<-function(a,b){
  B<-setdiff(a,b)
  C<-setdiff(b,a)
  lB<-length(B)
  lC<-length(C)
  almostD<-as.vector(rbind(a,b))
  aD<-unique(almostD)
  D<-length(aD)
  Simplematch<-(lB+lC)/D
  return(Simplematch)}

#merge data
datinfo<-full_join(dat,info, by="sample")

type<-unique(datinfo$type)
#For loop to look at each group
OUT<-NULL   
for (i in type) {
  eachtype<-datinfo[datinfo$type==i,] #get each group by itself
  eachtyp1<-eachtype[,1:13]
  only1<-toString(unique(eachtyp1$p1))
  only2<-toString(unique(eachtyp1$p2))
  only3<-toString(unique(eachtyp1$p3))
  only4<-toString(unique(eachtyp1$p4))
  only7<-toString(unique(eachtyp1$p7))
  only9<-toString(unique(eachtyp1$p9))
  only11<-toString(unique(eachtyp1$p11))
  only12<-toString(unique(eachtyp1$p12))
  only16<-toString(unique(eachtyp1$p16))
  only17<-toString(unique(eachtyp1$p17))
  only18<-toString(unique(eachtyp1$p18))
  only19<-toString(unique(eachtyp1$p19))

  p1<-c(i, "p1", only1)
  p2<-c(i, "p2", only2)
  p3<-c(i,"p3",only3)
  p4<-c(i,"p4",only4)
  p7<-c(i,"p7", only7)
  p9<-c(i,"p9", only9)
  p11<-c(i,"p11",only11)
  p12<-c(i, "p12", only12)
  p16<-c(i,"p16", only16)
  p17<-c(i,"p17", only17)
  p18<-c(i,"p18",only18)
  p19<-c(i,"p19", only19)
  OUT<-rbind(OUT, p1,p2,p3,p4,p7,p9,p11,p12,p16,p17,p18,p19)
}
colnames(OUT)<-c("type","primer","alleles")
write.table(OUT, "bygroups.txt", col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE)

dat1<- read.table("bygroups.txt", header=TRUE, sep="\t", dec=".", strip.white=TRUE) # loads the data
prim<-unique(dat1$primer)
type<-unique(dat1$type)
listy<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)

OUT<-NULL
for (i in prim){
  eachprimer<-dat1[dat1$primer==i,]
for (j in type){
  eachtype<-eachprimer[eachprimer$type==j,] #get each group by itself
    L<-as.character(eachtype$alleles)
    a<-unlist(strsplit(L, split=", "))
    a<-a[a != "NA"]
for (k in listy){
      othertype<-eachprimer[eachprimer$type==k,]
      M<-as.character(othertype$alleles)
      b<-unlist(strsplit(M, split=", "))
      b<-b[b != "NA"]
      val<-SimMatch(a,b)
      output<-c(i, j, k, val)
      OUT<-rbind(OUT, output)
      }
  }
}

colnames(OUT)<-c("primer","type","type2","SimMatch")
write.table(OUT, "simmatch.txt", col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE)

sims<- read.table("simmatch.txt", header=TRUE, sep="\t", dec=".", strip.white=TRUE) # loads the data
shortinfo<-info[,-1]
shortinfo<-unique(shortinfo)
simsid<-inner_join(sims, shortinfo, by=c("type2"="type"))
newshortinfo<-shortinfo
newshortinfo$species1<-newshortinfo$species
simsid1<-inner_join(simsid, newshortinfo, by="type")
notconfusing<-select(simsid1, primer, SimMatch, species.x, lake.x, date.x, lake.y, date.y, species.y)
notconfusing<-rename(notconfusing, c("species.x"="compare.species", "date.x"="compare.date", "lake.x"="compare.lake", "date.y"="date", "species.y"="species", "lake.y"="lake"))
final<-notconfusing[,c(1,8,6,7,3,4,5,2)]

#alldentifera  
Dent<-simsid[simsid$species=="dentifera",]
cDent<-Dent[Dent$type2==c(1,2,5,6,7,8,13,14),]
dentavg<-mean(cDent$SimMatch)

#dent to retro
cDRetro<-Dent[Dent$type2==c(3, 11, 12,15,16),]
cDRavg<-mean(cDRetro$SimMatch)

#dent to parvula
cDParv<-Dent[Dent$type2==c(10, 17),]
cDPavg<-mean(cDParv$SimMatch)

#early to late same lake
Tm<-cDent[cDent$]

exact<-function(a,b){
  mean(a==b, na.rm=TRUE)
}

dat<- read.table("pastmatrix3.31fiddles.txt", header=TRUE, sep="\t", dec=".", strip.white=TRUE) # loads the data
info<- read.table("info.txt", header=TRUE, sep="\t", dec=".", strip.white=TRUE) # loads the data

datinfo<-full_join(dat,info, by="sample")

type<-unique(datinfo$type)

#For loop for looking at similarity within a group
OUT<-NULL
for (i in type){
  eachtype<-datinfo[datinfo$type==i,]
  p1lines<-NULL
  p2lines<-NULL
  p3lines<-NULL
  p4lines<-NULL
  p7lines<-NULL
  p9lines<-NULL
  p11lines<-NULL
  p12lines<-NULL
  p16lines<-NULL
  p17lines<-NULL
  p18lines<-NULL
  p19lines<-NULL
    for (j in 1:length(na.omit(eachtype$p1))){
      for (m in 1:length(na.omit(eachtype$p1))){
      thisgroup<-select(eachtype, sample, p1)  
      thisgroup<-na.omit(thisgroup)
      ori<-thisgroup$sample[j]
      comp<-thisgroup$sample[m] 
      a<-thisgroup$p1[j]
      a[is.na(a)] <- 0 
      b<-thisgroup$p1[m]
      b[is.na(b)] <- 0 
      if (a==b){
        p1.match<-1 
        } else {
          p1.match<-0
        }
      p1line<-c(ori, comp, p1.match)
      p1lines<-rbind(p1lines, p1line)
      p1lines<-na.omit(p1lines)
      colnames(p1lines)<-c("ori","comp","match")
      }
      }
  for (j in 1:length(na.omit(eachtype$p2))){
    for (m in 1:length(na.omit(eachtype$p2))){
      thisgroup<-select(eachtype, sample, p2)  
      thisgroup<-na.omit(thisgroup)
      ori<-thisgroup$sample[j]
      comp<-thisgroup$sample[m] 
      a<-thisgroup$p2[j]
      a[is.na(a)] <- 0 
      b<-thisgroup$p2[m]
      b[is.na(b)] <- 0 
      if (a==b){
        p2.match<-1 
      } else {
        p2.match<-0
      }
      p2line<-c(ori, comp, p2.match)
      p2lines<-rbind(p2lines, p2line)
      p2lines<-na.omit(p2lines)
      colnames(p2lines)<-c("ori","comp","match")
    }
    }
  for (j in 1:length(na.omit(eachtype$p3))){ 
    for (m in 1:length(na.omit(eachtype$p3))){
      thisgroup<-select(eachtype, sample, p3)  
      thisgroup<-na.omit(thisgroup)
      ori<-thisgroup$sample[j]
      comp<-thisgroup$sample[m] 
      a<-thisgroup$p3[j]
      a[is.na(a)] <- 0 
      b<-thisgroup$p3[m]
      b[is.na(b)] <- 0 
      if (a==b){
        p3.match<-1 
        } else {
        p3.match<-0
      }
      p3line<-c(ori, comp, p3.match)
      p3lines<-rbind(p3lines, p3line)
      p3lines<-na.omit(p3lines)
      colnames(p3lines)<-c("ori","comp","match")
    }}
  for (j in 1:length(na.omit(eachtype$p4))){ #problem here for i=12 because none have a value
    for (m in 1:length(na.omit(eachtype$p4))){
      thisgroup<-select(eachtype, sample, p4)  
      thisgroup<-na.omit(thisgroup)
      ori<-thisgroup$sample[j]
      comp<-thisgroup$sample[m] 
      a<-thisgroup$p4[j]
      a[is.na(a)] <- 0 
      b<-thisgroup$p4[m]
      b[is.na(b)] <- 0 
      if (a==b){
        p4.match<-1 
      } else{
        p4.match<-0
      }
      p4line<-c(ori, comp, p4.match)
      p4lines<-rbind(p4lines, p4line)
      p4lines<-na.omit(p4lines)
      colnames(p4lines)<-c("ori","comp","match")
    }}
  for (j in 1:length(na.omit(eachtype$p7))){
    for (m in 1:length(na.omit(eachtype$p7))){
      thisgroup<-select(eachtype, sample, p7)  
      thisgroup<-na.omit(thisgroup)
      ori<-thisgroup$sample[j]
      comp<-thisgroup$sample[m] 
      a<-thisgroup$p7[j]
      a[is.na(a)] <- 0 
      b<-thisgroup$p7[m]
      b[is.na(b)] <- 0 
      if (a==b){
        p7.match<-1 
        } else {
        p7.match<-0
      }
      p7line<-c(ori, comp, p7.match)
      p7lines<-rbind(p7lines, p7line)
      p7lines<-na.omit(p7lines)
      colnames(p7lines)<-c("ori","comp","match")
    }}
  for (j in 1:length(na.omit(eachtype$p9))){
    for (m in 1:length(na.omit(eachtype$p9))){
      thisgroup<-select(eachtype, sample, p9)  
      thisgroup<-na.omit(thisgroup)
      ori<-thisgroup$sample[j]
      comp<-thisgroup$sample[m] 
      a<-thisgroup$p9[j]
      a[is.na(a)] <- 0 
      b<-thisgroup$p9[m]
      b[is.na(b)] <- 0 
      if (a==b){
        p9.match<-1 }
      else {
        p9.match<-0
      }
      p9line<-c(ori, comp, p9.match)
      p9lines<-rbind(p9lines, p9line)
      p9lines<-na.omit(p9lines)
      colnames(p9lines)<-c("ori","comp","match")
    }}
  for (j in 1:length(na.omit(eachtype$p11))){
    for (m in 1:length(na.omit(eachtype$p11))){
      thisgroup<-select(eachtype, sample, p11)  
      thisgroup<-na.omit(thisgroup)
      ori<-thisgroup$sample[j]
      comp<-thisgroup$sample[m] 
      a<-thisgroup$p11[j]
      a[is.na(a)] <- 0 
      b<-thisgroup$p11[m]
      b[is.na(b)] <- 0 
      if (a==b){
        p11.match<-1 }
      else {
        p11.match<-0
      }
      p11line<-c(ori, comp, p11.match)
      p11lines<-rbind(p11lines, p11line)
      p11lines<-na.omit(p11lines)
      colnames(p11lines)<-c("ori","comp","match")
    }}
  for (j in 1:length(na.omit(eachtype$p12))){
    for (m in 1:length(na.omit(eachtype$p12))){
      thisgroup<-select(eachtype, sample, p12)  
      thisgroup<-na.omit(thisgroup)
      ori<-thisgroup$sample[j]
      comp<-thisgroup$sample[m] 
      a<-thisgroup$p12[j]
      a[is.na(a)] <- 0 
      b<-thisgroup$p12[m]
      b[is.na(b)] <- 0 
      if (a==b){
        p12.match<-1 }
      else {
        p12.match<-0
      }
      p12line<-c(ori, comp, p12.match)
      p12lines<-rbind(p12lines, p12line)
      p12lines<-na.omit(p12lines)
      colnames(p12lines)<-c("ori","comp","match")
    }}
  for (j in 1:length(na.omit(eachtype$p16))){
    for (m in 1:length(na.omit(eachtype$p16))){
      thisgroup<-select(eachtype, sample, p16)  
      thisgroup<-na.omit(thisgroup)
      ori<-thisgroup$sample[j]
      comp<-thisgroup$sample[m] 
      a<-thisgroup$p16[j]
      a[is.na(a)] <- 0 
      b<-thisgroup$p16[m]
      b[is.na(b)] <- 0 
      if (a==b){
        p16.match<-1 
        } else {
        p16.match<-0
      }
      p16line<-c(ori, comp, p16.match)
      p16lines<-rbind(p16lines, p16line)
      p16lines<-na.omit(p16lines)
      colnames(p16lines)<-c("ori","comp","match")
    }}
  for (j in 1:length(na.omit(eachtype$p17))){
    for (m in 1:length(na.omit(eachtype$p17))){
      thisgroup<-select(eachtype, sample, p17)  
      thisgroup<-na.omit(thisgroup)
      ori<-thisgroup$sample[j]
      comp<-thisgroup$sample[m] 
      a<-thisgroup$p17[j]
      a[is.na(a)] <- 0 
      b<-thisgroup$p17[m]
      b[is.na(b)] <- 0 
      if (a==b){
        p17.match<-1 }
      else {
        p17.match<-0
      }
      p17line<-c(ori, comp, p17.match)
      p17lines<-rbind(p17lines, p17line)
      p17lines<-na.omit(p17lines)
      colnames(p17lines)<-c("ori","comp","match")
    }}
  for (j in 1:length(na.omit(eachtype$p18))){
    for (m in 1:length(na.omit(eachtype$p18))){
      thisgroup<-select(eachtype, sample, p18)  
      thisgroup<-na.omit(thisgroup)
      ori<-thisgroup$sample[j]
      comp<-thisgroup$sample[m] 
      a<-thisgroup$p18[j]
      a[is.na(a)] <- 0 
      b<-thisgroup$p18[m]
      b[is.na(b)] <- 0 
      if (a==b){
        p18.match<-1 }
      else {
        p18.match<-0
      }
      p18line<-c(ori, comp, p18.match)
      p18lines<-rbind(p18lines, p18line)
      p18lines<-na.omit(p18lines)
      colnames(p18lines)<-c("ori","comp","match")
    }}
  for (j in 1:length(na.omit(eachtype$p19))){
    for (m in 1:length(na.omit(eachtype$p19))){
      thisgroup<-select(eachtype, sample, p19)  
      thisgroup<-na.omit(thisgroup)
      ori<-thisgroup$sample[j]
      comp<-thisgroup$sample[m] 
      a<-thisgroup$p19[j]
      a[is.na(a)] <- 0 
      b<-thisgroup$p19[m]
      b[is.na(b)] <- 0 
      if (a==b){
        p19.match<-1 }
      else {
        p19.match<-0
      }
      p19line<-c(ori, comp, p19.match)
      p19lines<-rbind(p19lines, p19line)
      p19lines<-na.omit(p19lines)
      colnames(p19lines)<-c("ori","comp","match")
    }}
  all.lines<-rbind(p1lines, p2lines, p3lines, p4lines, p7lines, p9lines, p11lines, p12lines, p16lines, p17lines, p18lines, p19lines) 
  colnames(all.lines)=c("ori","comp","match")
  all.lines <-as.data.frame(all.lines)
  all.linesrm<-all.lines[all.lines$ori != all.lines$comp, ]
  avg.rel <- all.linesrm %>% group_by(ori, comp) %>% summarise(relate=mean(match))
  totalavg<-mean(avg.rel$relate)
  output<-c(i, totalavg)
  OUT<-rbind(OUT, output)
}

colnames(OUT)<-c("type","relatedness")
write.table(OUT, "sim.within.txt", col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE)

simwithin<- read.table("sim.within.txt", header=TRUE, sep="\t", dec=".", strip.white=TRUE) # loads the data
shortinfo<-info[,-1]
shortinfo<-unique(shortinfo)
simwithin1<-inner_join(simwithin, shortinfo, by="type")
simwithin2<-select(simwithin1, relatedness, species, lake, date)
simwithin3<-simwithin2[,c(2,3,4,1)]



#the thing that CSCAR wrote for me
rownames(eachtyp1)<-eachtyp1[,1]
  eachtyp1<-eachtyp1[,-1]
  sim<-simil(eachtyp1,method=exact, by_rows=TRUE)
  simmat<-as.matrix(sim)
  round(simmat,2)
  avgrel<-mean(simmat, na.rm=TRUE)





