#Daniel Ence and Jeremy Brawner
#June 28, 2020

library(tidyverse)
library(rrBLUP)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

plate_layout_filename <- snakemake@params[["data_plate_layout"]]
plate_extractions_filename <- snakemake@params[["data_plate_extractions"]]

#Plate_layout <- read.table("./UFL_104301_Plate_Layout_Dec2016.txt",sep="\t",header=T)
#plate_extractions <- read.table("./DNA_plate_layout.txt",sep="\t",header=T)

Plate_layout <- read.table(plate_layout_filename,sep="\t",header=T)
plate_extractions <- read.table(plate_extractions_filename,sep="\t",header=T)

head(plate_extractions)
tail(plate_extractions)

#tmp_plate_extractions <- plate_extractions %>% filter(Gall != "BLANK") %>% filter(Gall != "RNAseq exp") %>% filter(Gall != "")
tmp_plate_extractions <- filter(plate_extractions, Gall != "BLANK")
tmp_plate_extractions <- filter(tmp_plate_extractions, Gall != "RNAseq exp")
tmp_plate_extractions <- filter(tmp_plate_extractions, Gall != "")
head(tmp_plate_extractions)
tail(tmp_plate_extractions)

# tmp_Plate_layout <- Plate_layout %>% select("RAPiD.Genomics.Sample.Code","Unique.Client.Code","Comments")
tmp_Plate_layout <- select(Plate_layout , "RAPiD.Genomics.Sample.Code","Unique.Client.Code","Comments")
head(tmp_Plate_layout)
head(tmp_plate_extractions)
tail(tmp_Plate_layout)
tail(tmp_plate_extractions)

# plate_and_ID <- inner_join(tmp_Plate_layout,tmp_plate_extractions,by=c("Unique.Client.Code"="plant.ID")) %>% select("RAPiD.Genomics.Sample.Code","Gall","plate","column","row")
plate_and_ID <- inner_join(tmp_Plate_layout,tmp_plate_extractions,by=c("Unique.Client.Code"="plant.ID"))

plate_and_ID <- select(plate_and_ID, "RAPiD.Genomics.Sample.Code","Gall","plate","column","row")
#plate_and_ID$Unique.Client.Code <- as.factor(plate_and_ID$Unique.Client.Code)
plate_and_ID<-droplevels(plate_and_ID)
str(plate_and_ID)
head(plate_and_ID)
tail(plate_and_ID)
table(plate_and_ID$Gall)
unique(plate_and_ID$Gall)
unique(plate_and_ID$RAPiD.Genomics.Sample.Code)
plate_and_ID$galled <- 1
sum(plate_and_ID$galled)
plate_and_ID$galled[plate_and_ID$Gall=="no gall"] <- 0
sum(plate_and_ID$galled)
#galled is 0 1, binary data
head(plate_and_ID)
tail(plate_and_ID)
dim(plate_and_ID)
plate_and_ID


dat.012 <- read.table(snakemake@input[["data_012"]],header=F,sep="\t")

#dat.012 <- read.table("./out.012",header=F,sep="\t")
dat.012[1:20,1:15]; dim(dat.012)
# sum(dat.012[1,])
length(unique(dat.012[,1]))

#
dat.012.indv <- read.table(snakemake@input[["indv"]],header=F)
dat.012.indv$fixed <- paste("UFL_104301_",dat.012.indv$V1,sep="")
dat.012.indv[]; dim(dat.012.indv)
head(dat.012.indv)
length(unique(dat.012.indv[,1]))



tmp_dat.012 <- t(dat.012)
colnames(tmp_dat.012) <- dat.012.indv$V1
tmp_dat.012 <- as.data.frame(tmp_dat.012)[-1,]
tmp_dat.012[1:20,1:5]; dim(tmp_dat.012)
dat.012 <- as.data.frame(dat.012)[,-1]
dat.012[1:6,1:15]; dim(dat.012)
#
sum(is.na(dat.012) ) # all there
#rm(dat.012)#

#How do the names from markers in dat.012 connect to the gall data in plate_and_ID
gen<-tmp_dat.012
#rm(tmp_dat.012)#

length(colnames(gen)) #names of samples with genotypes colnames(gen)



phen <- as.data.frame(matrix(0,length(colnames(gen)),1))
head(plate_and_ID)
plate_and_ID$id <- substr(plate_and_ID$RAPiD.Genomics.Sample.Code,12,99)
phen$id  <-  colnames(gen)

idphen <- merge(phen, plate_and_ID[,c("id","galled")], by.y="id")
#idphen

mean(idphen$galled) #50.17


# mom and dad?
str(idphen)
class(idphen)
str(gen)
class(gen)

gent<-t(gen)
str(gent)
class(gent)



#library(asreml)
#source('pin.R')
#library(rrBLUP)
#library(GGally)
#asreml.license.offline(28)
#asreml.license.online()

# Marker data
# gen<-read.delim("out.012.pos", stringsAsFactors = FALSE, quote = "", sep="\t", as.is=T)  # No missing values, same order of individuals as data
# tmp_dat.012[1:5,1:15]; dim(tmp_dat.012)
#gen<-read.delim("out.012", stringsAsFactors = FALSE, quote = "", sep="\t", as.is=T)  # No missing values, same order of individuals as data

#gen[1:5,1:5]; dim(gen)
nmarkers <- length(gent[1,])   #495381 cols
nseedlings <- length(gent[,1]) #291    rows
gent<-gent-1; min(gent); max(gent) # recode -1 0 1
gen<-gen-1
dim(gent)
gent[1:5,1:5]
gen[1:5,1:5]
length(rownames(gent)) == length(unique(rownames(gent)))
length(colnames(gent)) == length(unique(colnames(gent)))

#rename things, no...
#rownames(gen)<-paste("PT",gen[,1],sep="") #ID is in 1
#gen<-gen[,-1] #remove ID column
#colnames(gent)[1:100]<-paste("M",1:nmarkers,sep="")
#gen<-as.matrix(gen)
class(gent)

#could use gapit
#source("http://zzlab.net/GAPIT/gapit_functions.txt")

##rrBLUP
amat<-rrBLUP::A.mat(gent)
amat[1:5,1:5] #top left
gened<-length(amat[1,])
amat[(gened-5):gened,(gened-5):gened] #bottom right
amat[1:5,(gened-5):gened] #top right





#Duplicates in phenotype... does it  matter which one? #this should move up to data cleansing...
dim(idphen)
class(idphen)
setdiff(idphen$id, rownames(amat))
setdiff(rownames(amat), idphen$id)
duplicated(rownames(amat))
dups<-droplevels(idphen[duplicated(idphen$id),])
str(dups)
idphen[idphen$id=="P01_WF12",] #duplicates have same phen, same gen only has one phen so -  it doesnt mater
str(idphen)
idphen[idphen$id %in% dups[,1],] #same, just make unique
head(plate_and_ID)
plate_and_ID[plate_and_ID$id %in% dups[,1],]
idphen<-unique(idphen) #drop repeat
idphen[idphen$id %in% dups[,1],] #same, just make unique

dim(idphen)
dim(gent)
rownames(idphen)
head(idphen)
idphen<-idphen[,c(1,3)]
head(idphen)


gen<-as.data.frame(gen)
#str(gent)
starts<-seq(from =1, to =(length(rownames(gen))*2),by=2)
stops<-seq(from =2, to =(length(rownames(gen))*2),by=2)
ss<-as.data.frame(cbind(starts,stops ))
mname<-paste("M",1:nmarkers,sep="")
idss<-as.data.frame(cbind(mname,ss))
colnames(idss)<-c("id","starts","stops")
tail(idss)
#idss[,colnames(idss)] <- as.character(idss[,colnames(idss)])
idss[,2] <- as.numeric(idss[,2])
idss[,3] <- as.numeric(idss[,3])
str(idss)
dim(idss)
dim(gen)
#rownames(gent)<-NULL
gen<-cbind(idss,gen)

gen[1:5,1:5]; dim(gen)
#gent[1:5,495377:495384]
#rownames(gent)<-NULL

dim(idphen)
dim(gen)
dim(amat)

#set the positions to keep track of markers
positions <- read.table(snakemake@input[["pos_012"]],header=F,sep="\t")
gen$id <- positions$V1
gen$starts <- positions$V2
gen$stops <- positions$V2 + 1

gwas1 <- rrBLUP::GWAS(idphen, gen, fixed=NULL, K=amat, n.PC=0,
                      min.MAF=0.05, n.core=1, P3D=TRUE, plot=TRUE)
head(gwas1)
tail(gwas1)


hist(gwas1$galled[gwas1$galled>.1])

damat<- amat
diag(damat)<-0
image(damat)
mean(diag(amat))

#gen<-gen-1 #code genotypes as -1,0,1
#gen[280:290,1:10]
#sum(is.na(gen[,1:nmarkers])) #no missing markers
length(rownames(gent))
dim(gent)
# phenotype data
dim(idphen)
