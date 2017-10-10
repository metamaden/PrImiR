#===================================================
# mapping miRNA promoter CPGs on methylation arrays
#
# online resource for genome mappings: ftp://mirbase.org/pub/mirbase/CURRENT/genomes/

# use bioconductor package, approx. latest available version for Ch37/hg19 => switched to Ch38 after..
#source("https://bioconductor.org/biocLite.R")
#biocLite("mirbase.db")

library(mirbase.db)

mirloc.start <- mirbaseCHRLOC
mirloc.end <- mirbaseCHRLOCEND

mkstart <- mappedkeys(mirloc.start)
xxstart <- as.list(mirloc.start[mkstart])
xxhsa.start <- xxstart[grep("hsa",names(xxstart))]

mkend <- mappedkeys(mirloc.end)
xxend <- as.list(mirloc.end[mkend])
xxhsa.end <- xxend[grep("hsa",names(xxend))]

length(intersect(names(xxhsa.start),names(xxhsa.end)))
identical(names(xxhsa.start),names(xxhsa.end))

hsacoor <- as.data.frame(matrix(nrow=length(xxhsa.start),ncol=5))
colnames(hsacoor) <- c("mirna","chr","start","end","strand")
for(i in 1:nrow(hsacoor)){
  hsacoor[i,1] <- names(xxhsa.start)[i]
  hsacoor[i,2] <- unique(names(xxhsa.start[[i]]))
  
  if(length(xxhsa.start[[i]])==1){
    hsacoor[i,3] <- gsub("-","",as.numeric(xxhsa.start[[i]]))
    hsacoor[i,4] <- gsub("-","",as.numeric(xxhsa.end[[i]]))
    hsacoor[i,5] <- ifelse(grepl("-",as.character(xxhsa.start[[i]])),"-","+")
    
  }
  if(length(xxhsa.start[[i]])>1){
    hsacoor[i,3] <- paste(gsub("-","",as.numeric(xxhsa.start[[i]])),collapse=";")
    hsacoor[i,4] <- paste(gsub("-","",as.numeric(xxhsa.end[[i]])),collapse=";")
    hsacoor[i,5] <- unique(ifelse(unique(grepl("-",as.character(xxhsa.start[[i]]))),"-","+"))
  }
  message(paste0(i," or ",round(100*(i/nrow(hsacoor)),3),"%"))
 
}

#=======================================#
# map hairpins to mature miR sequences

x <- mirbaseMATURE

hsacoor$miR <- NA
for(i in 1:nrow(hsacoor)){
  hsacoor[i,6] <- paste(matureName(get(hsacoor[i,1],x)),collapse=";")
  message(round(100*(i/nrow(hsacoor))),"%")
}

save(hsacoor,file="mirna-hairpin_chrcoor-df.rda")

#============================================================
# take intersect of granges objects, CpG anno and hsacoor

# new hsa miRNA df, every row is a unique range
a <- data.frame(mirna=c(hsacoor[1:529,1],hsacoor[530,1],hsacoor[530,1],hsacoor[531:nrow(hsacoor),1]),
                chr=c(hsacoor[1:529,2],hsacoor[530,2],hsacoor[530,2],hsacoor[531:nrow(hsacoor),2]),
                start=as.numeric(unlist(strsplit(paste(hsacoor[,3],collapse=";"),";"))),
                end=as.numeric(unlist(strsplit(paste(hsacoor[,4],collapse=";"),";"))),
                strand=c(hsacoor[1:529,5],hsacoor[530,5],hsacoor[530,5],hsacoor[531:nrow(hsacoor),5]),
                stringsAsFactors=FALSE);

a$start.promoter <- NA; a$end.promoter <- NA

pend <- ifelse(a$strand=="+",a$start-1500,a$start+1500) # designate promoter region ranges w.in 1.5KB relative to start, using strandedness

for(i in 1:nrow(a)){
  pcoori <- c(a$start[i],pend[i])
  a$start.promoter[i] <- min(pcoori); a$end.promoter[i] <- max(pcoori)
  message(round(100*(i/nrow(a)),3),"%")
}; 

# GRanges object for hsa miRNA hairpin seq's
#grhsa <- makeGRangesFromDataFrame(a,keep.extra.columns = TRUE)
x <- a; colnames(x) <- c(colnames(x)[1:2],"start.mirna","end.mirna",colnames(x)[5],"start","end")
grhsap <- makeGRangesFromDataFrame(x,keep.extra.columns = TRUE)

# GRanges object for CpG annotations from EPIC and HM450 arrays
cpg <- epic450.anno; cpg$start <- cpg$pos; cpg$end <- ifelse(cpg$strand=="+",cpg$pos+1,cpg$pos-1)
grcpg <- makeGRangesFromDataFrame(cpg,keep.extra.columns = TRUE)

# GRanges object for only 450k CpGs
cpg450 <- r1; cpg450$start <- cpg450$pos; cpg450$end <- ifelse(cpg450$strand=="+",cpg450$pos+1,cpg450$pos-1)
grcpg450 <- makeGRangesFromDataFrame(cpg450,keep.extra.columns = TRUE)

# GRanges object for only EPIC CpGs
cpgepic <- r3; cpgepic$start <- cpgepic$pos; cpgepic$end <- ifelse(cpgepic$strand=="+",cpgepic$pos+1,cpgepic$pos-1)
grcpgepic <- makeGRangesFromDataFrame(cpgepic,keep.extra.columns = TRUE)

mcpg <- a

xmirna <- mpcpg
xmirna$cpg <- NA; xmirna$cpg450 <- NA; xmirna$cpgepic <- NA
for(i in 1:nrow(xmirna)){
  # All BeadChip overlapping CpGs
  which <- queryHits(findOverlaps(grcpg,
                                  grhsap[mcols(grhsap)$mirna==xmirna$mirna[i]],
                                  ignore.strand=TRUE))
  
  xmirna$cpg[i] <- paste(mcols(grcpg[which])$Name,
                         collapse=";")
  
  # HM450 overlapping CpGs
  which450 <- queryHits(findOverlaps(grcpg450,
                                  grhsap[mcols(grhsap)$mirna==xmirna$mirna[i]],
                                  ignore.strand=TRUE))
  xmirna$cpg450[i] <- paste(mcols(grcpg450[which450])$Name,
                         collapse=";")
  
  # EPIC overlapping CpGS
  whichepic <- queryHits(findOverlaps(grcpgepic,
                                     grhsap[mcols(grhsap)$mirna==xmirna$mirna[i]],
                                     ignore.strand=TRUE))
  xmirna$cpgepic[i] <- paste(mcols(grcpgepic[whichepic])$Name,
                            collapse=";")
  
  message(round(100*(i/nrow(xmirna)),3),"%")
}

xmirna$miR <- NA
for(i in 1:nrow(xmirna)){
  xmirna$miR[i] <- hsacoor[hsacoor$mirna==xmirna$mirna[i],]$miR
  message(round(100*(i/nrow(xmirna)),3),"%")
}

save(xmirna,file="mirna-anno_CpG-450-epic-miR.rda")

###
