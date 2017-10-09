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

#================================
# take overlapping CpGs from BeadChip array annotations

#dim(getAnnotation(rg.epic.sub))
#epicanno <- as.data.frame(getAnnotation(rg.epic.sub))
#epic450.anno <- epicanno[!rownames(epicanno) %in% rownames(anno.hm450),]
#r1 <- anno.hm450[,c("chr","pos","strand","Name","Islands_Name","Relation_to_Island","UCSC_RefGene_Name","UCSC_RefGene_Group")]
#r2 <- epic450.anno[,c("chr","pos","strand","Name","Islands_Name","Relation_to_Island","UCSC_RefGene_Name","UCSC_RefGene_Group")]

#epic450.anno <- rbind(r1,r2)
#head(epic450.anno)

#save(epic450.anno,file="CpG_BeadChip_anno_EPIC450.rda")


#============================================================
# take intersect of granges objects, CpG anno and hsacoor

#


# new hsa miRNA df, every row is a unique range
a <- data.frame(mirna=c(hsacoor[1:529,1],hsacoor[530,1],hsacoor[530,1],hsacoor[531:nrow(hsacoor),1]),
                chr=c(hsacoor[1:529,2],hsacoor[530,2],hsacoor[530,2],hsacoor[531:nrow(hsacoor),2]),
                start=as.numeric(unlist(strsplit(paste(hsacoor[,3],collapse=";"),";"))),
                end=as.numeric(unlist(strsplit(paste(hsacoor[,4],collapse=";"),";"))),
                strand=c(hsacoor[1:529,5],hsacoor[530,5],hsacoor[530,5],hsacoor[531:nrow(hsacoor),5]),
                stringsAsFactors=FALSE);

# promoter df
ap <- a
pend <- ifelse(a$strand=="+",a$start-1000,a$start+1000) # designate promoter region ranges w.in 1KB relative to start, using strandedness
for(i in 1:nrow(ap)){
  pcoori <- c(ap$start[i],pend[i])
  ap$start[i] <- min(pcoori); ap$end[i] <- max(pcoori)
  message(round(100*(i/nrow(ap)),3),"%")
}; 

# GRanges object for hsa miRNA hairpin seq's
#grhsa <- makeGRangesFromDataFrame(a,keep.extra.columns = TRUE)
grhsap <- makeGRangesFromDataFrame(ap,keep.extra.columns = TRUE)
# GRanges object for CpG annotations from EPIC and HM450 arrays
cpg <- epic450.anno; cpg$start <- cpg$pos; cpg$end <- ifelse(cpg$strand=="+",cpg$pos+1,cpg$pos-1)
grcpg <- makeGRangesFromDataFrame(cpg,keep.extra.columns = TRUE)

#ai <- intersect(grhsa,grcpg,ignore.strand=TRUE)
#cpg.mirna <- grcpg[start(grcpg) %in% start(ai)]

api <- findOverlaps(grcpg,grhsap,ignore.strand=TRUE)
cpg.mirnaprom <- cpg[cpg$Name %in% mcols(grcpg[unique(queryHits(api))])$Name,]
mpcpg <- ap[ap$mirna %in% mcols(grhsap[unique(subjectHits(api))])$mirna,]
colnames(mpcpg) <- c("mirna","chr","start.prom","end.prom","strand")
xa <- a[a$mirna %in% mpcpg$mirna,]; identical(xa$mirna,mpcpg$mirna)
mpcpg$start.mirna <- xa$start; mpcpg$end.mirna <- xa$end

xcpg <- cpg.mirnaprom; xmirna <- mpcpg
xmirna$cpg <- NA
for(i in 1:nrow(xmirna)){
  which <- queryHits(findOverlaps(grcpg,
                                  grhsap[mcols(grhsap)$mirna==xmirna$mirna[i]],
                                  ignore.strand=TRUE))
  
  xmirna$cpg[i] <- paste(mcols(grcpg[which])$Name,
                         collapse=";")
  message(round(100*(i/nrow(xmirna)),3),"%")
}

xmirna$miR <- NA
for(i in 1:nrow(xmirna)){
  xmirna$miR[i] <- hsacoor[hsacoor$mirna==xmirna$mirna[i],]$miR
  message(round(100*(i/nrow(xmirna)),3),"%")
}

save(xmirna,file="mirna-anno_CpG-miR.rda")
