# PrImiR
Pairwise predicted Interactions for miRNAs 

### Description
Database of pairwise consensuses for miR-mRNA target interactions

This package is downloadable using the following in R:
```require(devtools); install_github("metamaden/PrImiR")```

### Background
The science of miR target prediction is actively evolving. There is no universal consensus for all computationally predicted miR-mRNA target interactions, and in addition most computationally predicted interactions have not been validated in the wet lab. To mitigate the limitations of conditions for individual computational algorithms, this package provides a consensus database of predicted miRs and miRNA-mRNA target interactions using the latest available versions of five established and widely used miR databases and target prediction algorithms:

1. Target Scan (v. 7.1)
2. Miranda (v)
3. miRDB (v)
4. miRTarBase (v)
5. miRmap (v)

## Supplemental Information
###1. Database Consensus Visualization: Simple Heatmap of random subset (N=2k miR-mRNA interactions) of PrImiR pairwise consensus miR-mRNA interactions database, showing presence of miR-mRNA interactions (rows) across databases (columns). blue=0/not present; red=1/present

![PrImiR pairwise consensus heatmap summary](https://github.com/metamaden/methyPre/blob/master/hm-consensus_interactions_2k.jpg "PrImiR Interactions")

The above heatmap was generated with the following code (cites objects in SI script 2, see below): 
```library(ComplexHeatmap)
xhm <- as.matrix(pci.df[sample(nrow(pci.df),2000),3:7])
class(xhm) <- "numeric"

Heatmap(xhm,show_column_names = TRUE,
        show_row_names = FALSE,
        cluster_columns = FALSE,
        cluster_rows=FALSE,
        name="Interaction\nin DB",
        row_title = "Predicted miR-mRNA\nInteraction (N=2000)",
        column_title="Databases and Algorithms")```

###2. Details for assembling pairwise consensus interaction databases.
The following script was used to assemble the pairwise consensus matrix of predicted miR-mRNA target interactions.
```#=========================================
# assemble miR target consensus database
# author: Sean Maden
library(org.Hs.eg.db); library(data.table)

# TargetScan 7.1 targets
# http://www.targetscan.org/cgi-bin/targetscan/data_download.vert71.cgi
{
  ts <- read.table("Conserved_Site_Context_Scores_targetscan71.txt",stringsAsFactors = FALSE,skip=1)
  ts <- ts[grep("hsa",ts[,5]),]
}

# miranda targets
# http://www.microrna.org/microrna/getDownloads.do
{
  mda <- fread("hg19_predictions_S_C_aug2010_miranda.txt",
               select=c("#mirbase_acc","mirna_name","gene_id","gene_symbol"),
               stringsAsFactors = FALSE,data.table=FALSE); class(mda)
  colnames(mda)[1] <- c("mirbase_acc")
}

# miRmap
# http://mirmap.ezlab.org/
{
  x <- fread("mirmap201301e_homsap_targets.csv",select=c(1,2,8,9),data.table=FALSE)
  grbg <- paste0(x$mirna_id,x$gene_name); length(grbg[duplicated(grbg)])
  mirmap.unique <- x2 <- x[!duplicated(grbg),]
}

# miRTarBase targets, high quality
# http://mirtarbase.mbc.nctu.edu.tw/php/download.php
{
  mtb <- read.csv("miRTarBase_SE_WR.csv",stringsAsFactors = FALSE)
  mtb <- mtb[grep("hsa",mtb$miRNA),]; mtb <- mtb[mtb$Species..miRNA.=="Homo sapiens",] 
  
}

# miRDB targets
# http://www.mirdb.org/download.html
{
  mdb <- read.table("miRDB_v5.0_prediction_result.txt",stringsAsFactors = FALSE)
  mdb <- mdb[grep("hsa",mdb[,1]),]
  colnames(mdb) <- c("mirid","target.refseqid","score")
  
  symbol.map <- as.list(org.Hs.egSYMBOL)
  acc.map <- as.list(org.Hs.egREFSEQ2EG) 
  x <- unique(mdb$target.refseqid)
  
  xfilt <- x[x %in% names(acc.map[as.character(acc.map) %in% names(symbol.map)])]
  
  df.map <- data.frame(db.acc=xfilt,ent=rep(NA,length(xfilt)),symbol=rep(NA,length(xfilt)))
  for(i in 1:nrow(df.map)){
    df.map[i,2] <- as.character(acc.map[names(acc.map)==df.map[i,1]])
    df.map[i,3] <- as.character(symbol.map[names(symbol.map)==as.character(acc.map[names(acc.map)==df.map[i,1]])])
    message(round(100*i/nrow(df.map),3),"%")
  }
  
  mdb.filt <- mdb[mdb$target.refseqid %in% df.map[,1],]; dim(mdb.filt)
  
  x2 <- mdb.filt$target.refseqid
  repx2 <- c()
  for(i in 1:nrow(df.map)){
    repx2[i]<-length(x2[x2==df.map[i,1]]);
    message(round(100*i/length(x),3),"%")
  }
  
  symrep <- c()
  for(i in 1:nrow(df.map)){
    symrep <- c(symrep,rep(df.map[i,3],repx2[i]))
    message(round(100*i/nrow(df.map),3),"%")
  }
  
  mdb.filt$target.symbol <- symrep; save(mdb.filt,file="mirdb_geneid-symbol.rda")
  
}

#==================================
# 
#   common miRs (in at least 2x db)
#
#==================================
head(mdb.filt); head(mda); head(ts); head(mirmap.unique); head(mtb)

# ts intersect
keep.mir.mda.ts <- intersect(unique(mda$mirna_name),unique(ts[,5]))
keep.mir.mdb.ts <- intersect(unique(mdb.filt$mirid),unique(ts[,5]))
keep.mir.mtb.ts <- intersect(unique(mtb$miRNA),unique(ts[,5]))
keep.mir.mmap.ts <- intersect(unique(mirmap.unique$mature_name),unique(ts[,5]))

# miranda intersect
keep.mir.mdb.mda <- intersect(unique(mdb.filt$mirid),unique(mda$mirna_name))
keep.mir.mtb.mda <- intersect(unique(mtb$miRNA),unique(mda$mirna_name))
keep.mir.mmap.mda <- intersect(unique(mirmap.unique$mature_name),unique(mda$mirna_name))

# mirDB intersect
keep.mir.mtb.mdb <- intersect(unique(mtb$miRNA),unique(mdb.filt$mirid))
keep.mir.mmap.mdb <- intersect(unique(mirmap.unique$mature_name),unique(mdb.filt$mirid))

# mirMap intersect
keep.mtb.mmap <- intersect(unique(mtb$miRNA),unique(mirmap.unique$mature_name))

keep.mir <- unique(c(keep.mir.mda.ts,keep.mir.mdb.ts,keep.mir.mtb.ts,keep.mir.mmap.ts,
                   keep.mir.mdb.mda,keep.mir.mtb.mda,keep.mir.mmap.mda,
                   keep.mir.mtb.mdb,keep.mir.mmap.mdb,
                   keep.mtb.mmap)) # 2035 miR IDs

#==============================
# 
# common miRs consensus matrix
#
#==============================

ts.df <- ts[ts[,5] %in% keep.mir,]; 
mda.df <- mda[mda$mirna_name %in% keep.mir,]; 
mdb.df <- mdb.filt[mdb.filt$mirid %in% keep.mir,]; 
mtb.df <- mtb[mtb$miRNA %in% keep.mir,]; 
mmap.df <- mirmap.unique[mirmap.unique$mature_name %in% keep.mir,]

cmir <- matrix(nrow=length(keep.mir),ncol=6)
colnames(cmir) <- c("targetscan","miranda","mirdb","mirtarbase","mirmap","total")
rownames(cmir) <- keep.mir

cmir[,1] <- ifelse(rownames(cmir) %in% ts.df[,5],1,0)
cmir[,2] <- ifelse(rownames(cmir) %in% mda.df$mirna_name,1,0)
cmir[,3] <- ifelse(rownames(cmir) %in% mdb.df$mirid,1,0)
cmir[,4] <- ifelse(rownames(cmir) %in% mtb.df$miRNA,1,0)
cmir[,5] <- ifelse(rownames(cmir) %in% mmap.df$mature_name,1,0)
cmir[,6] <- apply(cmir[,1:5],1,sum)

save(cmir,file="df-consensus-miR_5db_ts-mda-mdb-mtb-mmap.rda")


#==================================
# 
#   Common miR Target Interactions
#
#===================================

tsx <- paste0(ts.df$V5,";",ts.df$V2)
mdax <- paste0(mda.df$mirna_name,";",mda.df$gene_symbol)
mdbx <- paste0(mdb.df$mirid,";",mdb.df$target.symbol)
mtbx <- paste0(mtb.df$miRNA,";",mtb.df$Target.Gene)
mmapx <- paste0(mmap.df$mature_name,";",mmap.df$gene_name)

# pairwise consensus interactions
pci.mdax.tsx <- intersect(tsx,mdax)
pci.mdbx.tsx <- intersect(tsx,mdbx)
pci.mtbx.tsx <- intersect(tsx,mtbx)
pci.mmapx.tsx <- intersect(tsx,mmapx)

pci.mdbx.mdax <- intersect(mdax,mdbx)
pci.mtbx.mdax <- intersect(mdax,mtbx)
pci.mmapx.mdax <- intersect(mdax,mmapx)

pci.mtbx.mdbx <- intersect(mdbx,mtbx)
pci.mmapx.mdbx <- intersect(mdbx,mmapx)

pci.mtbx.mmapx <- intersect(mtbx,mmapx)

pci.keep <- unique(c(pci.mdax.tsx,pci.mdbx.tsx,pci.mtbx.tsx,pci.mmapx.tsx,
                     pci.mdbx.mdax,pci.mtbx.mdax,pci.mmapx.mdax,
                     pci.mtbx.mdbx,pci.mmapx.mdbx,
                     pci.mtbx.mmapx))

pci.df <- as.data.frame(matrix(nrow=length(pci.keep),ncol=8))
colnames(pci.df) <- c("miR","target.symbol","targetscan","miranda","mirdb","mirtarbase","mirmap","total")
pci.df[,1] <- gsub(";.*","",pci.keep)
pci.df[,2] <- gsub(".*;","",pci.keep)
rownames(pci.df) <- pci.keep

pci.df[,3] <- ifelse(rownames(pci.df) %in% tsx,1,0)
pci.df[,4] <- ifelse(rownames(pci.df) %in% mdax,1,0)
pci.df[,5] <- ifelse(rownames(pci.df) %in% mdbx,1,0)
pci.df[,6] <- ifelse(rownames(pci.df) %in% mtbx,1,0)
pci.df[,7] <- ifelse(rownames(pci.df) %in% mmapx,1,0)
pci.df[,8] <- apply(pci.df[,3:7],1,sum)

#==============
#
#   Citations
#
#==============

# 1. miRmap
# Charles E. Vejnar and Evgeny M. Zdobnov
# miRmap: Comprehensive prediction of microRNA target repression strength
# Nucleic Acids Research 2012 Dec 1;40(22):11673-83. doi: 10.1093/nar/gks901

# 2. TargetScan
# Vikram Agarwal George W Bell Jin-Wu Nam David P Bartel
# Predicting effective microRNA target sites in mammalian mRNAs
# eLife 2015;4:e05005

# 3a. miRDB 1 - database
# Nathan Wong and Xiaowei Wang (2015) miRDB: an online resource for microRNA target prediction and functional annotations. Nucleic Acids Research. 43(D1):D146-152.

# 3b. miRDB 2 - target finding algorithm
# Xiaowei Wang (2016) Improving microRNA target prediction by modeling with unambiguously identified microRNA-target pairs from CLIP-Ligation studies. Bioinformatics. 32(9):1316-1322.  

# 4. miRTarBase
# miRTarBase 2016: updates to the experimentally validated miRNA-target interactions database.
# Chou CH, Chang NW, Shrestha S, Hsu SD, Lin YL, Lee WH, Yang CD, Hong HC, Wei TY, Tu SJ, Tsai TR, Ho SY, Jian TY, Wu HY, Chen PR, Lin NC, Huang HT, Yang TL, Pai CY, Tai CS, Chen WL, Huang CY, Liu CC, Weng SL, Liao KW, Hsu WL, Huang HD. 
# (2016) Nucleic acids research.

# 5. Miranda target predictions
# The microRNA.org resource: targets and expression. 
# Betel D, Wilson M, Gabow A, Marks DS, Sander C., 
# Nucleic Acids Res. 2008 Jan; 36(Database Issue): D149-53. 

```

