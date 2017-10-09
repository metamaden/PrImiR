# PrImiR
## Description
Database of pairwise consensuses for computationally predicted miR-mRNA target interactions.

## Background
This package is downloadable using the following in R:
```require(devtools); install_github("metamaden/PrImiR")```

The science of miR target prediction is actively evolving. There is no universal consensus for all computationally predicted miR-mRNA target interactions, and in addition most computationally predicted interactions have not been validated in the wet lab. 

To mitigate the limitations of individual computational algorithms, this package provides a database of consensus predicted miRs (mature miRNAs) and miR-mRNA target interactions (present in at least two databases) using the latest available versions of five established and widely used miRNA databases and target prediction algorithms:

1. [Target Scan (v.7.1)](http://www.targetscan.org/vert_71/)
2. [Miranda/mirSVR (update 2010-08)](http://www.microrna.org/microrna/getDownloads.do)
3. [miRDB (v.5.0)](http://www.mirdb.org/download.html)
4. [miRTarBase (v.6.1)](http://mirtarbase.mbc.nctu.edu.tw/php/download.php)
5. [miRmap (update 2013-01-09)](http://mirmap.ezlab.org/)

#### 1. Database Consensus Visualization: Simple Heatmap of random subset (N=2k miR-mRNA interactions) of PrImiR pairwise consensus miR-mRNA interactions database, showing presence of miR-mRNA interactions (rows) across databases (columns). blue=0/not present; red=1/present

![PrImiR pairwise consensus heatmap summary](https://github.com/metamaden/PrImiR/blob/master/hm-consensus_interactions_2k.jpeg "PrImiR Interactions")

#### 2. Interactions (miR-mRNA targets) Predicted across databases

![PrImiR Venn Diagram Consensus of Predicted Interactions Across Five Databases](https://github.com/metamaden/PrImiR/blob/master/venn-quintuple_PrImiR-db.jpg "PrImiR Predicted Interaction Consensus Across Five Databases")

#### 3. Scripts for Assembling And Annotating miRNA/miR Database
Check scripts folder for information. Check data folder for database and annotation files, including CpGs overlapping putative promoter regions (+/- 1KB from hairpin start site).

## Citations
Agarwal, Vikram, George W Bell, Jin-Wu Nam, David P Bartel. _Predicting effective microRNA target sites in mammalian mRNAs._ eLife 2015;4:e05005

Betel D, Wilson M, Gabow A, Marks DS, Sander C. _The microRNA.org resource: targets and expression._ Nucleic Acids Res. 2008 Jan; 36(Database Issue): D149-53. 

Chou CH, Chang NW, Shrestha S, et. al. _miRTarBase 2016: updates to the experimentally validated miRNA-target interactions database._ (2016) Nucleic acids research.

Vejnar, Charles E. and Evgeny M. Zdobnov. _miRmap: Comprehensive prediction of microRNA target repression strength._ Nucleic Acids Research 2012 Dec 1;40(22):11673-83. doi: 10.1093/nar/gks901

Wang, Xiaowei _Improving microRNA target prediction by modeling with unambiguously identified microRNA-target pairs from CLIP-Ligation studies_. 2016 Bioinformatics. 32(9):1316-1322.  

Wong, Nathan and Xiaowei Wang. _miRDB: an online resource for microRNA target prediction and functional annotations._ 2015 Nucleic Acids Research. 43(D1):D146-152.

