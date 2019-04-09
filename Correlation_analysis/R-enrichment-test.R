##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##                       enrichment test of gene lists
##
##  @Bin Zhang, 2005 - 2016
## 
##  Input: tab-delimited files, in which the last column indocates groups or categories
##  Output: pairwised overlap analysis
##          intersection of all files
##               overlapped genes
##               overlapped number and significance
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#gene number
#16387 rosmap
#23201 bm
#39579 hb

library(class)
library(rpart)

#setwd("/sc/orga/projects/zhangb03a/lei_guo/YWHAZ/Corr_Results/")
source("R-functions.R")
  
  filelist = c("/sc/orga/projects/zhangb03a/lei_guo/YWHAZ/Corr_Results/YWHAZ_KD_KO_DEG.tsv", "/sc/orga/projects/zhangb03a/lei_guo/YWHAZ/Corr_Results/Genes_sig_cor_YWHAZ_consensus.tsv")

  fileshortnames = c("DEG", "Dataset")
  muniqueIdCol   = c(1, 1);
  mygenesymbolIdx= NULL # gene symbol column is different from "muniqueIdCol"

  ioutputDir="Module_Overlap/"
  nfiles    = length(filelist);

  mygeneInforCols=2

  ijtotals = matrix(48300, nfiles, nfiles);


ioutputDir="Module_Overlap/"
if(ioutputDir != ""){
  dir.create(ioutputDir)
}

###### pair-wised comparison ##########
#
no.files=length(filelist)

for (i in c(1:(no.files-1)) ){
    for (j in c(no.files) ){

        ijkeywords=paste(fileshortnames[i], "-vs-", fileshortnames[j], sep="")

        print(ijkeywords)
        moduleBasedIntersectionMatrix(filelist[i], filelist[j], 
                                outputDir=ioutputDir, keywords=ijkeywords, 
                                uniqueIdCol=muniqueIdCol[c(i,j)], 
                                geneInforCols=mygeneInforCols,
                                genesymbolIdx=mygenesymbolIdx, 
                                itotalGenes = ijtotals[i,j], signifpvalue=2, latex=F, removeDuplicate=TRUE, 
                                restrict_networkB_to_A=F)


    }
}



############################## END #################################################
