return_original=function(vect, topN=10){
  if(length(vect) <topN){ return(vect)
  } else {
    return(vect[1:topN])
  }
}


getAllParts=function(fullfnames, sep="-", max_len=NULL, ret_len=F){

  splitted=unlist( strsplit(fullfnames[1], sep) )
  nn  = length(fullfnames)
  nf  = length(splitted)
  if(!is.null(max_len)){nf=max(nf, max_len)}

  ret = matrix("", nn, nf)
  lens= rep(0, nn)
  for(i in c(1:nn) ) {
    each = fullfnames[i]
    splitted=unlist( strsplit(each, sep) )
    ino     = length(splitted)
    if(ino >=nf) {
       ret[i,] = splitted[1:nf]
    }else{
       ret[i,] = c(splitted, rep("",nf-ino ))
    }

    lens[i] = ino
  }

  if ( ret_len ){
     return (lens)
  } else {
     return (ret) 
  }
  
}


concatenateSelectedElelments=function(boolvect, mnames, sep="-"){
    return(concatenate(mnames[boolvect], mysep=sep))
}


concatenateSelectedElelmentNames = function(vect, vnames, selected_value, mysep="_")
{
   sel = vect==selected_value
   if(sum(sel)==1){return(vnames[sel]) }

   res = concatenate(vnames[sel], mysep=mysep)   
   res
}


concatenateSelectedElelments2vect = function(boolvect, boolvect2, mnames, sep="-"){
    return(concatenate(mnames[boolvect&boolvect2], mysep=sep))
}

concatenateSelectedElelments_vectMatrix =function(boolvect, boolMatrix, mnames, sep="-"){
    res = apply(boolMatrix, 2, concatenateSelectedElelments2vect, boolvect2=boolvect, mnames=mnames, sep=sep)
    return (res)
}

getAllPartsVect = function(fullfname, sep="-"){

  splitted=unlist( strsplit(fullfname, sep) )
  return (splitted) 
   
}


# Test the modules in the files for enrichment of gene ontology terms
moduleBasedIntersectionMatrix_GeneOntology = function(fnames, fnameB, outputDir="", uniqueIdCol=1, signifLevel=1, 
   removeDuplicate=TRUE, removeGrey=TRUE,totalbackground=NULL)

{
    ########################### GO Matrix ##################################################
    #
    #System	GeneCategory	PopulationHits	PopulationTotal	GeneSymbol
    #GO Biological Process	cell-cell signaling	933	16773	NAALAD2; NAALADL1; ...
    #GO Biological Process	protein biosynthesis	912	16773	PIGK; FARSLB; NR1H3; ...
    #
    bbMatrixAll <- read.delim(fnameB, sep="\t", header=T)
    bbMatrixAll <- as.matrix(bbMatrixAll)
    dim(bbMatrixAll)

    print("Build up GO term membership matrix ....")

    no.colsB = ncol(bbMatrixAll)
    newIDs      <- paste(bbMatrixAll[, 1], bbMatrixAll[, 2], sep=";;;")
    termGenes   <- tapply(bbMatrixAll[, no.colsB], INDEX=newIDs, FUN=getAllPartsVect, sep="; ")
    termGenes   <- lapply(termGenes, unique)
    termGenesLen<- unlist(lapply(termGenes, length))

    if(is.null(totalbackground)) totalbackground = unique(unlist(termGenes))
    no.allnodes     = length(totalbackground)
    allnodesIndex = c(1:no.allnodes); names(allnodesIndex) = totalbackground 

    empty    = rep(FALSE, no.allnodes)

    mres     = lapply(termGenes, fillInVector, indexvector=allnodesIndex, emptyvector=empty, allnames=totalbackground)
    datMtrxB = data.frame(do.call(rbind, mres))
    datMtrxB = t(datMtrxB)
    dim(datMtrxB)

    modlevelsB= names(termGenes)
    no.modulesB = length(modlevelsB)
    sizes.modulesB = apply(datMtrxB, 2, sum)

    #---------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------------
    #
    for(fnameA in fnames) {

    print(fnameA)

    fname         = getFileName(fnameA)
    outfnameTop   = paste(outputDir,  fname, "_Ontology.xls", sep='')

    aaMatrixAll <- read.delim(fnameA, sep="\t", header=T)
    aaMatrixAll <- as.matrix(aaMatrixAll)
    dim(aaMatrixAll)
    modProbeSizes <- tapply(aaMatrixAll[,1], INDEX=aaMatrixAll[,ncol(aaMatrixAll)], FUN=length)

    if(removeGrey) { # remove grey module
        gsel = aaMatrixAll[, ncol(aaMatrixAll)] != "grey"
        aaMatrixAll = aaMatrixAll[gsel,]
    }

    #aaMatrixAll[, uniqueIdCol[1]] = toupper(aaMatrixAll[, uniqueIdCol[1]])

    no.colsA = dim(aaMatrixAll)[2]

    # include gene infor in the aaMatrix but not bbMatrix, in this way the merged matrix
    #  will not have duplicated columns
    #
    aaMatrix <- aaMatrixAll[, c(uniqueIdCol[1], no.colsA)]
    colnames(aaMatrix) <- colnames(aaMatrixAll)[c(uniqueIdCol,no.colsA)]


    #------ remove duplicates in each module -------------------
    #
    imoduleColA=dim(aaMatrix)[2]
    strA = paste(aaMatrix[,1], aaMatrix[,imoduleColA])
    idxA = findUniqueIndex(strA)
    if(length(idxA)<nrow(aaMatrix)) {
       aaMatrix = aaMatrix[idxA, ]
    }
    
    #----------------- module-based overlap -------------------------
    #
    # restrict A to B (go terms)
    #
    restricted     = setElementInSet_Fast(aaMatrix[,1], totalbackground)
    aaMatrixRstrct = aaMatrix[restricted, ]

    modlevelsA= sort(unique(aaMatrixRstrct[,imoduleColA]) )
    no.modulesA = length(modlevelsA)

    #-------------------make probe-categories matrix -----------------
    #
    # TRUE means a probe belongs to the corresponding category
    #
    ncA = ncol(aaMatrixRstrct); 

    print("Build up module membership matrix ....")
    mres     = tapply(aaMatrixRstrct[, 1], INDEX=aaMatrixRstrct[, ncA], fillInVector, indexvector=allnodesIndex, emptyvector=empty, allnames=totalbackground )
    datMtrxA = data.frame(do.call(rbind, mres))
    datMtrxA = t(datMtrxA)
    dim(datMtrxA)

    # find out module sizes in A, B, B restricted to A 
    sizes.modulesA = apply(datMtrxA, 2, sum)
    sum(termGenesLen!= sizes.modulesB)

    #****************** overlap **************************************
    print("Perform FET test ....")

    ovlpMtrx = t(datMtrxA) %*% datMtrxB
    
    # match column names and rownames to the A module names and B module names
    midx           = getMatchedIndexFast(names(sizes.modulesA), rownames(ovlpMtrx))
    sizes.modulesA = sizes.modulesA[midx]    
    midx           = getMatchedIndexFast(names(sizes.modulesB), colnames(ovlpMtrx))
    sizes.modulesB = sizes.modulesB[midx]

    # prepare FET test vectors
    idxy       = getMatrixIndex(size=dim(ovlpMtrx), symmetric=FALSE, diagonal=TRUE)
    no.tests   = nrow(idxy)
    fetMtrx    = cbind(rep(no.allnodes, no.tests), sizes.modulesA[idxy[,1]], sizes.modulesB[idxy[,2]],  ovlpMtrx[idxy])
    foldEnrich = (fetMtrx[,4]/fetMtrx[,3])*(fetMtrx[,1]/fetMtrx[,2])

    # turn vector into matrix
    foldEnrichMtrx = matrix(foldEnrich, nrow=no.modulesA)

    #*************** fisher's test *********************************
    #
    fetP     = apply(fetMtrx, 1, fisherTest, minPopulationHits=0)
    od       = order(fetP)

    # turn vector into matrix, be careful of the the order of vector into matrix
    fetPMtrx = t(matrix(fetP,  nrow=no.modulesB))

    # output pair
    correction    = no.modulesA * no.modulesB
    fetPcorrected = fetP*correction
    fetPcorrected = ifelse(fetPcorrected >1, 1, fetPcorrected)
    fetPcorrected = signif(fetPcorrected, 2)

    #-------------------------------------------------------------------------------------
    #------------ find overlap genes for each group pair          ------------------------
    #

    # to pull out only the overlaps from the significant tests
    selMtrx = fetPMtrx <= signifLevel
    selMod  = (apply(selMtrx, 1, sum))>=1;
    selGot  = (apply(selMtrx, 2, sum))>=1;

    print(sprintf("Get overlap genes for %s modules and %s terms ....", sum(selMod), sum(selGot) ) )

    #tmp=concatenateSelectedElelments_vectMatrix(datMtrxA[,178], datMtrxB, mnames=totalbackground, sep=",")
    #mergedMembersGS = apply(datMtrxA, 2, concatenateSelectedElelments_vectMatrix, boolMatrix=datMtrxB, mnames=totalbackground, sep=",")

    mergedMembersGS = apply(datMtrxA[,selMod,drop=FALSE], 2, concatenateSelectedElelments_vectMatrix, boolMatrix=datMtrxB[,selGot], mnames=totalbackground, sep=",")
    mergedMembersGS = t(mergedMembersGS )

    idxy2      = getMatrixIndex(size=dim(mergedMembersGS), symmetric=FALSE, diagonal=TRUE)
    membersVect= mergedMembersGS[idxy2]
    names.membersVect = paste(rownames(mergedMembersGS)[idxy2[,1]], colnames(mergedMembersGS)[idxy2[,2]])

    # coimbine results from all the pair-wise tests
    xmodnames  = names(sizes.modulesA)[idxy[,1]]
    xmodSizesP = modProbeSizes[xmodnames]
    goSrcTerms = getAllParts(names(sizes.modulesB), ";;;")
    xgoSrcTerms= goSrcTerms[idxy[,2],]
    xgoTermSzs = cbind(sizes.modulesB[idxy[,2]], rep(no.allnodes, no.tests) )

    # matched the significant pairs to the all apirs
    names.alltests = paste(xmodnames, names(sizes.modulesB)[idxy[,2]])
    midx = getMatchedIndexFast(names.alltests, names.membersVect)
    sum(is.na(midx))
    membersVect2 = rep("", no.tests)
    membersVect2[midx] = membersVect

    ocompares = cbind(xmodnames,xmodSizesP, xgoSrcTerms)
    opairMtrx = cbind(ocompares, ovlpMtrx[idxy], sizes.modulesA[idxy[,1]], xgoTermSzs, fetP,  fetPcorrected, foldEnrich,  membersVect2)
    colnames(opairMtrx) = c("Module", "OriginalModuleSize", "System", "Gene.Category", "OverlapSize", "ModuleSize", 
                            "CategoryTermSize", "PopulationSize", "FET_P", "Corrected_P", "Fold_Enrichment", "Genes")
    opairMtrx = opairMtrx[od,]

    pidx = getMatchedIndex(colnames(opairMtrx), "FET_P")
    tmpP = as.numeric(opairMtrx[, pidx]);     
    xsel = tmpP<signifLevel
    opairMtrx2 = opairMtrx[xsel,]
    write.table(rbind(colnames(opairMtrx2)), outfnameTop, sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(opairMtrx2,                  outfnameTop, sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
    }

}



# assume that subvect is a subset of cvector
getMatchedIndex=function(cvector, subvect){
  fullindex = c(1:length(cvector) )
  orgIdx    = cbind(cvector, fullindex)

  index2    = c(1:length(subvect))
  subIdex   = cbind(subvect, index2)

  merged    = merge(subIdex, orgIdx, by.x=1, by.y=1, all.x=T)
  if( dim(merged)[1]==0 ) {return(NULL)}

  merged    = as.matrix(merged)

  if(dim(merged)[1]>1){
    od        = order(as.integer(merged[,2]))  # restore the original order of subvect
    merged    = merged[od, ]
  }
  
  outIndex  = as.integer(merged[,3])

  return (outIndex)
}


# assume that subvect is a subset of cvector
getMatchedIndexFast=function(cvector, subvect){
  fullindex = c(1:length(cvector) )
  orgIdx    = cbind(cvector, fullindex)

  index2    = c(1:length(subvect))
  subIdex   = cbind(subvect, index2)

  merged    = merge(subIdex, orgIdx, by.x=1, by.y=1, all.x=T)
  if( dim(merged)[1]==0 ) {return(NULL)}

  merged    = as.matrix(merged)

  if(dim(merged)[1]>1){
    od        = order(as.integer(merged[,2]))  # restore the original order of subvect
    merged    = merged[od, ]
  }
  
  outIndex  = as.integer(merged[,3])

  return (outIndex)
}



setElementInSet_Fast = function(setA,setB){
    found = rep(FALSE, length(setA))
    axed = cbind(setA, c(1:length(setA) ))
    merged=merge(axed, cbind(setB), by.x=1, by.y=1, all=F)
    if (dim(merged)[1]==0) {return(found) } 
    merged=as.matrix(merged)
   
    found[as.integer(merged[,2])] = TRUE

    return (found)
}



# for unsorted vectors
findUniqueIdx=function(vect){

   no.elements = length(vect)
   if(no.elements<=1){return (1)}

   vectIdxed   = cbind(vect, c(1:no.elements) )
   mo = order(vectIdxed[,1])
   merged = vectIdxed[mo,]

   bounds = findBoundary(merged[,1])

   uniqueIdx = as.integer(merged[,2])[bounds]

   # keep original order
   return ( sort(uniqueIdx) )
}

findUniqueIndex=function(vect){
 res = findUniqueIdx(vect)
 return (res)
}

# a= 1 1 2 2 3 3 3
# findBoundary(a) => [1] 1 3 5
#
findBoundary=function(vect){
   no.elements = length(vect)
   esel = vect[no.elements]==vect[1]
   esel = ifelse(is.na(esel), F, esel)

   #if(esel){
   #  return (1);
   #}

   shifted= c(vect[no.elements], vect[1:(no.elements-1)] )
   sel = shifted != vect

   # handle NA
   sel = ifelse(is.na(sel), F, sel)

   sel[1] = TRUE

   return ( c(1:no.elements)[sel] )
}



# to prevent genes that are not in allnames, we use intersection of input names and all the names 
#
fillInVector = function(mynames, indexvector, emptyvector, allnames)
{
   mynamesFound = intersect(mynames, allnames)
   iidx = indexvector[mynamesFound]
   ovect= emptyvector  
   ovect[iidx] = TRUE
   return (ovect)
}


#************ example, Fisher exact Test & Hypergeometric test ********************
#
# 196 of 673 genes overlap the 948 genes, total geens 24000
#a=rbind( c(673-196, 24000-(673-196)-(948-196)-196), 
#          c(196,     948-196))
#
#  signif(fisher.test(a)$p.value,2)
# hypergeometric test: 
#    phyper(196, 673, 24000-673, 948, lower.tail=F)
#
# By chance you expect 196/23756*948 = 8 overlap
#
# if population hits < minPopulationHits, then the test doesn't continue
#
#input vector:=[population, population_hit, list, list_hit]
fisherTest = function(poplisthits, minPopulationHits=5)
{
   if ( (poplisthits[2]<minPopulationHits) | (poplisthits[4]==0) )
      return (1)

   q= poplisthits[4] #list hit
   m= poplisthits[2] #population hit
   n= poplisthits[1]-poplisthits[2] #population non-hit
   k= poplisthits[3] #list total

   myp=phyper(q-1, m, n, k, lower.tail=F)
   signif(myp,3)
}


# restrict_networkB_to_A ==T: take networkB as background when mapping modulesB to network A
# restrict_networkB_to_A ==F: consider modulesB as independent, so map individual modules to network A
#
moduleBasedIntersectionMatrix=function(fnameA, fnameB,outputDir="", keywords="", uniqueIdCol=c(1,1), geneInforCols=8, genesymbolIdx=NULL, itotalGenes = 3600, signifpvalue=0.001, latex=F, removeDuplicate=F, restrict_networkB_to_A=T)
{
    #fkeyB       =getFileName(fnameB)
    #outfname    paste(outputDir, "intersect.csv", sep='')

    # get unique labels for each input files
    keys= unlist(strsplit(keywords, "-vs-"))
    fkeyA       =keys[1]
    fkeyB       =keys[2]

    #outfname     = paste(outputDir,  "intersect_", keywords, ".xls", sep='')
    outfname     = paste(outputDir,  "intersect_", keywords, ".xls", sep='')

    nofname     = paste(outputDir,  "intersect_", keywords, "_number.xls", sep='')
    pvalfname   = paste(outputDir,  "intersect_", keywords, "_pvalue.xls", sep='')
    freqfname   = paste(outputDir,  "intersect_", keywords, "_frequence.xls", sep='')

    pairfname   = paste(outputDir,  "intersect_", keywords, "_pair.xls", sep='')
    pairfnameR  = paste(outputDir,  "intersect_", keywords, "_pair-restricted.xls", sep='')


    latexnofname  = paste(outputDir, "intersect_", keywords, "_number.tex", sep='')
    latexpvalfname= paste(outputDir, "intersect_", keywords, "_pvalue.tex", sep='')

    aaMatrixAll <- read.delim(fnameA, sep="\t", header=T)
    dim(aaMatrixAll)
    aaMatrixAll <- as.matrix(aaMatrixAll)

    bbMatrixAll <- read.delim(fnameB, sep="\t", header=T)
    dim(bbMatrixAll)
    bbMatrixAll <- as.matrix(bbMatrixAll)

    aaMatrixAll[, uniqueIdCol[1]] = toupper(aaMatrixAll[, uniqueIdCol[1]])
    bbMatrixAll[, uniqueIdCol[2]] = toupper(bbMatrixAll[, uniqueIdCol[2]])

    no.colsA = dim(aaMatrixAll)[2]
    no.colsB = dim(bbMatrixAll)[2]

    # include gene infor in the aaMatrix but not bbMatrix, in this way the merged matrix
    #  will not have duplicated columns
    #
    aaMatrix <- aaMatrixAll[, c(uniqueIdCol[1], 1:geneInforCols, no.colsA)]
    colnames(aaMatrix) <- c( colnames(aaMatrixAll)[c(uniqueIdCol[1], 1:geneInforCols)], fkeyA)

    bbMatrix <- bbMatrixAll[, c(uniqueIdCol[2], no.colsB)]
    colnames(bbMatrix) <- c( colnames(bbMatrixAll)[1], fkeyB)

    imoduleColA=dim(aaMatrix)[2]
    imoduleColB=dim(bbMatrix)[2]

    # remove duplicates
    strA = paste(aaMatrix[,1], aaMatrix[,imoduleColA])
    idxA = findUniqueIndex(strA)
    if(length(idxA)<nrow(aaMatrix)) {
       aaMatrix = aaMatrix[idxA, ]
    }

    strB = paste(bbMatrix[,1], bbMatrix[,imoduleColB])
    idxB = findUniqueIndex(strB)
    if(length(idxB)<nrow(bbMatrix)) {
      bbMatrix = bbMatrix[idxB, ]
    }

    #--------------  total overlap --------------------
    AA=unique(as.character(aaMatrix[, 1]))
    BB=unique(as.character(bbMatrix[, 1]))
    #totaloverlap    = merge(AA,BB, by.x=1, by.y=1, sort=F, all=F)
    #no.totalOverlap = dim(totaloverlap)[1]
    totaloverlap    = intersect(AA, BB)
    no.totalOverlap = length(totaloverlap)

    totalbackground = length(union(AA, BB))
    if (!is.na(itotalGenes)){
        allnodes = sort(union(AA, BB))
        totalbackground = itotalGenes
    }else{
        totalbackground = length(AA)
        allnodes = sort(AA)
    }
    no.allnodes = length(allnodes)

    # global enrichment test
    #
    pval = phyper(no.totalOverlap-1, length(AA), itotalGenes-length(AA), length(BB), lower.tail=F)
    globalEnrich = c(no.totalOverlap , pval, length(AA), length(BB))
    
    # matched to aaMatrix (more conservative test)
    #
    ccMatrix <- merge(bbMatrix, cbind(totaloverlap), by.x=1, by.y=1, all=F)
    ccMatrix <- as.matrix(ccMatrix)

    #----------------- module-based overlap -------------------------
    #
    # restrict B to A
    restricted     = setElementInSet_Fast(bbMatrix[,1], allnodes)
    bbMatrixRstrct = bbMatrix[restricted, ]

    modlevelsA= names(table(aaMatrix[,imoduleColA]) )
    modlevelsB= names(table(bbMatrix[,imoduleColB]) )
    no.modulesA = length(modlevelsA)
    no.modulesB = length(modlevelsB)

    #-------------------make probe-categories matrix -----------------
    #
    # TRUE means a probe belongs to the corresponding category
    #
    allnodesIndex = c(1:no.allnodes); names(allnodesIndex) = allnodes
    ncA = ncol(aaMatrix); ncB = ncol(bbMatrix)

    empty    = rep(FALSE, no.allnodes)
    mres     = tapply(aaMatrix[, 1], INDEX=aaMatrix[, ncA], fillInVector, indexvector=allnodesIndex, emptyvector=empty, allnames=allnodes)
    datMtrxA = data.frame(do.call(rbind, mres))
    datMtrxA = t(datMtrxA)
    dim(datMtrxA)

    mres     = tapply(bbMatrix[, 1], INDEX=bbMatrix[, ncB], fillInVector, indexvector=allnodesIndex, emptyvector=empty, allnames=allnodes)
    datMtrxB = data.frame(do.call(rbind, mres))
    datMtrxB = t(datMtrxB)
    dim(datMtrxB)

    # find out module sizes in A, B, B restricted to A 
    sizes.modulesA = apply(datMtrxA, 2, sum)
    sizes.modulesBr= apply(datMtrxB, 2, sum)
    sizes.modulesB = table(bbMatrix[,imoduleColB])

    #****************** overlap **************************************
    ovlpMtrx = t(datMtrxA) %*% datMtrxB
    
    # match column names and rownames to the A module names and B module names
    midx = getMatchedIndex(names(sizes.modulesA), rownames(ovlpMtrx))
    sizes.modulesA = sizes.modulesA[midx]    
    midx = getMatchedIndex(names(sizes.modulesB), colnames(ovlpMtrx))
    sizes.modulesB = sizes.modulesB[midx]
    midx = getMatchedIndex(names(sizes.modulesBr), colnames(ovlpMtrx))
    sizes.modulesBr = sizes.modulesBr[midx]

    # prepare FET test vectors
    idxy = getMatrixIndex(size=dim(ovlpMtrx), symmetric=FALSE, diagonal=TRUE)
    no.tests = nrow(idxy)
    fetMtrx  = cbind(rep(totalbackground, no.tests), sizes.modulesA[idxy[,1]], sizes.modulesB[idxy[,2]],  ovlpMtrx[idxy])
    fetMtrxR = cbind(rep(no.allnodes,     no.tests), sizes.modulesA[idxy[,1]], sizes.modulesBr[idxy[,2]], ovlpMtrx[idxy])

    foldEnrich = (fetMtrx[,4]/fetMtrx[,3])*(fetMtrx[,1]/fetMtrx[,2])
    foldEnrichR= (fetMtrxR[,4]/fetMtrxR[,3])*(fetMtrxR[,1]/fetMtrxR[,2])

    # turn vector into matrix
    foldEnrichMtrx = matrix(foldEnrich, nrow=no.modulesA)
    foldEnrichMtrxR= matrix(foldEnrichR, nrow=no.modulesA)

    #************ fisher's test *********************************
    #
    fetP     = apply(fetMtrx,  1, fisherTest, minPopulationHits=0)
    fetPR    = apply(fetMtrxR, 1, fisherTest, minPopulationHits=0)

    # turn vector into matrix, be careful of the the order of vector into matrix
    fetPMtrx = t(matrix(fetP,  nrow=no.modulesB))
    fetPMtrxR= t(matrix(fetPR, nrow=no.modulesB))

    # output pair
    correction    = length(setdiff(names(sizes.modulesA), "grey")) *length(setdiff(names(sizes.modulesB), "grey")) 
    fetPcorrected =  fetP*correction; fetPcorrected = ifelse(fetPcorrected >1, 1, fetPcorrected)
    fetPcorrectedR= fetPR*correction; fetPcorrectedR= ifelse(fetPcorrectedR>1, 1, fetPcorrectedR)
    od = order(fetP)

    ocompares = paste(names(sizes.modulesA)[idxy[,1]], names(sizes.modulesB)[idxy[,2]])
    opairMtrx = cbind(ocompares, names(sizes.modulesA)[idxy[,1]], names(sizes.modulesB)[idxy[,2]],  fetMtrx,  foldEnrich,  fetP,  fetPcorrected )[od,]
    opairMtrxR= cbind(ocompares, names(sizes.modulesA)[idxy[,1]], names(sizes.modulesBr)[idxy[,2]], fetMtrxR, foldEnrichR, fetPR, fetPcorrectedR)[od,]
    tmpStrs = paste("module_size (", c(fkeyA, fkeyB), ")", sep="")
    colnames(opairMtrx) = c("Comparison", fkeyA, fkeyB,"population", tmpStrs, "Overlap", "Fold_Enrichment", "FET_P", "FET_P_corrected")
    colnames(opairMtrxR)= c("Comparison", fkeyA, fkeyB,"population", tmpStrs, "Overlap", "Fold_Enrichment", "FET_P", "FET_P_corrected")

    # -------- global FET test ----------------------
    globalFetMtrx = cbind(rep(totalbackground, no.modulesB), rep(length(AA)), 
                      rep(length(BB), no.modulesB),  sizes.modulesB)
    fetPglobal    = apply(globalFetMtrx, 1, fisherTest, minPopulationHits=0)
    foldEglobal   = (globalFetMtrx[,4]/globalFetMtrx[,3])*(globalFetMtrx[,1]/globalFetMtrx[,2])
    freqGlobal    = globalFetMtrx[,3]/totalbackground
    freqGlobalR   = globalFetMtrx[,3]/no.allnodes

    #-------------------------------------------------------------------------------------
    #------------ output AxB matrix results ----------------------------------------------
    #
    oNumberMtrx = rbind(sizes.modulesB, sizes.modulesBr, ovlpMtrx)
    oNumberMtrx = cbind(c("total", "All", names(sizes.modulesA)),
                        c(totalbackground, no.allnodes, sizes.modulesA), oNumberMtrx)
    colnames(oNumberMtrx) = c("overlap", "total", names(sizes.modulesB))

    #--- FET fold enrichment
    oFoldMtrx = rbind(freqGlobal, foldEnrichMtrx)
    oFoldMtrx = cbind(c("ALL_freq", names(sizes.modulesA)),oFoldMtrx)
    colnames(oFoldMtrx) = c("freq/fold_enrichment", names(sizes.modulesB))

    oFoldMtrxR= rbind(freqGlobalR, foldEnrichMtrxR)
    oFoldMtrxR= cbind(c("ALL_freq", names(sizes.modulesA)),oFoldMtrxR)
    colnames(oFoldMtrxR) = c("freq/fold_enrichment", names(sizes.modulesB))

    #--- FET p
    oFetPMtrx = rbind(fetPglobal, fetPMtrx)
    oFetPMtrx = cbind(c("All", names(sizes.modulesA)),oFetPMtrx)
    colnames(oFetPMtrx) = c("pvalue", names(sizes.modulesB))

    oFetPMtrxR= fetPMtrxR
    oFetPMtrxR= cbind(names(sizes.modulesA),oFetPMtrxR)
    colnames(oFetPMtrxR) = c("pvalue", names(sizes.modulesBr))

    #-------------------------------------------------------------------------------------
    #------------ save Matrix form results into files ------------------------------------
    #
    write.table(rbind(colnames(oNumberMtrx)), nofname,   sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(oNumberMtrx, nofname,   sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)

    write.table(rbind(colnames(oFetPMtrx)), pvalfname, sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(oFetPMtrx, pvalfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)
    if(restrict_networkB_to_A){
      appendStringToFile(pvalfname, "\nrestricted network B to network A\n")
    }else{appendStringToFile(pvalfname, "\nrestricted network B modules to network A\n")}
    write.table(rbind(colnames(oFetPMtrx)), pvalfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)
    write.table(oFetPMtrxR, pvalfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)

    write.table(rbind(colnames(oFoldMtrx)), freqfname, sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(oFoldMtrx, freqfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)
    if(restrict_networkB_to_A){
      appendStringToFile(pvalfname, "\nrestricted network B to network A\n")
    }else{appendStringToFile(freqfname, "\nrestricted network B modules to network A\n")}
    write.table(rbind(colnames(oFoldMtrx)), freqfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
    write.table(oFoldMtrxR, freqfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)

    #-------------------------------------------------------------------------------------
    #------------ find overlap genes for each group pair          ------------------------
    #
    #print("Merging two matrices")
    selA = aaMatrix[,ncol(aaMatrix)] != "grey"
    selB = bbMatrix[,ncol(bbMatrix)] != "grey"

    abMtrx = merge(aaMatrix[selA, ], bbMatrix[selB,c(1, ncol(bbMatrix))], by.x=1, by.y=1, all=FALSE)
    abMtrx = as.matrix(abMtrx); ncx = ncol(abMtrx)
    comparison = paste(abMtrx[, ncx-1], abMtrx[, ncx])
    abMtrx2= cbind(abMtrx, comparison)
    write.table(abMtrx2, outfname, sep="\t", quote=FALSE, col.names=TRUE, row.names=F)

    mergedMembers = tapply(abMtrx[,1], INDEX=comparison, concatenate, mysep=",", do_union=TRUE, do_sort=TRUE)
    selGS = !is.null(genesymbolIdx)
    if(selGS){selGS = genesymbolIdx!=uniqueIdCol[1] }

    if(selGS){
       mergedMembersGS= tapply(abMtrx[,genesymbolIdx+1], INDEX=comparison, concatenate, mysep=",", do_union=TRUE, do_sort=TRUE)
       mergedMembersMtrx = cbind(names(mergedMembers), mergedMembers, mergedMembersGS)
       colnames(mergedMembersMtrx) = c("Comparison", "SharedProbes", "SharedGenes")
    } else {
       mergedMembersMtrx = cbind(names(mergedMembers), mergedMembers)
       colnames(mergedMembersMtrx) = c("Comparison", "SharedMembers")
    }
    #mergedMembersModNames = getAllParts(names(mergedMembers), " ")
    #indexnamesA = c(1:length(sizes.modulesA)); names(indexnamesA) = names(sizes.modulesA)
    #indexnamesB = c(1:length(sizes.modulesB)); names(indexnamesB) = names(sizes.modulesB)
    #mergedModNamesIdxy = cbind(indexnamesA[mergedMembersModNames[,1]], indexnamesB[mergedMembersModNames[,2]])
    
    opairMtrx2 = merge(opairMtrx,  mergedMembersMtrx, by.x=1, by.y=1, all.x=TRUE, sort=FALSE)
    opairMtrx2R= merge(opairMtrxR, mergedMembersMtrx, by.x=1, by.y=1, all.x=TRUE, sort=FALSE)

    pidx = ncol(opairMtrx) - (ncol(mergedMembersMtrx)-1) -1
    tmpP = as.numeric(opairMtrx[, pidx]); od = order(tmpP)

    write.table(rbind(colnames(opairMtrx2)), pairfname, sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(opairMtrx2, pairfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)

    write.table(rbind(colnames(opairMtrx2R)), pairfnameR, sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(opairMtrx2R, pairfnameR, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)

}



# restrict_networkB_to_A ==T: take networkB as background when mapping modulesB to network A
# restrict_networkB_to_A ==F: consider modulesB as independent, so map individual modules to network A
#
moduleBasedIntersection=function(fnameA, fnameB,outputDir="", keywords="", uniqueIdCol=c(1,1), geneInforCols=8, itotalGenes = 3600, signifpvalue=0.001, latex=F, removeDuplicate=F, restrict_networkB_to_A=T)
{
    #fkeyB       =getFileName(fnameB)
    #outfname    paste(outputDir, "intersect.csv", sep='')

    # get unique labels for each input files
    keys= unlist(strsplit(keywords, "-vs-"))
    fkeyA       =keys[1]
    fkeyB       =keys[2]

    #outfname     = paste(outputDir,  "intersect_", keywords, ".xls", sep='')
    outfname     = paste(outputDir,  "intersect_", keywords, ".xls", sep='')

    nofname     = paste(outputDir,  "intersect_", keywords, "_number.xls", sep='')
    pvalfname   = paste(outputDir,  "intersect_", keywords, "_pvalue.xls", sep='')
    freqfname   = paste(outputDir,  "intersect_", keywords, "_frequence.xls", sep='')


    latexnofname  = paste(outputDir, "intersect_", keywords, "_number.tex", sep='')
    latexpvalfname= paste(outputDir, "intersect_", keywords, "_pvalue.tex", sep='')

    aaMatrixAll <- read.delim(fnameA, sep="\t", header=T)
    dim(aaMatrixAll)
    aaMatrixAll <- as.matrix(aaMatrixAll)

    bbMatrixAll <- read.delim(fnameB, sep="\t", header=T)
    dim(bbMatrixAll)
    bbMatrixAll <- as.matrix(bbMatrixAll)

    aaMatrixAll[, uniqueIdCol[1]] = toupper(aaMatrixAll[, uniqueIdCol[1]])
    bbMatrixAll[, uniqueIdCol[2]] = toupper(bbMatrixAll[, uniqueIdCol[2]])

    no.colsA = dim(aaMatrixAll)[2]
    no.colsB = dim(bbMatrixAll)[2]

    # include gene infor in the aaMatrix but not bbMatrix, in this way the merged matrix
    #  will not have duplicated columns
    aaMatrix <- aaMatrixAll[, c(uniqueIdCol[1], 1:geneInforCols, no.colsA)]
    colnames(aaMatrix) <- c( colnames(aaMatrixAll)[c(uniqueIdCol[1], 1:geneInforCols)], fkeyA)

    bbMatrix <- bbMatrixAll[, c(uniqueIdCol[2], no.colsB)]
    colnames(bbMatrix) <- c( colnames(bbMatrixAll)[1], fkeyB)


    #--------------  total overlap --------------------
    AA=union(as.character(as.matrix(aaMatrix[, 1])), NULL)
    BB=union(as.character(as.matrix(bbMatrix[, 1])), NULL)
    #totaloverlap    = merge(AA,BB, by.x=1, by.y=1, sort=F, all=F)
    #no.totalOverlap = dim(totaloverlap)[1]
    totaloverlap    = intersect(AA, BB)
    no.totalOverlap = length(totaloverlap)

    totalbackground = length(union(AA, BB))
    if (!is.na(itotalGenes)){
        totalbackground = itotalGenes
    }else{
        totalbackground = length(AA)
    }
    
    # global enrichment test
    #
    pval = phyper(no.totalOverlap-1, length(AA), itotalGenes-length(AA), length(BB), lower.tail=F)
    globalEnrich = c(no.totalOverlap , pval, length(AA), length(BB))
    
    # matched to aaMatrix (mnore conservative test)
    #
    ccMatrix <- merge(bbMatrix, cbind(totaloverlap), by.x=1, by.y=1, all=F)
    ccMatrix <- as.matrix(ccMatrix)

    #-------------- module-based overlap --------------------
    #
    imoduleColA=dim(aaMatrix)[2]
    imoduleColB=dim(bbMatrix)[2]

    modlevelsA= names(table(aaMatrix[,imoduleColA]) )
    modlevelsB= names(table(bbMatrix[,imoduleColB]) )
    
    #modlevelsA= levels(aaMatrix[,imoduleColA]) #levels(aaMatrix$module)
    #modlevelsA
    #modlevelsB= levels(bbMatrix[,imoduleColB]) #levels(bbMatrix$module)
    #modlevelsB

    no.modulesA = length(modlevelsA)
    no.modulesB = length(modlevelsB)

    latexchline2 = rep("\\\\\\hline", no.modulesA + 2)
    latexchline  = rep("\\\\\\hline", no.modulesA + 1)

    internoMatrix = matrix(NA, no.modulesA+1, no.modulesB)
    pvalMatrix    = matrix(NA, no.modulesA+1, no.modulesB)
    pvalMatrix_Rst= matrix(NA, no.modulesA,   no.modulesB) # restricted test (to AA)

    freqMatrix     = matrix(0, no.modulesA+1, no.modulesB)
    freqMatrix_Rst = matrix(0, no.modulesA+1,   no.modulesB)

    genenoInModuelsA =c()
    genenoInModuelsB =c()
    allmergedmatrix=NULL

    # overlap between AA and individual BB modules
    #
    #overlap_bbmod2AA = rep(0, no.modulesB+ 1); overlap_bbmod2AA[1]=no.totalOverlap;
    AAX = cbind(AA, rep("ALL", length(AA)) )
    for (j in c(1:no.modulesB) ) {
        comrslt = compareTwoModules(gifA=AAX, gifB=bbMatrix, moduleNameA="ALL", moduleNameB=modlevelsB[j], uniqueIdCol=1, 
                     moduleColA=2,moduleColB=2, totalGenes=totalbackground, removeDuplicate=removeDuplicate)
        interrslt =comrslt[[1]];
        internoMatrix[1,j] = interrslt[1]
        pvalMatrix [1,j]   = interrslt[2]
                
        if(restrict_networkB_to_A) {
           freqMatrix[1, j]    = length(BB)/totalbackground # frequency of this moduleB in the whole population
           freqMatrix_Rst[1, j]= no.totalOverlap/length(AA)  # frequency of this moduleB in the network A
        }else{
           freqMatrix_Rst[1, j]= interrslt[1]/length(AA)  # frequency of this moduleB in the network A
           freqMatrix[1, j]    = interrslt[4]/totalbackground # frequency of this moduleB in the whole population
        }
    }
    
    for (i in c(1:no.modulesA) ) {
      for (j in c(1:no.modulesB) ) {
        
        comrslt = compareTwoModules(aaMatrix, bbMatrix, modlevelsA[i], modlevelsB[j], uniqueIdCol=1, 
                        moduleColA=imoduleColA,moduleColB=imoduleColB, totalGenes=totalbackground, removeDuplicate=removeDuplicate)

        comrsltRst= compareTwoModules(gifA=aaMatrix, gifB=ccMatrix, moduleNameA=modlevelsA[i], 
                        moduleNameB=modlevelsB[j], uniqueIdCol=1, 
                        moduleColA=imoduleColA,moduleColB=imoduleColB, 
                        totalGenes=length(AA), removeDuplicate=removeDuplicate)

        interrslt =comrslt[[1]]; interrslt_Rst =comrsltRst[[1]]
        internoMatrix[i+1,j] = interrslt[1]
        pvalMatrix [i+1,j]   = interrslt[2]
        pvalMatrix_Rst[i,j]= interrslt_Rst[2]

        freqMatrix[i+1, j] =(interrslt[1]/interrslt[3])/freqMatrix[1, j];
        freqMatrix_Rst[i+1, j] =(interrslt_Rst[1]/interrslt_Rst[3])/freqMatrix_Rst[1, j];

        if (i== 1) {
           genenoInModuelsB =c(genenoInModuelsB, interrslt[4])
        }

        #if ( (interrslt[2] <= signifpvalue) & (modlevelsA[i]!="grey") & (modlevelsB[j]!="grey") ){
            allmergedmatrix=rbind(allmergedmatrix, comrslt[[2]])
        #}
      }
      genenoInModuelsA =c(genenoInModuelsA, interrslt[3])
    }

    #add in gene no in second column, names first column
    xinternoMatrix = cbind(c(length(AA),as.character(genenoInModuelsA)), internoMatrix)
    yinternoMatrix = cbind(c("ALL", as.character(modlevelsA)),       xinternoMatrix)

    #add in gene no in second row, names first row
    yinternoMatrix = rbind( c("total", as.character(no.totalOverlap), as.character(genenoInModuelsB)), yinternoMatrix)
    zinternoMatrix = rbind( c("overlap", "total", as.character(modlevelsB)), yinternoMatrix)

    write.table(zinternoMatrix, nofname,   sep="\t", quote=FALSE, col.names=F, row.names=F)

    if(latex) {
    zinternoMatrixLatex = zinternoMatrix
    for (i in c(1: (no.modulesA+2) ))
         zinternoMatrixLatex[i, no.modulesB+2] = paste(zinternoMatrixLatex[i, no.modulesB+2], "\\\\\\hline", sep="")
    write.table(zinternoMatrixLatex, latexnofname,   sep=" &", quote=FALSE, col.names=F, row.names=F)
    }

    #xpvalMatrix = cbind(as.character(genenoInModuelsA), pvalMatrix)
    #ypvalMatrix = cbind(as.character(modlevelsA),       xpvalMatrix)
    #ypvalMatrix = rbind( c("total", as.character(no.totalOverlap), as.character(genenoInModuelsB)), ypvalMatrix)
    #zpvalMatrix = rbind( c("overlap", "total", as.character(modlevelsB)), ypvalMatrix)

    # --------------------- pvalue --------------------------------
    #
    sigMatrix = pvalMatrix < signifpvalue
    sigMatrix_Rst = pvalMatrix_Rst < signifpvalue

    threshpMatrix = ifelse(sigMatrix, signif(pvalMatrix,3), "")
    threshpMatrix_Rst= ifelse(sigMatrix_Rst, signif(pvalMatrix_Rst,3), "")

    #zpvalMatrix = cbind(as.character(modlevelsA), signif(pvalMatrix,3))
    zpvalMatrix = cbind(c("ALL", as.character(modlevelsA)), threshpMatrix)
    zpvalMatrix = rbind( c("p.value", as.character(modlevelsB)), zpvalMatrix)
    write.table(zpvalMatrix,    pvalfname, sep="\t", quote=FALSE, col.names=F, row.names=F)

    if(restrict_networkB_to_A){
      appendStringToFile(pvalfname, "\nrestricted network B to network A\n")
    }else{appendStringToFile(pvalfname, "\nrestricted network B modules to network A\n")}

    zpvalMatrix = cbind(as.character(modlevelsA), threshpMatrix_Rst)
    zpvalMatrix = rbind( c("p.value", as.character(modlevelsB)), zpvalMatrix)
    write.table(zpvalMatrix, pvalfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)


    # --------------------- frequency/fold-change --------------------------------
    #
    
    zpvalMatrix = cbind(c("ALL_freq", as.character(modlevelsA)),  freqMatrix)
    zpvalMatrix = rbind( c("freq/fold_change", as.character(modlevelsB)), zpvalMatrix)
    write.table(zpvalMatrix,    freqfname, sep="\t", quote=FALSE, col.names=F, row.names=F)

    if(restrict_networkB_to_A){
      appendStringToFile(freqfname, "\nrestricted network B to network A\n")
    }else{appendStringToFile(freqfname, "\nrestricted network B modules to network A\n")}

    zpvalMatrix = cbind(c("ALL_freq", as.character(modlevelsA)), freqMatrix_Rst)
    zpvalMatrix = rbind( c("freq/fold_change", as.character(modlevelsB)), zpvalMatrix)
    write.table(zpvalMatrix, freqfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)


    if(latex) {
    zpvalMatrixLatex = zpvalMatrix
    for (i in c(1: (no.modulesA+1)) )
         zpvalMatrixLatex[i, no.modulesB+1] = paste(zpvalMatrixLatex[i, no.modulesB+1], "\\\\\\hline", sep="")
    write.table(zpvalMatrixLatex,    latexpvalfname, sep=" &", quote=FALSE, col.names=F, row.names=F)
    }

    if ( !is.null(allmergedmatrix)){
      colsmerged=dim(allmergedmatrix)[2]
      colnames(allmergedmatrix) <- c( colnames(allmergedmatrix)[c(1:(colsmerged-1))], keywords)
      write.table(allmergedmatrix, outfname, sep="\t", quote=FALSE, col.names=T, row.names=F)
      retname= outfname
    }else{
      retname=NULL
    }

    fn = paste("intersect_", keywords, sep="")
    #prefixes=fn; indir=outputDir; fname=fn
    merge_GO_PvalFreqNumb(prefixes=fn, indir=outputDir, fname=fn)

    retname
}


# convert multiple matrices
#
gomatrices_to_pairs = function (ymtrxNumb, ymtrxPval, ymtrxFreq) {

  otitleNumb <- colnames(ymtrxNumb)
  otitlePval <- colnames(ymtrxPval)
  otitleFreq <- colnames(ymtrxFreq)

  colnames(ymtrxNumb) <- paste("Number",      otitleNumb)
  colnames(ymtrxFreq) <- paste("Fold_Change", otitleFreq)
  colnames(ymtrxPval) <- paste("Pvalue",      otitlePval)

  ymtrxNumb[1:3,]
  ymtrxFreq[1:3,]
  ymtrxPval[1:3,]

  dim(ymtrxNumb)
  dim(ymtrxFreq)
  dim(ymtrxPval)

  # merge
  final = merge(ymtrxNumb, ymtrxFreq, by.x=1, by.y=1, all=T)
  final = merge(final,    ymtrxPval, by.x=1, by.y=1, all=T)
  final = as.matrix(final)
  colnames(final) <- c("Module", colnames(final)[-1])

  # matrix form
  #write.table(final, foutymtrx, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

  # pairwise form
  if(nrow(ymtrxPval)>2 & ncol(ymtrxPval)>2) {
    yymtrxNumb = as.matrix(ymtrxNumb)[-1,-c(1,2)]
    yymtrxFreq = as.matrix(ymtrxFreq)[-1,-1]
    yymtrxPval = as.matrix(ymtrxPval)[-1,-1]
  } else if(nrow(ymtrxPval)==2) {
    yymtrxNumb = rbind(as.matrix(ymtrxNumb)[-1,-c(1,2)])
    yymtrxFreq = rbind(as.matrix(ymtrxFreq)[-1,-1])
    yymtrxPval = rbind(as.matrix(ymtrxPval)[-1,-1])
  } else {
    yymtrxNumb = cbind(as.matrix(ymtrxNumb)[-1,-c(1,2)])
    yymtrxFreq = cbind(as.matrix(ymtrxFreq)[-1,-1])
    yymtrxPval = cbind(as.matrix(ymtrxPval)[-1,-1])
  }
  dim(yymtrxNumb)
  dim(yymtrxFreq)
  dim(yymtrxPval)

  modulesizes   = ymtrxNumb[-1,2]
  singaturesize = ymtrxNumb[1,-c(1,2)]

  if(is.null(dim(yymtrxNumb)) ) {
     index     = c(1: length(yymtrxNumb) )
  } else {
     index     = getMatrixIndex(dim(yymtrxNumb), symmetric=F, diagonal=T)
  }

  yrownames = as.character(as.matrix(ymtrxFreq)[-1,1])
  ycolnames = getSecondPart(colnames(ymtrxFreq)[-1], " ", 2)

  datMatrix = cbind(yymtrxNumb[index], yymtrxFreq[index], yymtrxPval[index])

  if(is.null(dim(yymtrxNumb)) ) {
    final2    = cbind(yrownames[index], modulesizes[index], 
                    rep(ycolnames, length(index)), rep(singaturesize,length(index)), datMatrix)     
  } else {
    final2    = cbind(yrownames[index[,1]], modulesizes[index[,1]], 
                    ycolnames[index[,2]], singaturesize[index[,2]], datMatrix)
  }
 
  sel = !is.na(datMatrix[,3])  
  final2 = final2[sel, ] 
  colnames(final2) <- c("module", "modulesize", "signature", "signature_size", "overlap","fold_change", "pvalue" )

  return (final2)
}

# merge three GO files
#
merge_GO_PvalFreqNumb = function(prefixes, indir="", fname="", specialHBTRC=F)
{

nocases = length(prefixes)

flag = "restricted network B modules to network A"

# output
#foutMtrx = paste(indir, fname, "_matrix.xls", sep="")
foutPair = paste(indir, fname, "_pair.xls",   sep="")
foutPair2= paste(indir, fname, "_pair-restricted.xls",   sep="")

xfinal = NULL
xfinal2= NULL

for (i in c(1:nocases) ) {

  fmemb = paste(indir, prefixes[i], ".xls",    sep="")
  fpval = paste(indir, prefixes[i], "_pvalue.xls",    sep="")
  ffreq = paste(indir, prefixes[i], "_frequence.xls", sep="")
  fnumb = paste(indir, prefixes[i], "_number.xls",    sep="")

  print( paste(i, prefixes[i]) )

  mtrxPval = getSecondTable(fpval, firstTable=T, endflag=flag, 
               blankline_before=1, blankline_after=1)

  mtrxFreq = getSecondTable(ffreq, firstTable=T, endflag=flag, 
               blankline_before=1, blankline_after=1)

  mtrxNumb <- read.delim(fnumb,sep="\t", header=T)
  mtrxNumbR<- as.matrix(mtrxNumb)[-1,] # restricted
  mtrxNumb <- as.matrix(mtrxNumb)[-2,]
  mtrxNumb[1,1]="ALL"
  mtrxNumbR[1,1]="ALL"

  mtrxFreq <- as.matrix(mtrxFreq)
  mtrxPval <- as.matrix(mtrxPval)

  mtrxMemb <- read.delim(fmemb,sep="\t", header=T)
  mtrxMemb <- as.matrix(mtrxMemb)
  mtrxMemb[1:3,]

  ncc= dim(mtrxMemb)[2]
  membBYmod = tapply(X=mtrxMemb[,1], INDEX=mtrxMemb[,ncc], FUN=concatenate, mysep=",")

  #[1,] "CR"  "A"  "N"  "all" 
  #[2,] "CR"  "A"  "N"  "down"
  #
  ncols = dim(mtrxNumb)[2]

  if(specialHBTRC) {
    parts = getAllParts( colnames(mtrxNumb)[-c(1:2)], "\\.")
    sel=parts[,4]!="all"; 
    selIdx1   = c(1,2, c(3:ncols)[sel])
    mtrxNumb  = mtrxNumb[, selIdx1]
    mtrxNumbR = mtrxNumbR[, selIdx1]
    ncols = dim(mtrxPval)[2]
    selIdx2   = c(1, c(2:ncols)[sel])
    mtrxPval  = mtrxPval[,selIdx2]
    mtrxFreq  = mtrxFreq[,selIdx2]
  }

  #ymtrxNumb=mtrxNumb; ymtrxPval=mtrxPval; ymtrxFreq=mtrxFreq 
  tmp=gomatrices_to_pairs(ymtrxNumb=mtrxNumb, ymtrxPval=mtrxPval, ymtrxFreq=mtrxFreq )
  rnTmp = paste(tmp[,1], tmp[,3], sep=".")
  tmp2  = cbind(tmp, membBYmod[rnTmp])
  colnames(tmp2) <- c(colnames(tmp), "members")
  xfinal = rbind(xfinal, tmp2)
  

  # 2nd table
  #
  mtrxPval2 = getSecondTable(fpval, firstTable=FALSE, endflag=flag, 
               blankline_before=1, blankline_after=1)

  mtrxFreq2 = getSecondTable(ffreq, firstTable=FALSE, endflag=flag, 
               blankline_before=1, blankline_after=1)
  mtrxFreq2 <- as.matrix(mtrxFreq2)
  mtrxPval2 <- as.matrix(mtrxPval2)

  if(specialHBTRC) {
     mtrxPval2  = mtrxPval2[,selIdx2]
     mtrxFreq2  = mtrxFreq2[,selIdx2]
  }

  nsize=dim(mtrxPval2)
  mtrxPval3  = rbind(rep(1,nsize[2]), mtrxPval2); colnames(mtrxPval3) <- colnames(mtrxPval2)
  tmp2  = gomatrices_to_pairs(ymtrxNumb=mtrxNumbR, ymtrxPval=mtrxPval3, ymtrxFreq=mtrxFreq2)

  rnTmp = paste(tmp2[,1], tmp2[,3], sep=".")
  tmp3  = cbind(tmp2, membBYmod[rnTmp])
  colnames(tmp3) <- c(colnames(tmp2), "members")

  xfinal2 = rbind(xfinal2, tmp3) 

}

ncols = dim(xfinal2)[2]-1
od = order(as.numeric(xfinal2[,ncols]))
xfinal2 = xfinal2[od, ] 

od = order(as.numeric(xfinal[,ncols]))
xfinal = xfinal[od, ] 

# matrix form
write.table(xfinal, foutPair, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)
write.table(xfinal2, foutPair2, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

}


# convert multiple matrices
#
gomatrices_to_pairs = function (ymtrxNumb, ymtrxPval, ymtrxFreq) {

  otitleNumb <- colnames(ymtrxNumb)
  otitlePval <- colnames(ymtrxPval)
  otitleFreq <- colnames(ymtrxFreq)

  colnames(ymtrxNumb) <- paste("Number",      otitleNumb)
  colnames(ymtrxFreq) <- paste("Fold_Change", otitleFreq)
  colnames(ymtrxPval) <- paste("Pvalue",      otitlePval)

  ymtrxNumb[1:3,]
  ymtrxFreq[1:3,]
  ymtrxPval[1:3,]

  dim(ymtrxNumb)
  dim(ymtrxFreq)
  dim(ymtrxPval)

  # merge
  final = merge(ymtrxNumb, ymtrxFreq, by.x=1, by.y=1, all=T)
  final = merge(final,    ymtrxPval, by.x=1, by.y=1, all=T)
  final = as.matrix(final)
  colnames(final) <- c("Module", colnames(final)[-1])

  # matrix form
  #write.table(final, foutymtrx, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

  # pairwise form
  if(nrow(ymtrxPval)>2 & ncol(ymtrxPval)>2) {
    yymtrxNumb = as.matrix(ymtrxNumb)[-1,-c(1,2)]
    yymtrxFreq = as.matrix(ymtrxFreq)[-1,-1]
    yymtrxPval = as.matrix(ymtrxPval)[-1,-1]
  } else if(nrow(ymtrxPval)==2) {
    yymtrxNumb = rbind(as.matrix(ymtrxNumb)[-1,-c(1,2)])
    yymtrxFreq = rbind(as.matrix(ymtrxFreq)[-1,-1])
    yymtrxPval = rbind(as.matrix(ymtrxPval)[-1,-1])
  } else {
    yymtrxNumb = cbind(as.matrix(ymtrxNumb)[-1,-c(1,2)])
    yymtrxFreq = cbind(as.matrix(ymtrxFreq)[-1,-1])
    yymtrxPval = cbind(as.matrix(ymtrxPval)[-1,-1])
  }
  dim(yymtrxNumb)
  dim(yymtrxFreq)
  dim(yymtrxPval)

  modulesizes   = ymtrxNumb[-1,2]
  singaturesize = ymtrxNumb[1,-c(1,2)]

  if(is.null(dim(yymtrxNumb)) ) {
     index     = c(1: length(yymtrxNumb) )
  } else {
     index     = getMatrixIndex(dim(yymtrxNumb), symmetric=F, diagonal=T)
  }

  yrownames = as.character(as.matrix(ymtrxFreq)[-1,1])
  ycolnames = getSecondPart(colnames(ymtrxFreq)[-1], " ", 2)

  datMatrix = cbind(yymtrxNumb[index], yymtrxFreq[index], yymtrxPval[index])

  if(is.null(dim(yymtrxNumb)) ) {
    final2    = cbind(yrownames[index], modulesizes[index], 
                    rep(ycolnames, length(index)), rep(singaturesize,length(index)), datMatrix)     
  } else {
    final2    = cbind(yrownames[index[,1]], modulesizes[index[,1]], 
                    ycolnames[index[,2]], singaturesize[index[,2]], datMatrix)
  }
 
  sel = !is.na(datMatrix[,3])  
  final2 = final2[sel, ] 
  colnames(final2) <- c("module", "modulesize", "signature", "signature_size", "overlap","fold_change", "pvalue" )

  return (final2)
}


# get all indices of a matrix
# for symetric matrix size is a number otherwise size=(row, col)
getMatrixIndex = function(size, symmetric=T, diagonal=F)
{
   allidx = NULL

   if(symmetric){
      for(i in c(1:(size-1)) ) {
         iv = cbind(i, (i+1):size)
         allidx = rbind(allidx, iv)
      }
      if(diagonal) {allidx = rbind(allidx, cbind(1:size, 1:size) )}

   } else {
      for(i in c(1:(size[1])) ) {
         iv = cbind(i, 1:(size[2]))
         allidx = rbind(allidx, iv)
      }
   }

   return (allidx)
}



# endflag="": end of the first table and start of the 2nd tabke
# blankline_before: a blank line before endflag?
# blankline_after=0: a blank line before endflag?
# skiplines: skip the first line

getSecondTable = function(inputfname, firstTable=T, endflag="", blankline_before=0, blankline_after=0, skiplines=0)
{
     lastrow = 1
     while (1==1) {
        iline <- read.delim(inputfname,sep="\t", header=F, skip=lastrow, nrow=1)
        iline <- as.character(as.matrix(iline))
        if (length(iline)==1) {
           if (iline==endflag) {break;}
        }
        if(!is.na(iline[1]) ) {
           if (iline[1]==endflag) {break;}
        }

        lastrow <- lastrow + 1    
     }

     allMatrix <- read.delim(inputfname,sep="\t", header=T, nrow=lastrow-blankline_before-skiplines, skip=skiplines)
     if (firstTable){
        return(allMatrix)
     }
     
     allMatrix <- read.delim(inputfname,sep="\t", header=T, skip=lastrow + 1 + blankline_after)
     return (allMatrix)
}

#get the filename without extension
#
getFileExtension=function(fullfname){
    splitted=unlist( strsplit(fullfname, "\\.") )
    
    if( length(splitted) >1){
      return (splitted[length(splitted)])
    } else{
      return ("")
    }
}

#get the filename without extension
getFileName=function(fullfname){
    ext=getFileExtension(fullfname)
    if(ext ==""){
       return (fullfname)
    }
    extd = paste(".", ext, sep="")
    splitted=splitString(fullfname, extd)

    splitted[1]
}

#get the filename without extension
getFileNames=function(fullfnames){

  final = NULL
  for(ef in fullfnames) {
     fn = getFileName(ef)
     final = c(final, fn)
  }
  return (final)
}


# get second part: 31357-31351 ==> 31351 
#
getSecondPart=function(fullfnames, sep="-", whichpart=-1, retBlank_ifNoMatch=F){

   n.elements = length(fullfnames)
   p1 = rep(T, n.elements);
   if(sep!="") {
     p1= strINstrings(sep, fullfnames)
   }
  
  ret = NULL
  for(each in fullfnames) {
    splitted=unlist( strsplit(each, sep) )
    if (whichpart==-1) {
       reti=splitted[ length(splitted) ]
    } else {
       if ( whichpart > length(splitted) ) {
          reti= each
       } else {
          reti=splitted[whichpart]
       }
    }

    ret = c(ret, reti)
  }

  if(retBlank_ifNoMatch) {
    ret[!p1] = ""
  }

  ret
}

compareTwoModules = function(gifA, gifB, moduleNameA, moduleNameB, uniqueIdCol=1, moduleColA, moduleColB, totalGenes, removeDuplicate=F)
{
  restrictA = gifA[, moduleColA]== moduleNameA 
  restrictB = gifB[, moduleColB]== moduleNameB
  
  restrictA = ifelse(is.na(gifA[,uniqueIdCol]), F, restrictA)
  restrictB = ifelse(is.na(gifB[,uniqueIdCol]), F, restrictB)

  noA= sum(restrictA)
  noB= sum(restrictB)
  moduleSetA = rbind(gifA[restrictA, ])
  moduleSetB = rbind(gifB[restrictB, ])

  if (noA==0 | noB==0){
     ret  = c(0, 1, noA, noB)
     # here we include merged matrix
     finalret= list(ret, NULL)
     return (finalret)
  }
  
  if (removeDuplicate) {
    orderA    = order(moduleSetA[,uniqueIdCol])
    moduleSetA=       rbind(moduleSetA[orderA,])
    boundA    = findBoundary(moduleSetA[,uniqueIdCol])
    moduleSetA= rbind(moduleSetA[boundA,])

    orderB = order(moduleSetB[,uniqueIdCol])
    moduleSetB=    rbind(moduleSetB[orderB,])
    boundB    = findBoundary(moduleSetB[,uniqueIdCol])
    moduleSetB=    rbind(moduleSetB[boundB,])

    noA= dim(moduleSetA)[1]
    noB= dim(moduleSetB)[1]
  }

  mergeMatrix = merge(moduleSetA, moduleSetB, by.x=uniqueIdCol, by.y=uniqueIdCol, sort=F,all=FALSE)
  intersectNo = dim(mergeMatrix)[1]

  # new module assignment
  combmodulename=paste(moduleNameA, ".", moduleNameB, sep="")
  newmergeMatrix = cbind(mergeMatrix, rep(combmodulename, intersectNo) )

  # phyper(89,702,4000-702,280, lower.tail=F)
  if(intersectNo>0) {
    pval = phyper(intersectNo-1, noA, totalGenes-noA, noB, lower.tail=F)

  }else{
    pval = 1
  }

  ret  = c(intersectNo, pval, noA, noB)

  # here we include merged matrix
  finalret= list(ret, newmergeMatrix)
  finalret
}


appendStringToFile=function(fname, mstring, newline=T){
    fp <- file(fname, "a")
    if(newline){
     cat(mstring, "\n", file=fp, sep="")
    }else{
     cat(mstring, file=fp)
    }
    close(fp)    
}

appendListToFile=function(fname, mlist, listTitle=""){
  fp <- file(fname, "a")
  if (length(listTitle)>0){
       cat(as.character(listTitle), "\n", file=fp)
  }  
      no.fields=length(mlist)
      #write column title
      for (z in 1:no.fields ){
         cat(as.character(names(mlist[z])),"\t", file=fp)
      }      
      cat("\n", file=fp)

      for (z in 1:no.fields ){
         itab=mlist[z]
         cat(as.character(itab[[1]]),"\t", file=fp)
      }      
      cat("\n", file=fp)
      close(fp)
 }


appendTableToFile=function(fname, mtable, tableTitle="", myappend=T){
   if ( is.null(mtable) ){
     return
   }

   if(myappend==T){
     fp <- file(fname, "a")    
   }else{
     fp <- file(fname, "w")
   }
   if (tableTitle != "" ){
     cat(as.character(tableTitle), "\n", file=fp)    
   }

  if ( (!is.na( dim(mtable)[1])) & is.na( dim(mtable)[2]) ) {#only one row in the table
      #write column title
      coltitles=names(mtable)
      for (z in 1:length(coltitles) ){
         cat(as.character(coltitles[z]),"\t", file=fp)
      }
      cat("\n", file=fp)

      for (i in 1:(dim(mtable)[1]) ){
          cat(as.character(mtable[i]), "\t", file=fp)
      }
      cat("\n", file=fp)

   }else{ # normal table
       cat(" \t", file=fp)
       #write column title
       coltitles=colnames(mtable)
       for (z in 1:length(coltitles) ){
          cat(as.character(coltitles[z]),"\t", file=fp)
       }
       cat("\n", file=fp)

       rowsname = rownames(mtable)
       for (i in 1:(dim(mtable)[1]) ){
           cat(as.character(rowsname[i]), "\t", file=fp)
           for(j in 1:(dim(mtable)[2])){
              cat(as.character(mtable[i, j]), "\t", file=fp)
            }
           cat("\n", file=fp)
       }   
   }
   cat("\n", file=fp)
   close(fp)
}

appendMultiTablesToFile=function(fname, multitables){
    #table of tables
    if(is.na( dim(multitables)[2]) ) {
        titles=names(multitables)
        for (i in 1:(dim(multitables)[1]) ){
            if ( is.null(multitables[[i]]) )
                 next
            appendTableToFile(fname, multitables[[i]], as.character(titles[i]))
         }
    }else{#single table
      appendTableToFile(fname, multitables)
    }
}

concatenate=function(myvect, mysep="", do_union=FALSE, do_sort=FALSE)
{
  noitems = length(myvect)
  if (noitems==0){
    return ("")
  }else if (noitems==1){
    return (as.character(myvect) )
  }
  
  if(do_union) {
     myvect2=unique(as.character(myvect))
  } else{
     myvect2 = myvect
  }
  if(do_sort){myvect2=sort(myvect2)}
  concatenated <- paste(myvect2, sep="", collapse=mysep)

  return (concatenated)

  #tmpfn = "tmp.txt"
  #write.table(t(as.character(myvect)),tmpfn,sep=mysep,quote=FALSE, col.names=F, row.names=FALSE)
  #concatenated <- read.delim(tmpfn, sep="!", header=F)
  #return (as.character(as.matrix(concatenated) ))
}



# to split "abc|123", use sep="\\|", "abc.123" use "\\."
splitString =function(mystring, separator="; "){
  splitted = NULL
  for (each in mystring){
     if (is.na(each) | is.null(each)){
        next
     }
     a=unlist( strsplit(each, separator) )
     splitted =c(splitted, a)
  }
  #a=unlist( strsplit(mystring, separator) )
  return(splitted )
}

getTopNChars = function(mystrs, N=4){
    res = NULL
    for(str in mystrs){
       chars = splitString(str, "")
       cstr = concatenate(chars[c(1:N)], "")
       res  = c(res, cstr)
    }
    return (res)
}

splitStringsAsLists =function(mystring, separator="; "){
  nelemts = length(mystring)
  splitted = as.list(rep(NA, nelemts))
  for (i in c(1:nelemts) ){
     each = mystring[i]
     if (is.na(each) | is.null(each)){
        next
     }
     a=unlist( strsplit(each, separator) )
     splitted[[i]]=a
  }
  #a=unlist( strsplit(mystring, separator) )
  return(splitted )
}

# "HOXA9K"                  "RTHB6 /// LOC650428"     "TAKEDA_NUP8_HOXA9_8D_DN"
#   ===>
#
# "HOXA9K" "RTHB6"     "TAKEDA_NUP8_HOXA9_8D_DN"
# "HOXA9K" "LOC650428" "TAKEDA_NUP8_HOXA9_8D_DN"
#
splitMiddleElementNames = function(vect, elementIdx)
{
    #print(vect)

    igenes = splitString(vect[elementIdx], " /// ")
    igenes = setdiff(igenes, c("T", "L", "I") )
    igenes = replaceString(igenes, " ///", "")
    igenes = replaceString(igenes, " //", "")
    if(length(igenes)==1){
        newvect = vect; newvect[elementIdx]=igenes
        return (newvect)
    }

    # matrix construction row based
    #
    newMatrix = t(matrix(rep(vect, length(igenes)), nrow=length(vect)))
    newMatrix[,elementIdx] = igenes
    
    return (newMatrix)
}

elementInSplitString = function(vect, sep=" ", element=c("breast"))
{
    #print(vect)
    igenes = replaceString(vect, ", ", " ")
    igenes = replaceString(igenes, ",", " ")
    igenes = toupper(splitString(igenes, sep))
    
    return (is.element(element, igenes))
}


getStringLength = function(mystring, separator="; "){
  nelemts = length(mystring)
  mylen = rep(0, nelemts)
  for (i in c(1:nelemts) ){
     each = mystring[i]
     if (is.na(each) | is.null(each)){
        next
     }
     a=unlist( strsplit(each, separator) )
     mylen[i]=length(a)
  }
  return(mylen)
}


strINstr = function(substr, istring){
   splitted = splitString(istring, substr)
   if(length(splitted)==1){return (F);
   }else{return(T);}
}

strINstrings = function(substr, strings){
   ns  = length(strings)
   res = rep(F, ns)
   for(i in c(1:ns) ){
      istring = strings[i]
      if(istring==""){next}
      splitted = splitString(istring, substr)
      if(length(splitted)>1){res[i]=T;}
   }
   return (res)
}

union_length =function(mvect)
{ 
  return (length(union(mvect,NULL)))
}


# a= 1 1 2 2 3 3 3
# findBoundary(a) => [1] 1 3 5
#
findBoundary=function(vect){
   no.elements = length(vect)
   esel = vect[no.elements]==vect[1]
   esel = ifelse(is.na(esel), F, esel)

   #if(esel){
   #  return (1);
   #}

   shifted= c(vect[no.elements], vect[1:(no.elements-1)] )
   sel = shifted != vect

   # handle NA
   sel = ifelse(is.na(sel), F, sel)

   sel[1] = TRUE

   return ( c(1:no.elements)[sel] )
}

