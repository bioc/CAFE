ProcessCels <- function(threshold.over=1.5, threshold.under=(2/3), remove_method=1, local_file=NULL){
    
    message("processing....\n This can take a while,you can grab a coffee in the meantime\n")
    
    t <- try(affyBatch<-ReadAffy())
    if("try-error" %in% class(t)){
        message("Oops! Reading cell files failed.Are you sure the working directory is correct?")
        #TODO: go back to beginning
    }
    
    eset <- rma(affyBatch)
    callsEset <- mas5calls(affyBatch)
    calls <- exprs(callsEset)
    rm(affyBatch)
    
    #calculate >20% absent exclusion vector
    message("Calculating absent probes\n")
    finalv <- data.frame(ID=rownames(calls),retain=FALSE,stringsAsFactors=FALSE)
    trufalsev <- 0
    trufalsev[1:length(calls[,1])] <- FALSE
    totl <- length(calls[1,])
    per20 <- totl * 0.2
    for(i in 1:length(calls[,1])){ # this is gonna take some time
        tmpcalls <- calls[i,]
        sel <- tmpcalls[tmpcalls %in% "A"]
        selL <- length(sel)
        if(selL > (as.integer(per20)+1)){
        }
        else{
            trufalsev[i] <- TRUE
        }
    }
    
    finalv$retain <- as.logical(trufalsev)
    
    write.exprs(eset,"Expression_values.xls") #is this necessary?
    
    #get chromosomal lomessageions
    anno <- eset@annotation
    annodb <- paste(anno,".db",sep="")
    biocLite <- NULL #prevent R CMD check from thinking this
    # variable doesnt exist
    f <- try(library(annodb,character.only=TRUE))
    if(class(f)=="try-error"){
        if (!requireNamespace("BiocManager", quietly=TRUE))
            install.packages("BiocManager")
        BiocManager::install(annodb)
        library(annodb,character.only=TRUE)
    }
    
    ex <- exprs(eset)
    #de-log (raises 2 to the power n) of all values in ex.
    ex <- 2^ex
    av <- rowMedians(ex)
    dex <- data.frame(ex)
    
    ov <- data.frame(a=dex[1])
    message('Calculating relative expressions\n')
    for(i in 1:length(dex)){
        subs <- dex[i]
        rel <- subs/av
        nam <- names(dex)[i]
        ov[nam] <- rel  	
    }
    IDs <- featureNames(eset)
    Symbols <- getSYMBOL(IDs,anno)
    ChrNum <- as.character(sapply(lookUp(IDs,annodb,"CHR"),"[",i=1))
    Bands <- as.character(sapply(lookUp(IDs,annodb,"MAP"),"[",i=1))
    findArms <- function(Bands){
        Arms <- 1:length(Bands)
        
        which_p <- grep("p",Bands)
        which_q <- grep("q",Bands)
        chr_p <- sapply(strsplit(grep("p",Bands,value=TRUE),"p"),"[",1)
        chr_q <- sapply(strsplit(grep("q",Bands,value=TRUE),"q"),"[",1)
        arm_p <- paste(chr_p,"p",sep="")
        arm_q <- paste(chr_q,"q",sep="")
        
        Arms[which_p] <- arm_p
        Arms[which_q] <- arm_q
        
        return(Arms)
    }
    Arms <- findArms(Bands)
    message("Removing absent probes\n")
    IDs <- IDs[finalv$retain]
    Symbols <- Symbols[finalv$retain]
    ChrNum <- ChrNum[finalv$retain]
    Bands <- Bands[finalv$retain]
    Arms <- Arms[finalv$retain]
    dex <- dex[finalv$retain,]
    ov <- ov[finalv$retain,]
    if(annodb=="hgu133plus2.db" || annodb=="hgu133a.db" || annodb=="hgu95a.db"){
        #new method to find locs
        tmpdf <- data.frame(ID=IDs,Sym=Symbols,Chr=ChrNum,Loc=NA,Band=Bands,
                            Arm=Arms,stringsAsFactors=FALSE)
        tmpdf <- cbind(tmpdf,dex)
        tmpdf <- cbind(tmpdf,ov)
        affy_table <- NULL
        if(length(local_file)==0){
            if(annodb=="hgu133plus2.db"){
                affy_table <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/affyU133Plus2.txt.gz"
            }
            else if(annodb=="hgu133a.db"){
                affy_table <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/affyU133.txt.gz"
            }
            else if(annodb=="hgu95a.db"){
                affy_table <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/affyU95.txt.gz"
            }
            tmpfile <- tempfile(tmpdir=getwd())
            file.create(tmpfile)
            download.file(affy_table,tmpfile)
            big_table <- read.table(tmpfile)
        }
        
        else if(length(local_file)!=0){
            big_table <- read.table(local_file)
        }
        
        if(annodb=="hgu133a.db"){
            big_table[,11] <- sapply(strsplit(as.character(big_table[,11]),":"),
                                     "[[",2)
            big_table[,11] <- sapply(strsplit(as.character(big_table[,11]),";"),
                                     "[[",1)
        }
        big_table <- big_table[!duplicated(big_table[,11]),]
        big_table <- big_table[order(as.character(big_table[,11])),]
        tmpdf <- tmpdf[which(IDs %in% as.character(big_table[,11])),]
        tmpdf$Loc <- big_table[,17][which(
            as.character(big_table[,11]) %in% IDs)]
        IDs <- tmpdf$ID
        Symbols <- tmpdf$Sym
        ChrNum <- tmpdf$Chr
        ChrLoc <- tmpdf$Loc
        Bands <- tmpdf$Band
        Arms <- tmpdf$Arm
        dex <- tmpdf[,7:(6+length(dex))]
        ov <- tmpdf[,(length(dex)+7):length(tmpdf)]
    }
    
    else{ #other annotation packages
        ChrLoc <- abs(as.integer(sapply(lookUp(IDs,annodb,"CHRLOC"),"[",i=1)))
    } 
    
    message("Removing duplicatees\n")
    fina <- 1:length(IDs)
    fina[1:length(IDs)] <- FALSE
    tmp <- data.frame(ID=IDs,Sym=Symbols,stringsAsFactors=FALSE)
    uniques <- unique(tmp$Sym)
    if(remove_method == 1){ #normal
        pb <- txtProgressBar(min=1,max=length(uniques),style=3)
        for(i in 1:length(uniques)){
            setTxtProgressBar(pb,i)
            sel <- tmp$Sym == uniques[i]
            sel <- which(sel)
            if(length(sel)==1){
                fina[sel] <- TRUE
            }
            if(length(sel)>1){
                ab <- grep("[[:digit:]_]_at",tmp$ID[sel])
                if(length(ab)!=0){ #has nnn_at
                    selsel <- sel[ab]
                    fina[selsel] <- TRUE
                    if(length(selsel)>1){
                        probes <- tmp$ID[selsel]
                        probes_n <- unlist(strsplit(probes,"_at"))
                        min_probe <- min(as.numeric(probes_n))
                        probes_n <- probes_n[!probes_n==min_probe]
                        probes_n <- paste(probes_n,"_at",sep="")
                        delete <- which(tmp$ID %in% probes_n)
                        fina[delete] <- FALSE
                    }
                    deselect <- 1:length(sel)
                    deselect[1:length(sel)] <- TRUE
                    deselect[ab] <- FALSE
                    deselect <- as.logical(deselect)
                    nonsel <- sel[deselect]
                    fina[nonsel] <- FALSE
                }
                
                else{
                    ab <- grep("[[:digit:]_]a_at",tmp$ID[sel])
                    if(length(ab)!=0){ #has nnn_a_at but not nnn_at
                        selsel <- sel[ab]
                        fina[selsel] <- TRUE
                        if(length(selsel)>1){
                            if(length(selsel)>1){
                                probes <- tmp$ID[selsel]
                                probes_n <- unlist(strsplit(probes,"_a_at"))
                                min_probe <- min(as.numeric(probes_n))
                                probes_n <- probes_n[!probes_n==min_probe]
                                probes_n <- paste(probes_n,"_a_at",sep="")
                                delete <- which(tmp$ID %in% probes_n)
                                fina[delete] <- FALSE
                            }
                        }
                        deselect <- 1:length(sel)
                        deselect[1:length(sel)] <- TRUE
                        deselect[ab] <- FALSE
                        deselect <- as.logical(deselect)
                        nonsel <- sel[deselect]          
                        fina[nonsel] <- FALSE
                    }
                    
                    else{
                        ab <- grep("[[:digit:]_]s_at",tmp$ID[sel])
                        if(length(ab)!=0){ #has nnn_s_at but not nnn_at
                            selsel <- sel[ab]
                            fina[selsel] <- TRUE
                            if(length(selsel)>1){
                                if(length(selsel)>1){
                                    probes <- tmp$ID[selsel]
                                    probes_n <- unlist(strsplit(probes,"_s_at"))
                                    min_probe <- min(as.numeric(probes_n))
                                    probes_n <- probes_n[!probes_n==min_probe]
                                    probes_n <- paste(probes_n,"_s_at",sep="")
                                    delete <- which(tmp$ID %in% probes_n)
                                    fina[delete] <- FALSE
                                }
                            }
                            deselect <- 1:length(sel)
                            deselect[1:length(sel)] <- TRUE
                            deselect[ab] <- FALSE
                            deselect <- as.logical(deselect)
                            nonsel <- sel[deselect]          
                            fina[nonsel] <- FALSE
                        }
                        
                        else{ #has nnn_x_at or nnn_a_at but not nnn_s_at or
                            # nnnn_at
                            ab <- grep("[[:digit:]_]x_at",tmp$ID[sel])
                            selsel <- sel[ab]
                            fina[selsel] <- TRUE
                            if(length(selsel)>1){
                                if(length(selsel)>1){
                                    probes <- tmp$ID[selsel]
                                    probes_n <- unlist(strsplit(probes,"_x_at"))
                                    min_probe <- min(as.numeric(probes_n))
                                    probes_n <- probes_n[!probes_n==min_probe]
                                    probes_n <- paste(probes_n,"_x_at",sep="")
                                    delete <- which(tmp$ID %in% probes_n)
                                    fina[delete] <- FALSE
                                }
                            }
                            deselect <- 1:length(sel)
                            deselect[1:length(sel)] <- TRUE
                            deselect[ab] <- FALSE
                            deselect <- as.logical(deselect)
                            nonsel <- sel[deselect]
                            fina[nonsel] <- FALSE
                        }
                    }
                }
            }
        }
    }
    
    else if(remove_method == 2){ #only remove probes if multiple map to same
        # location
        pb <- txtProgressBar(min=1,max=length(uniques),style=3)
        tmp <- data.frame(ID=IDs,Sym=Symbols,Loc=ChrLoc,stringsAsFactors=FALSE)
        for(i in 1:length(uniques)){      
            setTxtProgressBar(pb,i)
            tmpresult <- tmp[which(tmp$Sym %in% uniques[i]),]
            if(length(unique(tmpresult$Loc))==length(tmpresult$Loc)){
                fina[which(tmp$Sym %in% tmpresult$Sym)] <- TRUE
            }
            else{
                dups <- tmpresult[which(
                    tmpresult$Loc %in% tmpresult$Loc[duplicated(
                        tmpresult$Loc)]),]
                tmpv <- 1:length(tmpresult[,1])
                tmpv[1:length(tmpresult[,1])] <- TRUE
                tmpv[which(tmpresult$Loc %in% dups$Loc)] <- FALSE
                tmpv <- as.logical(tmpv)
                nondups <- tmpresult[tmpv,]
                fina[which(tmp$ID %in% nondups$ID)] <- TRUE
                dupID <- dups$ID
                ab <- grep("[[:digit:]_]_at",dupID)
                if(length(ab)!=0){
                    ats <- dupID[ab]
                    if(length(ats)==1){
                        fina[which(IDs %in% ats)] <- TRUE
                    }
                    else if(length(ats)>1){
                        probes_n <- unlist(strsplit(ats,"_at"))
                        min_probe <- min(probes_n)
                        el <- grep(min_probe,ats)
                        fina[which(IDs %in% ats[el])] <- TRUE
                    }
                    
                }
                
                ab <- grep("[[:digit:]_]a_at",dupID)
                if(length(ab)!=0){
                    aats <- dupID[ab]
                    if(length(aats)==1){
                        fina[which(IDs %in% aats)] <- TRUE
                    }
                    else if(length(aats)>1){
                        probes_n <- unlist(strsplit(aats,"_a_at"))
                        min_probe <- min(probes_n)
                        el <- grep(min_probe,aats)
                        fina[which(IDs %in% aats[el])] <- TRUE
                    }
                    
                }
                
                ab <- grep("[[:digit:]_]s_at",dupID)
                if(length(ab)!=0){
                    sats <- dupID[ab]
                    if(length(sats)==1){
                        fina[which(IDs %in% sats)] <- TRUE
                    }
                    else if(length(sats)>1){
                        probes_n <- unlist(strsplit(sats,"_s_at"))
                        min_probe <- min(probes_n)
                        el <- grep(min_probe,sats)
                        fina[which(IDs %in% sats[el])] <- TRUE
                    }
                    
                }
                
                ab <- grep("[[:digit:]_]x_at",dupID)
                if(length(ab)!=0){
                    xats <- dupID[ab]
                    if(length(xats)==1){
                        fina[which(IDs %in% xats)] <- TRUE
                    }
                    
                    else if(length(xats)>1){
                        probes_n <- unlist(strsplit(xats,"_x_at"))
                        min_probe <- min(probes_n)
                        el <- grep(min_probe,xats)
                        fina[which(IDs %in% xats[el])] <- TRUE
                    }                
                }        
            }
        }    
    }
    
    else if(remove_method == 0){ #dont remove at all
        fina[1:length(fina)] <- TRUE
    } 
    
    fina <- as.logical(fina)
    IDs <- IDs[fina]
    Symbols <- Symbols[fina]
    ChrLoc <- ChrLoc[fina]
    ChrNum <- ChrNum[fina]
    Bands <- Bands[fina]
    Arms <- Arms[fina]
    dex <- dex[fina,]
    ov <- ov[fina,]
    
    
    
    nSamples <- length(sampleNames(eset))
    exprsList <- list()
    for(o in 1:nSamples){
        expressions <- data.frame(ID=IDs,Sym=Symbols,Value=dex[,o],Rel=ov[,o],
                                  Loc=ChrLoc,Chr=ChrNum,Band=Bands,Arm=Arms,
                                  stringsAsFactors=FALSE)
        expressions <- na.omit(expressions) # REMOVES NA chromlocs
        names(expressions) <- c("ID","Sym","Value","LogRel","Loc","Chr",
                                "Band","Arm")
        exprsList[[o]] <- expressions
    }
    names(exprsList) <- sampleNames(eset)
    
    overList <- list()
    underList <- list()
    message("\n")
    message('Finding over- and underexpressed probes\n')
    for(v in 1:length(exprsList)){
        nam <- names(exprsList[v])
        tmp <- exprsList[[v]]
        oversel <- tmp[,4] > threshold.over
        undersel <- tmp[,4] < threshold.under
        tmp2 <- tmp[oversel,]
        tmp3 <- tmp[undersel,]
        overList[[nam]] <- tmp2
        underList[[nam]] <- tmp3
    }
    
    message('Log2 transforming\n')
    for(x in 1:length(exprsList)){
        exprsList[[x]][,3] <- log2(exprsList[[x]][,3])
        exprsList[[x]][,4] <- log2(exprsList[[x]][,4])
        overList[[x]][,4] <- log2(overList[[x]][,4])
        overList[[x]][,3] <- log2(overList[[x]][,3])
        underList[[x]][,3] <- log2(underList[[x]][,3])
        underList[[x]][,4] <- log2(underList[[x]][,4])
    }
    
    combine <- list(whole=exprsList,over=overList,under=underList)
    return(combine)
    
}
