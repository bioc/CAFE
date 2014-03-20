fisher.method <- function(pvals){
    p <- pchisq(-2*sum(log(pvals)),df=(2*length(pvals)),lower.tail=FALSE)
    return(p)  
}

makeContingency <- function(chromNum,overExpressedSample,wholeSample){
    wholeSample <- wholeSample[!(wholeSample$ID %in% overExpressedSample$ID),]
    overExist <- 0
    overFalse <- 0
    overExist <- overExist + 
        sum(overExpressedSample$Chr==as.character(chromNum))
    overFalse <- overFalse + 
        sum(overExpressedSample$Chr!=as.character(chromNum))
    
    wholeExist <- 0
    wholeFalse <- 0
    wholeExist <- wholeExist + sum(wholeSample$Chr==as.character(chromNum))
    wholeFalse <- wholeFalse + sum(wholeSample$Chr!=as.character(chromNum))
    
    d <- data.frame(fal=c(wholeFalse,overFalse),tru=c(wholeExist,overExist))
    rownames(d) <- c("Whole","Over")
    return(d)  
}

fisherExact <- function(datalist,chromNum,samples=c(1,2),method="cli",
                        enrichment){
    result <- NULL
    subset <- NULL
    if(length(samples)==0){
        if(method == "cli"){
            subset <- cliSubset(datalist,enrichment)
        }
        
        if(method == "gui"){
            guiSubset(datalist,enrichment)
            message("Hit <ENTER> when finished selecting")
            nonsense <- readline()
            subset <- get("guiSelectedSet",envir=as.environment(globalenv()))
        }
    }
    
    else{
        subs_o <- NULL
        if(enrichment=="greater"){
            subs_o <- datalist$over[samples]
        }
        else if(enrichment=="less"){
            subs_o <- datalist$under[samples]
        }      
        subset <- list(samples=subs_o,whole=datalist$whole[samples])
    }
    
    if(chromNum == "ALL"){
        #COMPUTE FOR ALL CHROMOSOMES (EXCEPT SEX)
        for(i in 1:suppressWarnings(max(as.numeric(subset$whole[[1]]$Chr),
                                        na.rm=TRUE))){
            comb <- NULL
            for(x in 1:length(samples)){
                ctable <- makeContingency(i,subset$samples[[x]],
                                          subset$whole[[x]])
                test <- fisher.test(ctable,alternative="greater")
                test <- test$p.value
                comb <- append(comb,test)
            }
            pval <- fisher.method(comb)
            result <- append(result,pval)
            names(result)[i] <- paste("Chr",i,sep="")
        }
    }
    
    else{
        #COMPUTE FOR JUST 1 CHROMOSOME      
        comb <- NULL
        for(x in 1:length(samples)){
            ctable <- makeContingency(chromNum,subset$samples[[x]],
                                      subset$whole[[x]])
            test <- fisher.test(ctable,alternative="greater")
            test <- test$p.value
            comb <- append(comb,test)
        }
        result <- fisher.method(comb)
        names(result) <- paste('Chr',as.character(chromNum),sep="")
    }
    
    return(result)
}

chisqr <- function(datalist,chromNum,samples=c(1,2),method="cli",
                   enrichment){
    result<-NULL
    subset <- NULL
    if(length(samples)==0){
        if(method == "cli"){
            subset <- cliSubset(datalist,enrichment)
        }
        
        if(method == "gui"){
            guiSubset(datalist,enrichment)
            message("Hit <ENTER> when finished selecting")
            nonsense <- readline()
            subset <- get("guiSelectedSet",envir=as.environment(globalenv()))
        }
    }
    
    else{
        subs_o <- NULL
        if(enrichment=="greater"){
            subs_o <- datalist$over[samples]
        }
        else if(enrichment=="less"){
            subs_o <- datalist$under[samples]
        }  
        subset <- list(samples=subs_o,whole=datalist$whole[samples])
    }
    
    if(chromNum == "ALL"){
        #COMPUTE FOR ALL CHROMOSOMES (EXCEPT SEX)
        comb <- NULL
        for(i in 1:suppressWarnings(max(as.numeric(subset$whole[[1]]$Chr),
                                        na.rm=TRUE))){
            for(x in 1:length(samples)){
                ctable <- makeContingency(i,subset$samples[[1]],
                                          subset$whole[[1]])
                test <- chisq.test(ctable)
                pval <- test$p.value
                comb <- append(comb,pval)
            }
            result <- fisher.method(comb)
            pval <- result$p.value
            result <- append(result,pval) 
            names(result)[i] <- paste("Chr",i,sep="")
        }     
    }
    
    else{
        #COMPUTE FOR JUST 1 CHROMOSOME
        comb <- NULL
        for(x in 1:length(samples)){
            ctable <- makeContingency(chromNum,subset$samples[[1]],
                                      subset$whole[[1]])
            test <- chisq.test(ctable)
            result <- test$p.value
            comb <- append(comb,result)        
        }
        result <- fisher.method(comb)
        names(result) <- paste('Chr',as.character(chromNum),sep="")
    }
    
    return(result)
}

makeBandContingency <- function(chromNum,band,overExpressedSample,wholeSample){
    # TODO Make this faster here. Runs for a full 1.5 seconds on each iteration
    # make a contingency table based on a particular bands on a particular 
    # chromosome. 
    wholeSample <- wholeSample[!(wholeSample$ID %in% overExpressedSample$ID),]
    overExist <- 0
    overFalse <- 0
    overSel <- overExpressedSample[overExpressedSample$Chr==as.character(chromNum),]
    overExist <- overExist + 
        sum(overSel$Band==as.character(band))
    overFalse <- overFalse + 
        sum(overSel$Band!=as.character(band))
    
    wholeExist <- 0
    wholeFalse <- 0
    wholeSel <- wholeSample[wholeSample$Chr==as.character(chromNum),]
    wholeExist <- wholeExist + sum(wholeSel$Band==as.character(band))
    wholeFalse <- wholeFalse + sum(wholeSel$Band!=as.character(band))
    
    d <- data.frame(fal=c(wholeFalse,overFalse),tru=c(wholeExist,overExist))
    rownames(d) <- c("Whole","Over")
    
    return(d)  
    
}

allBandsOnChrom <- function(datalist,chromNum,samples=c(1,2),method,enrichment){
    #make NAs appear as 1.0
    subset <- NULL
    if(length(samples)==0){
        if(method == "cli"){
            subset <- cliSubset(datalist,enrichment)
        }
        
        if(method == "gui"){
            guiSubset(datalist,enrichment)
            message("Hit <ENTER> when finished selecting")
            nonsense <- readline()
            subset <- get("guiSelectedSet",envir=as.environment(globalenv()))
        }
    }
    
    else{
        subs_o <- NULL
        if(enrichment=="greater"){
            subs_o <- datalist$over[samples]
        }
        else if(enrichment=="less"){
            subs_o <- datalist$under[samples]
        }      
        subset <- list(samples=subs_o,whole=datalist$whole[samples])
    }
    ultimo <- NULL
    if(chromNum == 'ALL'){
        ultimo <- data.frame(Chr=NULL,Band=NULL,P=NULL)
        
        for(o in 1:suppressWarnings(max(as.numeric(subset$whole[[1]]$Chr),
                                        na.rm=TRUE))){
            result <- NULL
            sel <- subset$whole[[1]]$Chr==as.character(o)
            selection <- subset$whole[[1]][sel,]
            totbands <- as.factor(selection$Band)
            bands <- levels(totbands)
            for(i in 1:length(bands)){
                comb <- NULL
                for(x in 1:length(samples)){
                    ctable <- makeBandContingency(o,bands[i],
                                                  subset$samples[[x]],
                                                  subset$whole[[x]])
                    p <- fisher.test(ctable,alternative="greater")$p.value
                    comb <- append(comb,p)
                }
                pval <- fisher.method(comb)
                result <- append(result,pval)
                names(result)[i] <- as.character(bands[i])
            }
            tmpdf <- data.frame(Chr=o,Band=bands,P=result)
            ultimo <- rbind(ultimo,tmpdf)
        }
        p_s <- ultimo[,3]
        ultimo$P <- p_s
    }
    
    else{
        result <- NULL
        sel <- subset$whole[[1]][order(subset$whole[[1]]$ID),]
        sell <- sel$Chr==as.character(chromNum)
        selection <- sel[sell,]
        totbands <- as.factor(selection$Band)
        bands <- levels(totbands)
        for(i in 1:length(bands)){
            comb <- NULL
            for(x in 1:length(samples)){
                ctable <- makeBandContingency(chromNum,bands[i],
                                              subset$samples[[x]],
                                              subset$whole[[x]])
                p <- fisher.test(ctable,alternative="greater")$p.value
                comb <- append(comb,p)
            }
            pval <- fisher.method(comb)               
            ultimo <- append(ultimo,pval)
            names(ultimo)[i] <- as.character(bands[i])
        }
    }
    return(ultimo)
}

allBandsOnChromChi <- function(datalist,chromNum,samples=c(1,2),
                               method,enrichment){
    subset <- NULL
    if(length(samples)==0){
        if(method == "cli"){
            subset <- cliSubset(datalist,enrichment)
        }
        
        if(method == "gui"){
            guiSubset(datalist,enrichment)
            message("Hit <ENTER> when finished selecting")
            nonsense <- readline()
            subset <- get("guiSelectedSet",envir=as.environment(globalenv()))
        }
    }
    
    else{
        subs_o <- NULL
        if(enrichment=="greater"){
            subs_o <- datalist$over[samples]
        }
        else if(enrichment=="less"){
            subs_o <- datalist$under[samples]
        }      
        subset <- list(samples=subs_o,whole=datalist$whole[samples])
    }
    ultimo <- NULL
    if(chromNum == 'ALL'){
        ultimo <- data.frame(Chr=NULL,Band=NULL,P=NULL)
        for(o in 1:suppressWarnings(max(as.numeric(subset$whole[[1]]$Chr),
                                        na.rm=TRUE))){
            result <- NULL
            sel <- subset$whole[[1]]$Chr==as.character(o)
            selection <- subset$whole[[1]][sel,]
            totbands <- as.factor(selection$Band)
            bands <- levels(totbands)
            for(i in 1:length(bands)){
                comb <- NULL
                for(x in 1:length(samples)){
                    ctable <- makeBandContingency(o,bands[i],
                                                  subset$samples[[x]],
                                                  subset$whole[[x]])
                    p <- chisq.test(ctable)$p.value
                    comb <- append(comb,p)
                }
                pval <- fisher.method(comb)
                result <- append(result,pval)
                names(result)[i] <- as.character(bands[i])
            }
            tmpdf <- data.frame(Chr=o,Band=bands,P=result)
            ultimo <- rbind(ultimo,tmpdf)
        }
        p_s <- ultimo[,3]
        ultimo$P <- p_s
    }
    
    else{
        sel <- subset$whole[[1]]$Chr==as.character(chromNum)
        selection <- subset$whole[[1]][sel,]
        totbands <- as.factor(selection$Band)
        bands <- levels(totbands)
        for(i in 1:length(bands)){
            comb <- NULL
            for(x in 1:length(samples)){
                ctable <- makeBandContingency(chromNum,bands[i],
                                              subset$samples[[x]],
                                              subset$whole[[x]])
                p <- chisq.test(ctable)$p.value
                comb <- append(comb,p)
            }
            pval <- fisher.method(comb)
            ultimo <- append(ultimo,pval)
            names(ultimo)[i] <- as.character(bands[i])
        }
    }  
    return(ultimo)
}

makeArmContingency <- function(chromNum,arm,overExpressedSample,wholeSample){
    arm <- paste(chromNum,arm,sep="")
    wholeSample <- wholeSample[!(wholeSample$ID %in% overExpressedSample$ID),]
    overExist <- 0
    overFalse <- 0
    overExist <- overExist + length(grep(arm,overExpressedSample$Arm))
    overFalse <- overFalse + length(grep(arm,overExpressedSample$Arm,
                                         invert=TRUE))
    
    wholeExist <- 0
    wholeFalse <- 0
    wholeExist <- wholeExist + length(grep(arm,wholeSample$Arm))
    wholeFalse <- wholeFalse + length(grep(arm,wholeSample$Arm,invert=TRUE))
    
    d <- data.frame(fal=c(wholeFalse,overFalse),tru=c(wholeExist,overExist))
    rownames(d) <- c("Whole","Over")
    return(d)  
}

ArmfisherExact <- function(datalist,chromNum,arm,samples=c(1,2),method="cli",
                           enrichment){
    result <- NULL
    subset <- NULL
    if(length(samples)==0){
        if(method == "cli"){
            subset <- cliSubset(datalist,enrichment)
        }
        
        if(method == "gui"){
            guiSubset(datalist,enrichment)
            message("Hit <ENTER> when finished selecting")
            nonsense <- readline()
            subset <- get("guiSelectedSet",envir=as.environment(globalenv()))
        }
    }
    
    else{
        subs_o <- NULL
        if(enrichment=="greater"){
            subs_o <- datalist$over[samples]
        }
        else if(enrichment=="less"){
            subs_o <- datalist$under[samples]
        }      
        subset <- list(samples=subs_o,whole=datalist$whole[samples])
    }
    
    if(chromNum == "ALL"){
        #COMPUTE FOR ALL CHROMOSOMES (EXCEPT SEX)
        for(i in 1:suppressWarnings(max(as.numeric(subset$whole[[1]]$Chr),
                                        na.rm=TRUE))){
            comb <- NULL
            for(x in 1:length(samples)){
                ctable <- makeArmContingency(i,arm=arm,subset$samples[[x]],
                                             subset$whole[[x]])
                p <- fisher.test(ctable,alternative="greater")$p.value
                comb <- append(comb,p)
            }
            pval <- fisher.method(comb)
            result <- append(result,pval)
            names(result)[i] <- paste("Chr",i,arm,sep="")
        }    
    }
    
    else{
        #COMPUTE FOR JUST 1 CHROMOSOME
        comb <- NULL
        for(x in 1:length(samples)){
            ctable <- makeArmContingency(chromNum,arm=arm,subset$samples[[x]],
                                         subset$whole[[x]])
            p <- fisher.test(ctable,alternative="greater")$p.value
            comb <- append(comb,p)
        }
        result <- fisher.method(comb)
        names(result) <- paste('Chr',as.character(chromNum),arm,sep="")
    }  
    return(result)
}

Armchisqr <- function(datalist,chromNum,arm,samples=c(1,2),method="cli",
                      enrichment){
    result<-NULL
    subset <- NULL
    if(length(samples)==0){
        if(method == "cli"){
            subset <- cliSubset(datalist,enrichment)
        }
        
        if(method == "gui"){
            guiSubset(datalist,enrichment)
            message("Hit <ENTER> when finished selecting")
            nonsense <- readline()
            subset <- get("guiSelectedSet",envir=as.environment(globalenv()))
        }
    }
    
    else{
        subs_o <- NULL
        if(enrichment=="greater"){
            subs_o <- datalist$over[samples]
        }
        else if(enrichment=="less"){
            subs_o <- datalist$under[samples]
        }  
        subset <- list(samples=subs_o,whole=datalist$whole[samples])
    }
    
    if(chromNum == "ALL"){
        #COMPUTE FOR ALL CHROMOSOMES (EXCEPT SEX)
        for(i in 1:suppressWarnings(max(as.numeric(subset$whole[[1]]$Chr),
                                        na.rm=TRUE))){
            comb <- NULL
            for(x in 1:length(samples)){
                ctable <- makeArmContingency(i,arm=arm,subset$samples[[x]],
                                             subset$whole[[x]])
                p <- chisq.test(ctable)
                comb <- append(comb,p)
            }
            pval <- fisher.method(comb)
            result <- append(result,pval) 
            names(result)[i] <- paste("Chr",i,arm,sep="")
        }    
    }
    
    else{
        #COMPUTE FOR JUST 1 CHROMOSOME
        comb <- NULL
        for(x in 1:length(samples)){
            ctable <- makeArmContingency(chromNum,arm=arm,subset$samples[[x]],
                                         subset$whole[[x]])
            p <- chisq.test(ctable)
            comb <- append(comb,p)
        }
        result <- fisher.method(comb)
        names(result) <- paste('Chr',as.character(chromNum),arm,sep="")
    }  
    return(result)
}

chromosomeStats <- function(datalist, chromNum=1, samples=NULL, select="cli",
                            test="fisher", bonferroni = TRUE,
                            enrichment = "greater"){
    
    result <- NULL
    if(test=="fisher"){
        if(length(samples)!=0){
            result <- fisherExact(datalist,chromNum=chromNum,samples=samples,
                                  enrichment=enrichment)
        }
        else if(length(samples)==0){
            result <- fisherExact(datalist,chromNum=chromNum,samples=NULL,
                                  method=select,enrichment=enrichment)
        }
    }
    
    else if(test=="chisqr"){
        if(length(samples)!=0){
            result <-  chisqr(datalist,chromNum=chromNum,samples=samples,
                              enrichment=enrichment)
        }
        
        else if(length(samples)==0){
            result <- chisqr(datalist,chromNum=chromNum,samples=NULL,
                             method=select,enrichment=enrichment)
        }
    }
    
    if(bonferroni==TRUE){
        result <- p.adjust(result,method="bonferroni")
    }
    return(result)
    
}

bandStats <- function(datalist, chromNum=1, samples=NULL, select="cli",
                      test="fisher", bonferroni = TRUE, enrichment = "greater"){
    result <- NULL
    if(test=="fisher"){
        if(length(samples)!=0){
            result <- allBandsOnChrom(datalist,chromNum,samples,method=select,
                                      enrichment=enrichment)
        }
        else if(length(samples)==0){
            result <- allBandsOnChrom(datalist,chromNum,samples=NULL,
                                      method=select,enrichment=enrichment)
        }
        
    }
    
    else if(test=="chisqr"){
        if(length(samples)!=0){
            result <- allBandsOnChromChi(datalist,chromNum,samples,
                                         method=select,enrichment=enrichment)
        }
        
        else if(length(samples)==0){
            result <- allBandsOnChromChi(datalist,chromNum,samples=NULL,
                                         method=select,enrichment=enrichment)
        }
        
    }
    
    if(bonferroni==TRUE){
        if(class(chromNum)=="numeric"){
            result <- p.adjust(result,method="bonferroni")
        }
        else if(class(chromNum)=="character"){
            adj <- p.adjust(result[,3],method="bonferroni")
            result[,3] <- adj
        }      
    }
    return(result)
    
}

armStats <- function(datalist, chromNum=1, arm="q", samples=NULL, select="cli",
                     test="fisher", bonferroni = TRUE, enrichment = "greater"){
    
    result <- NULL
    if(test=="fisher"){
        if(length(samples)!=0){
            result <- ArmfisherExact(datalist,chromNum=chromNum,arm=arm,
                                     samples=samples,enrichment=enrichment)
        }
        else if(length(samples)==0){
            result <- ArmfisherExact(datalist,chromNum=chromNum,arm=arm,
                                     samples=NULL,method=select,
                                     enrichment=enrichment)
        }
    }
    
    else if(test=="chisqr"){
        if(length(samples)!=0){
            result <-  Armchisqr(datalist,chromNum=chromNum,arm=arm,
                                 samples=samples,enrichment=enrichment)
        }
        
        else if(length(samples)==0){
            result <- Armchisqr(datalist,chromNum=chromNum,arm=arm,
                                samples=NULL,method=select,
                                enrichment=enrichment)
        }
    }
    
    if(bonferroni==TRUE){
        result <- p.adjust(result,method="bonferroni")
    }
    return(result)
}
