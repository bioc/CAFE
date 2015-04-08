slidSmooth <- function(x,k){
    k_begin <- k+1
    k_end <- length(x) - k
    x_smooth <- NULL
    x_smooth[1:length(x)] <- NA
    for(i in k_begin:k_end){
        tmp <- sum(x[(i):(i+k)])/((k)+1)
        x_smooth[i] <- tmp 
    }
    
    for(i in 1:k){
        s <- max(c(1,(i-k)))
        t <- min(c(length(x),(i+k)))
        tmp <- sum(x[s:t])/((t-s+1))
        x_smooth[i] <- tmp
    }
    
    for(i in (length(x)-k+1):length(x)){
        s <- max(c(1,(i-k)))
        t <- min(c(length(x),(i+k)))
        tmp <- sum(x[s:t])/((t-s+1))
        x_smooth[i] <- tmp
    }
    return(x_smooth)
}

slidWithGaps <- function(sliding,len){
    c <- 0
    
    for(w in 1:length(sliding[,1])){
        tmploc1 <- sliding$Loc[w]
        tmploc2 <- sliding$Loc[(w+1)]
        if(!is.na(tmploc2[1])){
            if((tmploc2-tmploc1) > 1000000){
                c <- c+1
            }
        }
    }
    
    slidFinLoc <- NULL
    slidFinLoc[1:(len+c)] <- NA
    slidFinRel <- NULL
    slidFinRel[1:(len+c)] <- NA
    slidFinLab <- NULL
    slidFinLab[1:(len+c)] <- NA          
    slidingFin <- data.frame(Loc=slidFinLoc,rel=slidFinRel,samples=slidFinLab)
    
    for(q in 1:length(sliding[,1])){
        tmploc1 <- sliding$Loc[q]
        tmploc2 <- sliding$Loc[(q+1)]
        tmpdf2 <- sliding[(q+1),]
        
        if(!is.na(tmpdf2[1])){
            
            if((tmploc2 - tmploc1) > 1000000){
                tmpdf2 <- data.frame(Loc=mean(tmploc1:tmploc2),rel=NA,
                                     lab=sliding$samples[q],
                                     stringsAsFactors=FALSE)
            }
        }        
        slidingFin[q,] <- tmpdf2        
    }
    
    return(slidingFin)
}

slidPlot <- function(datalist,samples=c(1,2),chromNum=1,combine=FALSE,k=1,
                     idiogram=FALSE,file="default"){
    
    Loc <- NULL
    val <- NULL
    hg19IdeogramCyto <- NULL
    rel <- NULL
    lab <- NULL
    nrel <- NULL
    if(length(samples)==1){
        if(combine==FALSE){      
            wholeSet <- datalist$whole
            
            if(chromNum == 'ALL'){
                returnlist <- list()
                length(returnlist) <- suppressWarnings(max(
                    as.numeric(wholeSet[[1]]$Chr),na.rm=TRUE))
                for(o in 1:suppressWarnings(max(as.numeric(wholeSet[[1]]$Chr),
                                                na.rm=TRUE))){
                    Locs <- NULL
                    rels <- NULL
                    daf <- wholeSet[[samples]]
                    daf  <- daf[order(daf$Chr),]
                    sel <- daf$Chr == as.character(o)
                    daf <- daf[sel,]
                    daf <- daf[order(daf$Loc),]
                    Locs <- daf$Loc
                    rels <- daf[,4]
                    nrels <- slidSmooth(rels,k)
                    
                    nam <- NULL
                    if(file=="default"){
                        nam <- paste("Chr",as.character(o),"slid",
                                     as.character(samples),sep="")
                    }
                    
                    else{
                        nam <- file
                    }
                    png(nam,1600,1200)
                    message(paste("Writing plot to ",nam,sep=""))
                    message("\n")
                    dfr <- data.frame(Loc=Locs,rel=rels,nrel=nrels,
                                      stringsAsFactors=FALSE)
                    dfr$Loc <- (dfr$Loc/1000000)
                    
                    if(idiogram == FALSE){        
                        pl <- qplot(Loc,nrel,data=dfr,geom=c("point"))
                        if(max(dfr$nrel)> 1.5 | min(dfr$nrel)< (-1.5)){
                            pl <- pl + coord_cartesian(
                                ylim = c(min(dfr$nrel), max(dfr$nrel))) +
                                xlab(paste("Chromosome ",o," (Mb)",sep="")) +
                                ylab("Log2 relative expression")
                        }
                        else{
                            pl <- pl + coord_cartesian(ylim = c(-0.75,0.75)) +
                                xlab(paste("Chromosome ",o," (Mb)",sep="")) +
                                ylab("Log2 relative expression")
                        }            
                        print(pl)
                        dev.off()
                        returnlist[[o]] <- pl
                    }
                    
                    if(idiogram == TRUE){
                        data(hg19IdeogramCyto, envir=environment())
                        chr <- paste("chr",as.character(o),sep="")
                        suppressMessages(p <- plotIdeogram(hg19IdeogramCyto,
                                                              chr,xlabel=FALSE))
                        p <- p + theme(strip.background=element_blank(),
                                       strip.text=element_blank()) +
                            theme(rect=element_blank())
                        p2 <- ggplot(data=dfr,aes(Loc,nrel)) +
                            geom_point(aes(alpha=0.5)) +
                            theme(legend.position='bottom')
                        if(max(dfr$nrel)>1.5 | min(dfr$nrel)< (-1.5)){
                            p2 <- p2 + coord_cartesian(
                                ylim=c(min(dfr$nrel),max(dfr$nrel))) +
                                xlab(paste('Chromosome',o,sep='')) +
                                ylab("Log2 relative expression")
                        }
                        else{
                            p2 <- p2 + coord_cartesian(ylim = c(-0.75,0.75)) +
                                xlab(paste('Chromosome',o,sep='')) +
                                ylab("Log2 relative expression")
                        }            
                        gA <- ggplot_gtable(ggplot_build(p@ggplot))
                        gB <- ggplot_gtable(ggplot_build(p2))
                        maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3])
                        gA$widths[2:3] <- as.list(maxWidth)
                        gB$widths[2:3] <- as.list(maxWidth)
                        arranged <- grid.arrange(gA, gB, ncol=1,heights=c(1,5))
                        print(arranged)
                        dev.off()
                        returnlist[[o]] <- p2
                    }          
                }
                return(returnlist)
            }
            
            else{
                Locs <- NULL
                rels <- NULL
                daf <- wholeSet[[samples]]
                daf  <- daf[order(daf$Chr),]
                sel <- daf$Chr == as.character(chromNum)
                daf <- daf[sel,]
                daf <- daf[order(daf$Loc),]
                Locs <- daf$Loc
                rels <- daf[,4]
                nrels <- slidSmooth(rels,k)
                
                nam <- NULL
                if(file=="default"){
                    nam <- paste("Chr",as.character(chromNum),"slid",
                                 as.character(samples),sep="")
                }
                
                else{
                    nam <- file
                }
                dfr <- data.frame(Loc=Locs,rel=rels,nrel=nrels,
                                  stringsAsFactors=FALSE)
                dfr$Loc <- (dfr$Loc/1000000)
                png(nam,1600,1200)
                message(paste("Writing plot to ",nam,sep=""))
                message("\n")
                if(idiogram==FALSE){      
                    pl <- qplot(Loc,nrel,data=dfr,geom=c("point"))
                    if(max(dfr$nrel)> 1.5 | min(dfr$nrel)< (-1.5)){
                        pl <- pl + coord_cartesian(
                            ylim = c(min(dfr$nrel), max(dfr$nrel))) +
                            xlab(paste("Chromosome ",chromNum," (Mb)",sep="")) +
                            ylab("Log2 relative expression")
                    }
                    else{
                        pl <- pl + coord_cartesian(ylim = c(-0.75,0.75)) +
                            xlab(paste("Chromosome ",chromNum," (Mb)",sep="")) +
                            ylab("Log2 relative expression")
                    }          
                    print(pl)
                    dev.off()
                    return(pl)
                }
                
                if(idiogram==TRUE){
                    data(hg19IdeogramCyto, envir=environment())
                    chr <- paste("chr",as.character(chromNum),sep="")
                    suppressMessages(p <- plotIdeogram(hg19IdeogramCyto,
                                                          chr,xlabel=FALSE))
                    p <- p + theme(strip.background=element_blank(),
                                   strip.text=element_blank()) +
                        theme(rect=element_blank())
                    p2 <- ggplot(data=dfr,aes(Loc,nrel)) +
                        geom_point(aes(alpha=0.5)) +
                        theme(legend.position='bottom')
                    if(max(dfr$nrel)> 1.5 | min(dfr$nrel)< (-1.5)){
                        p2 <- p2 + coord_cartesian(
                            ylim=c(min(dfr$nrel),max(dfr$nrel))) +
                            xlab(paste('Chromosome',chromNum," (Mb)",sep='')) +
                            ylab("Log2 relative expression")
                    }
                    else{
                        p2 <- p2 + coord_cartesian(ylim = c(-0.75,0.75)) +
                            xlab(paste('Chromosome',chromNum," (Mb)",sep='')) +
                            ylab("Log2 relative expression")
                    }          
                    gA <- ggplot_gtable(ggplot_build(p@ggplot))
                    gB <- ggplot_gtable(ggplot_build(p2))
                    maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3])
                    gA$widths[2:3] <- as.list(maxWidth)
                    gB$widths[2:3] <- as.list(maxWidth)
                    arranged <- grid.arrange(gA, gB, ncol=1,heights=c(1,5))
                    print(arranged)
                    dev.off()
                    return(p2)
                }
            } 
        }
        
        if(combine==TRUE){
            wholeSet <- datalist$whole
            
            if(chromNum == 'ALL'){
                returnlist <- list()
                length(returnlist) <- suppressWarnings(
                    max(as.numeric(wholeSet[[1]]$Chr),na.rm=TRUE))
                for(o in 1:suppressWarnings(max(as.numeric(wholeSet[[1]]$Chr),
                                                na.rm=TRUE))){
                    Locs <- NULL
                    rels <- NULL
                    daf <- wholeSet[[samples]]
                    daf  <- daf[order(daf$Chr),]
                    sel <- daf$Chr == as.character(o)
                    daf <- daf[sel,]
                    daf <- daf[order(daf$Loc),]
                    Locs <- daf$Loc
                    rels <- daf[,4]
                    
                    nrels <- slidSmooth(rels,k)
                    
                    len <- length(rels)
                    labls <- NULL
                    labls2 <- NULL
                    labls[1:len] <- 'RAW'
                    labls2[1:len] <- 'SLIDING'
                    raw <- data.frame(Loc=Locs,rel=rels,samples=labls)
                    sliding <- data.frame(Loc=Locs,rel=nrels,samples=labls2)
                    slidingFin <- slidWithGaps(sliding,len)
                    
                    sam <- names(wholeSet)[[samples]]
                    nam <- NULL
                    if(file=="default"){
                        nam <- paste("Chr",as.character(o),"slid",
                                     as.character(samples),sep="")
                    }
                    
                    else{
                        nam <- file
                    }
                    png(nam,1600,1200)
                    message(paste("Writing plot to ",nam,sep=""))
                    message("\n")
                    raw$Loc <- (raw$Loc/1000000)
                    slidingFin$Loc <- (slidingFin$Loc/1000000)
                    if(idiogram==FALSE){
                        pl <- ggplot(raw,aes(Loc,rel))
                        pl <- pl + geom_point()
                        pl <- pl + geom_line(data=slidingFin,colour='green')
                        if(max(slidingFin$rel,na.rm=TRUE)>1.5 | min(
                            slidingFin$rel,na.rm=TRUE)< (-1.5)){
                            pl <- pl + coord_cartesian(ylim = c(
                                min(slidingFin$rel,na.rm=TRUE),
                                max(slidingFin$rel,na.rm=TRUE)))
                        }
                        else{
                            pl <- pl + coord_cartesian(ylim = c(-0.75,0.75))
                        }
                        pl <- pl + xlab(paste('Chromosome ',as.character(o),
                                              " (Mb)",sep="")) +
                            ylab('Relative expression to the median')
                        print(pl)
                        dev.off()
                        returnlist[[o]] <- pl
                    }
                    
                    if(idiogram==TRUE){
                        data(hg19IdeogramCyto, envir=environment())
                        chr <- paste("chr",as.character(o),sep="")
                        suppressMessages(p <- plotIdeogram(hg19IdeogramCyto,
                                                              chr,xlabel=FALSE))
                        p <- p + theme(strip.background=element_blank(),
                                       strip.text=element_blank()) +
                            theme(rect=element_blank())
                        p2 <- ggplot(data=raw,aes(Loc,rel)) +
                            geom_point(aes(color=samples),alpha=0.5) +
                            geom_line(data=slidingFin,size=1) +
                            theme(legend.position='bottom')
                        if(max(slidingFin$rel,na.rm=TRUE)> 1.5 | min(
                            slidingFin$rel,na.rm=TRUE)< (-1.5)){
                            p2 <- p2 + coord_cartesian(
                                ylim=c(min(slidingFin$rel,na.rm=TRUE),
                                       max(slidingFin$rel,na.rm=TRUE))) +
                                xlab(paste('Chromosome',o," (Mb)",sep='')) +
                                ylab("Log2 relative expression")
                        }
                        else{
                            p2 <- p2 + coord_cartesian(ylim = c(-0.75,0.75)) +
                                xlab(paste('Chromosome',o," (Mb)",sep='')) +
                                ylab("Log2 relative expression")
                        }            
                        gA <- ggplot_gtable(ggplot_build(p@ggplot))
                        gB <- ggplot_gtable(ggplot_build(p2))
                        maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3])
                        gA$widths[2:3] <- as.list(maxWidth)
                        gB$widths[2:3] <- as.list(maxWidth)
                        arranged <- grid.arrange(gA, gB, ncol=1,heights=c(1,5))
                        print(arranged)
                        dev.off()
                        returnlist[[o]] <- p2
                    }
                }
                return(returnlist)
            }
            
            else{
                Locs <- NULL
                rels <- NULL
                daf <- wholeSet[[samples]]
                daf  <- daf[order(daf$Chr),]
                sel <- daf$Chr == as.character(chromNum)
                daf <- daf[sel,]
                daf <- daf[order(daf$Loc),]
                Locs <- daf$Loc
                rels <- daf[,4]
                
                nrels <- slidSmooth(rels,k)
                len <- length(rels)
                labls <- NULL
                labls2 <- NULL
                labls[1:len] <- 'RAW'
                labls2[1:len] <- 'SLIDING'
                raw <- data.frame(Loc=Locs,rel=rels,samples=labls)
                sliding <- data.frame(Loc=Locs,rel=nrels,samples=labls2)
                slidingFin <- slidWithGaps(sliding,len)
                
                sam <- names(wholeSet)[[samples]]
                nam <- NULL
                if(file=="default"){
                    nam <- paste("Chr",as.character(chromNum),"slid",
                                 as.character(samples),sep="")
                }
                
                else{
                    nam <- file
                }
                png(nam,1600,1200)
                message(paste("Writing plot to ",nam,sep=""))
                message("\n")
                raw$Loc <- (raw$Loc/1000000)
                slidingFin$Loc <- (slidingFin$Loc/1000000)
                if(idiogram==FALSE){
                    pl <- ggplot(raw,aes(Loc,rel))
                    pl <- pl + geom_point()
                    pl <- pl + geom_line(data=slidingFin,colour='green')
                    if(max(slidingFin$rel,na.rm=TRUE)> 1.5 | min(
                        slidingFin$rel,na.rm=TRUE)< (-1.5)){
                        pl <- pl + coord_cartesian(
                            ylim = c(min(slidingFin$rel,na.rm=TRUE),
                                     max(slidingFin$rel,na.rm=TRUE)))
                    }
                    else{
                        pl <- pl + coord_cartesian(ylim = c(-0.75,0.75))
                    }          
                    pl <- pl + labs(title="Sliding plot") +
                        xlab(paste("Chromosome ",chromNum," (Mb)",sep="")) +
                        ylab("Log2 relative expression")
                    print(pl)
                    dev.off()
                    return(pl)
                }
                
                if(idiogram==TRUE){
                    data(hg19IdeogramCyto, envir=environment())
                    chr <- paste("chr",as.character(chromNum),sep="")
                    suppressMessages(p <- plotIdeogram(hg19IdeogramCyto,
                                                          chr,xlabel=FALSE))
                    p <- p + theme(strip.background=element_blank(),
                                   strip.text=element_blank()) +
                        theme(rect=element_blank())
                    p2 <- ggplot(data=raw,aes(Loc,rel)) +
                        geom_point(aes(color=samples),alpha=0.5) +
                        geom_line(data=slidingFin,size=1) +
                        theme(legend.position='bottom')
                    if(max(slidingFin$rel,na.rm=TRUE)> 1.5 | min(
                        slidingFin$rel,na.rm=TRUE)< (-1.5)){
                        p2 <- p2 + coord_cartesian(
                            ylim=c(min(slidingFin$rel,na.rm=TRUE),
                                   max(slidingFin$rel,na.rm=TRUE))) +
                            xlab(paste('Chromosome',chromNum," (Mb)",sep='')) +
                            ylab("Log2 relative expression")
                    }
                    else{
                        p2 <- p2 + coord_cartesian(ylim = c(-0.75,0.75)) +
                            xlab(paste('Chromosome',chromNum," (Mb)",sep='')) +
                            ylab("Log2 relative expression")
                    }          
                    gA <- ggplot_gtable(ggplot_build(p@ggplot))
                    gB <- ggplot_gtable(ggplot_build(p2))
                    maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3])
                    gA$widths[2:3] <- as.list(maxWidth)
                    gB$widths[2:3] <- as.list(maxWidth)
                    arranged <- grid.arrange(gA, gB, ncol=1,heights=c(1,5))
                    print(arranged)
                    dev.off()
                    return(p2)
                }        
            }
        }
    }
    
    if(length(samples)>1){
        if(combine==FALSE){
            wholeSet <- datalist$whole
            selSet <- wholeSet[as.numeric(samples)]
            
            if(chromNum == "ALL"){
                returnlist <- list()
                length(returnlist) <- suppressWarnings(
                    max(as.numeric(wholeSet[[1]]$Chr),na.rm=TRUE))
                for(g in 1:suppressWarnings(
                    max(as.numeric(wholeSet[[1]]$Chr),na.rm=TRUE))){
                    Locs <- NULL
                    rels <- NULL
                    labls <- NULL
                    nrels <- NULL
                    labls2 <- NULL
                    
                    for(v in 1:length(selSet)){
                        daf <- selSet[[v]]
                        daf  <- daf[order(daf$Chr),]
                        sel <- daf$Chr == as.character(g)
                        daf <- daf[sel,]
                        daf <- daf[order(daf$Loc),]
                        Locs <- append(Locs,daf$Loc)
                        rels <- append(rels,daf[,4])
                        
                        nrels <- slidSmooth(rels,k)
                        len <- length(daf$Chr)
                        tmp <- NULL
                        tmp[1:len] <- names(selSet)[v]
                        labls <- append(labls,tmp)
                        labls2 <- append(labls2,tmp)
                        
                    }
                    
                    raw <- data.frame(Loc=Locs,rel=rels,samples=labls,
                                      stringsAsFactors=FALSE)
                    sliding <- data.frame(Loc=Locs,rel=nrels,samples=labls2,
                                          stringsAsFactors=FALSE)
                    nam <- NULL
                    if(file=="default"){
                        nam <- paste("Chr",as.character(g),"slidMul",sep="")
                    }
                    
                    else{
                        nam <- file
                    }
                    png(nam,1600,1200)
                    message(paste("Writing plot to ",nam,sep=""))
                    message("\n")
                    dfr <- data.frame(Loc=Locs,rel=rels,nrel=nrels,
                                      samples=labls,stringsAsFactors=FALSE)
                    dfr$Loc <- (dfr$Loc/1000000)
                    if(idiogram == FALSE){
                        pl <- qplot(Loc,nrel,data=dfr,colour=samples,
                                    geom=c("point"))
                        if(max(dfr$nrel)> 1.5 | min(dfr$nrel)< (-1.5)){
                            pl <- pl + coord_cartesian(
                                ylim=c(min(dfr$nrel),max(dfr$nrel))) +
                                xlab(paste('Chromosome',g," (Mb)",sep='')) +
                                ylab("Log2 relative expression")
                        }
                        else{
                            pl <- pl + coord_cartesian(ylim = c(-0.75,0.75)) +
                                xlab(paste('Chromosome',g," (Mb)",sep='')) +
                                ylab("Log2 relative expression")
                        }
                        pl <- pl + guides(fill=guide_legend(title="Samples"))
                        print(pl)
                        dev.off()
                        returnlist[[g]] <- pl
                    }
                    
                    if(idiogram == TRUE){
                        data(hg19IdeogramCyto, envir=environment())
                        chr <- paste("chr",as.character(g),sep="")
                        suppressMessages(p <- plotIdeogram(hg19IdeogramCyto,
                                                              chr,xlabel=FALSE))
                        p <- p + theme(strip.background=element_blank(),
                                       strip.text=element_blank()) +
                            theme(rect=element_blank())
                        p2 <- ggplot(data=dfr,aes(Loc,nrel)) +
                            geom_point(aes(color=samples),alpha=I(1/2)) +
                            theme(legend.position='bottom')
                        if(max(dfr$nrel)>1.5 | min(dfr$nrel)< (-1.5)){
                            p2 <- p2 + coord_cartesian(
                                ylim=c(min(dfr$nrel),max(dfr$nrel))) +
                                xlab(paste('Chromosome',g," (Mb)",sep='')) +
                                ylab("Log2 relative expression")
                        }
                        else{
                            p2 <- p2 + coord_cartesian(ylim = c(-0.75,0.75)) +
                                xlab(paste('Chromosome',g," (Mb)",sep='')) +
                                ylab("Log2 relative expression")
                        }
                        p2 <- p2 + guides(fill=guide_legend(title="Samples"))
                        gA <- ggplot_gtable(ggplot_build(p@ggplot))
                        gB <- ggplot_gtable(ggplot_build(p2))
                        maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3])
                        gA$widths[2:3] <- as.list(maxWidth)
                        gB$widths[2:3] <- as.list(maxWidth)
                        arranged <- grid.arrange(gA, gB, ncol=1,heights=c(1,5))
                        print(arranged)
                        dev.off()
                        returnlist[[g]] <- p2
                    }
                }
                return(returnlist)
            }
            
            else{
                Locs <- NULL
                rels <- NULL
                nrels <- NULL
                labls <- NULL
                labls2 <- NULL
                for(v in 1:length(selSet)){
                    daf <- selSet[[v]]
                    daf  <- daf[order(daf$Chr),]
                    sel <- daf$Chr == as.character(chromNum)
                    daf <- daf[sel,]
                    daf <- daf[order(daf$Loc),]
                    Locs <- append(Locs,daf$Loc)
                    rels <- append(rels,daf[,4])
                    
                    nrels <- slidSmooth(rels,k)
                    tmp <- NULL
                    len <- length(daf$Chr)
                    tmp[1:len] <- names(selSet)[v]
                    labls <- append(labls,tmp)
                    labls2 <- append(labls2,tmp)
                }
                
                nam <- NULL
                if(file=="default"){
                    nam <- paste("Chr",as.character(chromNum),"slidMul",sep="")
                }
                
                else{
                    nam <- file
                }
                png(nam,1600,1200)
                message(paste("Writing plot to ",nam,sep=""))
                message("\n")
                dfr <- data.frame(Loc=Locs,rel=rels,nrel=nrels,samples=labls,
                                  stringsAsFactors=FALSE)
                dfr$Loc <- (dfr$Loc/1000000)
                if(idiogram == FALSE){
                    pl <- qplot(Loc,nrel,data=dfr,colour=samples,
                                geom=c("point"))
                    if(max(dfr$nrel)>1.5 | min(dfr$nrel)< (-1.5)){
                        pl <- pl + coord_cartesian(
                            ylim = c(min(dfr$nrel), max(dfr$nrel))) +
                            xlab(paste("Chromosome ",chromNum," (Mb)",sep="")) +
                            ylab("Log2 relative expression") 
                    }
                    else{
                        pl <- pl + coord_cartesian(ylim = c(-0.75,0.75)) +
                            xlab(paste("Chromosome ",chromNum," (Mb)",sep="")) +
                            ylab("Log2 relative expression")
                    }
                    pl <- pl + guides(fill=guide_legend(title="Samples"))
                    print(pl)
                    dev.off()
                    return(pl)
                }
                
                if(idiogram == TRUE){
                    data(hg19IdeogramCyto, envir=environment())
                    chr <- paste("chr",as.character(chromNum),sep="")
                    suppressMessages(p <- plotIdeogram(hg19IdeogramCyto,
                                                          chr,xlabel=FALSE))
                    p <- p + theme(strip.background=element_blank(),
                                   strip.text=element_blank()) +
                        theme(rect=element_blank())
                    p2 <- ggplot(data=dfr,aes(Loc,nrel)) +
                        geom_point(aes(color=samples),alpha=I(1/2)) +
                        theme(legend.position='bottom')
                    if(max(dfr$nrel)>1.5 | min(dfr$nrel)< (-1.5)){
                        p2 <- p2 + coord_cartesian(
                            ylim=c(min(dfr$nrel),max(dfr$nrel))) +
                            xlab(paste('Chromosome',chromNum," (Mb)",sep='')) +
                            ylab("Log2 relative expression")
                    }
                    else{
                        p2 <- p2 + coord_cartesian(ylim = c(-0.75,0.75)) +
                            xlab(paste('Chromosome',chromNum, "(Mb)",sep='')) +
                            ylab("Log2 relative expression")
                    }
                    p2 <- p2 + guides(fill=guide_legend(title="Samples"))
                    gA <- ggplot_gtable(ggplot_build(p@ggplot))
                    gB <- ggplot_gtable(ggplot_build(p2))
                    maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3])
                    gA$widths[2:3] <- as.list(maxWidth)
                    gB$widths[2:3] <- as.list(maxWidth)
                    arranged <- grid.arrange(gA, gB, ncol=1,heights=c(1,5))
                    print(arranged)
                    dev.off()
                    return(p2)
                }
            }  
        }
        
        if(combine==TRUE){
            wholeSet <- datalist$whole
            selSet <- wholeSet[as.numeric(samples)]
            
            if(chromNum == "ALL"){
                returnlist <- list()
                length(returnlist) <- suppressWarnings(
                    max(as.numeric(wholeSet[[1]]$Chr),na.rm=TRUE))
                for(g in 1:suppressWarnings(
                    max(as.numeric(wholeSet[[1]]$Chr),na.rm=TRUE))){
                    Locs <- NULL
                    rels <- NULL
                    labls <- NULL
                    nrels <- NULL
                    labls2 <- NULL
                    for(v in 1:length(selSet)){
                        daf <- selSet[[v]]
                        daf  <- daf[order(daf$Chr),]
                        sel <- daf$Chr == as.character(g)
                        daf <- daf[sel,]
                        daf <- daf[order(daf$Loc),]
                        Locs <- append(Locs,daf$Loc)
                        rels <- append(rels,daf[,4])
                        
                        nrels <- slidSmooth(rels,k)
                        len <- length(daf$Chr)
                        tmp <- NULL
                        tmp[1:len] <- names(selSet)[v]
                        labls <- append(labls,tmp)
                        labls2 <- append(labls2,tmp)
                        
                    }          
                    
                    raw <- data.frame(Loc=Locs,rel=rels,samples=labls,
                                      stringsAsFactors=FALSE)
                    sliding <- data.frame(Loc=Locs,rel=nrels,samples=labls2,
                                          stringsAsFactors=FALSE)
                    slidingFin <- slidWithGaps(sliding,len)
                    
                    sam <- names(wholeSet)
                    nam <- NULL
                    if(file=="default"){
                        nam <- paste("Chr",as.character(g),"slidMul",sep="")
                    }
                    
                    else{
                        nam <- file
                    }
                    png(nam,1600,1200)
                    message(paste("Writing plot to ",nam,sep=""))
                    message("\n")
                    raw$Loc <- (raw$Loc/1000000)
                    slidingFin$Loc <- (slidingFin$Loc/1000000)
                    if(idiogram==FALSE){
                        pl <- ggplot(raw,aes(Loc,rel,colour=samples)) 
                        pl <- pl + geom_point(alpha=0.5)
                        pl <- pl + geom_line(data=slidingFin,size=2)
                        if(max(slidingFin$rel,na.rm=TRUE)>1.5 | min(
                            slidingFin$rel,na.rm=TRUE)< (-1.5)){
                            pl <- pl + coord_cartesian(
                                ylim = c(min(slidingFin$rel,na.rm=TRUE),
                                         max(slidingFin$rel,na.rm=TRUE)))
                        }
                        else{
                            pl <- pl + coord_cartesian(ylim = c(-0.75,0.75))
                        }
                        pl <- pl + guides(fill=guide_legend(title="Samples"))
                        pl <- pl + xlab(paste('Chromosome ',as.character(g),
                                              " (Mb)",sep="")) +
                            ylab('Relative expression to the median') +
                            labs(title="Sliding plot")
                        print(pl)
                        dev.off()
                        returnlist[[g]] <- pl
                    }
                    
                    if(idiogram==TRUE){
                        data(hg19IdeogramCyto, envir=environment())
                        chr <- paste("chr",as.character(v),sep="")
                        suppressMessages(p <- plotIdeogram(hg19IdeogramCyto,
                                                              chr,xlabel=FALSE))
                        p <- p + theme(strip.background=element_blank(),
                                       strip.text=element_blank()) +
                            theme(rect=element_blank())
                        p2 <- ggplot(data=raw,aes(Loc,rel)) +
                            geom_point(aes(color=samples),alpha=0.5) +
                            geom_line(data=slidingFin,size=2,
                                      aes(colour=samples)) +
                            theme(legend.position='bottom')
                        if(max(slidingFin$rel,na.rm=TRUE)>1.5 | min(
                            slidingFin$rel,na.rm=TRUE)< (-1.5)){
                            p2 <- p2 + coord_cartesian(
                                ylim=c(min(slidingFin$rel,na.rm=TRUE),
                                       max(slidingFin$rel,na.rm=TRUE))) +
                                xlab(paste('Chromosome',g," (Mb)",sep='')) +
                                ylab("Log2 relative expression")
                        }
                        else{
                            p2 <- p2 + coord_cartesian(ylim = c(-0.75,0.75)) +
                                xlab(paste('Chromosome',g," (Mb)",sep='')) +
                                ylab("Log2 relative expression")
                        }
                        p2 <- p2 + guides(fill=guide_legend(title="Samples"))
                        gA <- ggplot_gtable(ggplot_build(p@ggplot))
                        gB <- ggplot_gtable(ggplot_build(p2))
                        maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3])
                        gA$widths[2:3] <- as.list(maxWidth)
                        gB$widths[2:3] <- as.list(maxWidth)
                        arranged <- grid.arrange(gA, gB, ncol=1,heights=c(1,5))
                        print(arranged)
                        dev.off()
                        returnlist[[g]] <- p2
                    }          
                }
                return(returnlist)
            }
            
            else{
                Locs <- NULL
                rels <- NULL
                nrels <- NULL
                labls <- NULL
                labls2 <- NULL
                for(v in 1:length(selSet)){
                    daf <- selSet[[v]]
                    daf  <- daf[order(daf$Chr),]
                    sel <- daf$Chr == as.character(chromNum)
                    daf <- daf[sel,]
                    daf <- daf[order(daf$Loc),]
                    Locs <- append(Locs,daf$Loc)
                    rels <- append(rels,daf[,4])
                    
                    nrels <- slidSmooth(rels,k)
                    
                    tmp <- NULL
                    len <- length(daf$Chr)
                    tmp[1:len] <- names(selSet)[v]
                    labls <- append(labls,tmp)
                    labls2 <- append(labls2,tmp)
                }
                
                raw <- data.frame(Loc=Locs,rel=rels,samples=labls,
                                  stringsAsFactors=FALSE)
                sliding <- data.frame(Loc=Locs,rel=nrels,samples=labls2,
                                      stringsAsFactors=FALSE)
                slidingFin <- slidWithGaps(sliding,len)
                
                sam <- names(selSet)
                nam <- NULL
                if(file=="default"){
                    nam <- paste("Chr",as.character(chromNum),"slidMul",sep="")
                }
                
                else{
                    nam <- file
                }
                png(nam,1600,1200)
                message(paste("Writing plot to ",nam,sep=""))
                message("\n")
                raw$Loc <- (raw$Loc/1000000)
                slidingFin$Loc <- (slidingFin$Loc/1000000)
                if(idiogram == FALSE){
                    pl <- ggplot(raw,aes(Loc,rel,colour=samples)) 
                    pl <- pl + geom_point(alpha=0.5)
                    pl <- pl + geom_line(data=slidingFin,size=2)
                    if(max(slidingFin$rel,na.rm=TRUE)> 1.5 | min(
                        slidingFin$rel,na.rm=TRUE)< (-1.5)){
                        pl <- pl + coord_cartesian(
                            ylim = c(min(slidingFin$rel,na.rm=TRUE),
                                     max(slidingFin$rel,na.rm=TRUE)))
                    }
                    else{
                        pl <- pl + coord_cartesian(ylim = c(-0.75,0.75))
                    } 
                    pl <- pl + guides(fill=guide_legend(title="Samples"))
                    pl <- pl + xlab(paste('Chromosome ',as.character(chromNum),
                                          " (Mb)",sep="")) +
                        ylab('Relative expression to the median')
                    pl <- pl + labs(title="Sliding plot") +
                        xlab(paste("Chromosome ",chromNum," (Mb)",sep="")) +
                        ylab("Log2 relative expression")
                    print(pl)
                    dev.off()
                    return(pl)
                }
                
                if(idiogram == TRUE){
                    data(hg19IdeogramCyto, envir=environment())
                    chr <- paste("chr",as.character(chromNum),sep="")
                    suppressMessages(p <- plotIdeogram(hg19IdeogramCyto,
                                                          chr,xlabel=FALSE))
                    p <- p + theme(strip.background=element_blank(),
                                   strip.text=element_blank()) +
                        theme(rect=element_blank())
                    p2 <- ggplot(data=raw,aes(Loc,rel)) +
                        geom_point(aes(color=samples),alpha=0.5) +
                        geom_line(data=slidingFin,size=1,aes(colour=samples)) +
                        theme(legend.position='bottom')
                    if(max(slidingFin$rel,na.rm=TRUE)>1.5 | min(
                        slidingFin$rel,na.rm=TRUE)< (-1.5)){
                        p2 <- p2 + coord_cartesian(
                            ylim=c(min(slidingFin$rel,na.rm=TRUE),
                                   max(slidingFin$rel,na.rm=TRUE))) +
                            xlab(paste('Chromosome',chromNum," (Mb)",sep='')) +
                            ylab("Log2 relative expression")
                    }
                    else{
                        p2 <- p2 + coord_cartesian(ylim = c(-1.5,1.5)) +
                            xlab(paste('Chromosome',chromNum," (Mb)",sep='')) +
                            ylab("Log2 relative expression")
                    }
                    p2 <- p2 + guides(fill=guide_legend(title="Samples"))
                    gA <- ggplot_gtable(ggplot_build(p@ggplot))
                    gB <- ggplot_gtable(ggplot_build(p2))
                    maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3])
                    gA$widths[2:3] <- as.list(maxWidth)
                    gB$widths[2:3] <- as.list(maxWidth)
                    arranged <- grid.arrange(gA, gB, ncol=1,heights=c(1,5))
                    print(arranged)
                    dev.off()
                    return(p2)
                }        
            }
        }    
    }
}

