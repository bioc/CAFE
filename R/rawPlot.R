rawPlot <- function(datalist,samples=c(1,2),chromNum=1,idiogram=FALSE,
                    file="default"){
    Loc <- NULL
    val <- NULL
    hg19IdeogramCyto <- NULL
    rel <- NULL
    lab <- NULL
    nrel <- NULL
    if(length(samples)==1){
        wholeSet <- datalist$whole
        
        if(chromNum=='ALL'){
            returnlist <- list()
            length(returnlist) <- suppressWarnings(max(as.numeric(
                wholeSet[[1]]$Chr),na.rm=TRUE))
            for(i in 1:suppressWarnings(max(as.numeric(wholeSet[[1]]$Chr),
                                            na.rm=TRUE))){
                daf <- wholeSet[[samples]]
                daf <- daf[order(daf$Chr),]
                sel <- daf$Chr == as.character(i)
                daf <- daf[sel,]
                nam <- NULL
                if(file=="default"){
                    nam <- paste("Chr",as.character(i),"-",
                                 as.character(samples),sep="")
                }
                
                else{
                    nam <- file
                }
                message(paste("Writing plot to ",nam,sep=""))
                message("\n")
                png(nam,1600,1200)
                daf$Loc <- (daf$Loc/1000000)
                if(idiogram==FALSE){
                    pl <- qplot(Loc,daf[,4],data=daf,geom=c("point"))
                    if(max(daf[,4])>1 | min(daf[,4])< (-1)){
                        pl <- pl + coord_cartesian(ylim = c(min(daf[,4]),
                                                            max(daf[,4]))) +
                            xlab(paste('Chromosome',i," (Mb)",sep="")) +
                            ylab("Log2 relative expression")
                    }
                    else{
                        pl <- pl + coord_cartesian(ylim = c(-1, 1)) +
                            xlab(paste('Chromosome',i," (Mb)",sep="")) +
                            ylab("Log2 relative expression")
                    }          
                    print(pl)
                    dev.off() 
                    returnlist[[i]] <- pl
                }
                
                if(idiogram==TRUE){
                    dfr <- data.frame(Loc=daf$Loc,rel=daf[,4])
                    data(hg19IdeogramCyto, envir=environment())
                    chr <- paste("chr",as.character(i),sep="")
                    suppressMessages(p <- plotIdeogram(hg19IdeogramCyto, 
                                                          chr,xlabel=FALSE))
                    p <- p + theme(strip.background=element_blank(),
                                   strip.text=element_blank()) +
                        theme(rect=element_blank())
                    p2 <- ggplot(data=dfr,aes(Loc,rel)) +
                        geom_point() +
                        theme(legend.position='bottom')
                    if(max(dfr$rel)>1.5 | min(dfr$rel) < (1.5)){
                        p2 <- p2 + coord_cartesian(ylim=c(min(dfr$rel),
                                                          max(dfr$rel))) +
                            xlab(paste('Chromosome',i," (Mb)",sep='')) +
                            ylab("Log2 relative expression")
                    }
                    else{
                        p2 <- p2 + coord_cartesian(ylim=c(-1.5,1.5)) +
                            xlab(paste('Chromosome',i," (Mb)",sep='')) +
                            ylab("Log2 relative expression")
                    }          
                    gA <- ggplot_gtable(ggplot_build(p@ggplot))
                    gB <- ggplot_gtable(ggplot_build(p2))
                    maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3])
                    gA$widths[2:3] <- as.list(maxWidth)
                    gB$widths[2:3] <- as.list(maxWidth)
                    arranged <- grid.arrange(gA, gB, ncol=1,heights=c(1,5))
                    dev.off()
                    
                    returnlist[[i]] <- p2
                }  
            }
            return(returnlist)
        }
        
        else{
            daf <- wholeSet[[samples]]
            daf <- daf[order(daf$Chr),]
            sel <- daf$Chr == as.character(chromNum)
            daf <- daf[sel,]
            nam <- NULL
            if(file=="default"){
                nam <- paste("Chr",as.character(chromNum),"-",
                             as.character(samples),sep="")
            }
            
            else{
                nam <- file
            }
            daf$Loc <- (daf$Loc/1000000)
            png(nam,1600,1200)
            message(paste("Writing plot to ",nam,sep=""))
            message("\n")
            if(idiogram==FALSE){
                pl <- qplot(Loc,daf[,4],data=daf,geom=c("point"))
                if(max(daf[,4])>1.5 | min(daf[,4])< (-1.5)){
                    pl <- pl + coord_cartesian(ylim=c(min(daf[,4]),
                                                      max(daf[,4]))) +
                        xlab(paste('Chromosome',chromNum," (Mb)",sep='')) +
                        ylab("Log2 relative expression")
                }
                pl <- pl + coord_cartesian(ylim=c(-1.5,1.5)) +
                    xlab(paste('Chromosome',chromNum," (Mb)",sep='')) +
                    ylab("Log2 relative expression")
                print(pl)
                dev.off()
                return(pl)
            }
            
            if(idiogram==TRUE){
                dfr <- data.frame(Loc=daf$Loc,rel=daf[,4])
                data(hg19IdeogramCyto, envir=environment())
                chr <- paste("chr",as.character(chromNum),sep="")
                suppressMessages(p <- plotIdeogram(hg19IdeogramCyto, chr,
                                                      xlabel=FALSE))
                p <- p + theme(strip.background=element_blank(),
                               strip.text=element_blank()) +
                    theme(rect=element_blank())
                p2 <- ggplot(data=dfr,aes(Loc,rel)) +geom_point() +
                    theme(legend.position='bottom') 
                if(max(dfr$rel)> 1.5 | min(dfr$rel)< (1.5)){
                    p2 <- p2 + coord_cartesian(ylim=c(min(dfr$rel),
                                                      max(dfr$rel))) +
                        xlab(paste('Chromosome',chromNum," (Mb)",sep='')) +
                        ylab("Log2 relative expression")
                }
                else{
                    p2 <- p2 + coord_cartesian(ylim=c(-1.5,1.5)) +
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
    
    if(length(samples)>1){
        wholeSet <- datalist$whole
        selSet <- wholeSet[as.numeric(samples)]
        if(chromNum=='ALL'){
            returnlist <- list()
            length(returnlist) <- suppressWarnings(max(
                as.numeric(wholeSet[[1]]$Chr),na.rm=TRUE))
            for(r in 1:suppressWarnings(max(
                as.numeric(wholeSet[[1]]$Chr),na.rm=TRUE))){
                Locs <- NULL
                rels <- NULL
                labls <- NULL     
                
                for(v in 1:length(selSet)){
                    daf <- selSet[[v]]
                    daf <- daf[order(daf$Chr),]
                    sel <- daf$Chr == as.character(r)
                    daf <- daf[sel,]
                    len <- length(daf$Chr)
                    Locs <- append(Locs,daf$Loc)
                    rels <- append(rels,daf[,4])
                    tmp <- NULL
                    tmp[1:len] <- names(selSet)[v]
                    labls <- append(labls,tmp)
                    nreltmp <- NULL
                }
                nam <- NULL
                if(file=="default"){
                    nam <- paste("Chr",as.character(r),"Mul",sep="")
                }
                
                else{
                    nam <- file
                }
                png(nam,1600,1200)
                message(paste("Writing plot to ",nam,sep=""))
                message("\n")
                dfr <- data.frame(Loc=Locs,rel=rels,samples=labls,
                                  stringsAsFactors=FALSE)
                dfr$Loc <-  (dfr$Loc/1000000)
                if(idiogram==FALSE){
                    
                    pl <- qplot(Loc,rel,data=dfr,colour=lab,geom=c("point"))
                    if(max(dfr$rel)> 1.5 | min(dfr$rel)< (1.5)){
                        pl <- pl + coord_cartesian(ylim = c(min(dfr$rel), 
                                                            max(dfr$rel))) +
                            xlab(paste('Chromosome',r," (Mb)",sep="")) +
                            ylab("Log2 relative expression")
                    }
                    else{
                        pl <- pl + coord_cartesian(ylim = c(-1.5, 1.5)) +
                            xlab(paste('Chromosome',r," (Mb)",sep="")) +
                            ylab("Log2 relative expression")
                    }
                    pl <- pl + guides(fill=guide_legend(title="Samples"))
                    print(pl)
                    dev.off()
                    returnlist[[r]] <- pl
                }
                
                if(idiogram==TRUE){
                    data(hg19IdeogramCyto, envir=environment())
                    chr <- paste("chr",as.character(r),sep="")
                    suppressMessages(p <- plotIdeogram(hg19IdeogramCyto,
                                                          chr,xlabel=FALSE))
                    p <- p + theme(strip.background=element_blank(),
                                   strip.text=element_blank()) +
                        theme(rect=element_blank())
                    p2 <- ggplot(data=dfr,aes(Loc,rel)) +
                        geom_point(aes(colour=samples),alpha=I(1/2)) +
                        theme(legend.position='bottom')
                    if(max(dfr$rel)>1.5 | min(dfr$rel)< (1.5)){
                        p2 <- p2 + coord_cartesian(ylim=c(min(dfr$rel),
                                                          max(dfr$rel))) +
                            xlab(paste('Chromosome',r," (Mb)",sep='')) +
                            ylab("Log2 relative expression")
                    }
                    else{
                        p2 <- p2 + coord_cartesian(ylim=c(-1.5,1.5)) +
                            xlab(paste('Chromosome',r," (Mb)",sep='')) +
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
                    returnlist[[r]] <- p2
                }
            }
            return(returnlist)
        }
        
        else{
            
            Locs <- NULL
            rels <- NULL
            labls <- NULL
            
            for(v in 1:length(selSet)){
                daf <- selSet[[v]]
                daf <- daf[order(daf$Chr),]
                sel <- daf$Chr == as.character(chromNum)
                daf <- daf[sel,]
                len <- length(daf$Chr)
                Locs <- append(Locs,daf$Loc)
                rels <- append(rels,daf[,4])
                tmp <- NULL
                tmp[1:len] <- names(selSet)[v]
                labls <- append(labls,tmp)
                nreltmp <- NULL
            }
            
            nam <- NULL
            if(file=="default"){
                nam <- paste("Chr",as.character(chromNum),"Mul",sep="")
            }
            
            else{
                nam <- file
            }
            png(nam,1600,1200)
            message(paste("Writing plot to ",nam,sep=""))
            message("\n")
            dfr <- data.frame(Loc=Locs,rel=rels,samples=labls,
                              stringsAsFactors=FALSE)
            dfr$Loc <- (dfr$Loc/1000000)
            if(idiogram==FALSE){        
                pl <- qplot(Loc,rel,data=dfr,colour=samples,geom=c("point"))
                if(max(dfr$rel)> 1.5 | min(dfr$rel)< (1.5)){
                    pl <- pl + coord_cartesian(ylim = c(min(dfr$rel),
                                                        max(dfr$rel))) +
                        xlab(paste('Chromosome ',chromNum," (Mb)",sep="")) +
                        ylab("Log2 relative expression")
                }
                else{
                    pl <- pl + coord_cartesian(ylim = c(-1.5, 1.5)) +
                        xlab(paste('Chromosome ',chromNum, " (Mb)",sep="")) +
                        ylab("Log2 relative expression")
                }
                pl <- pl + guides(fill=guide_legend(title="Samples"))
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
                p2 <- ggplot(data=dfr,aes(Loc,rel)) +
                    geom_point(aes(colour=samples),alpha=I(1/2)) +
                    theme(legend.position='bottom')
                if(max(dfr$rel)>1.5 | min(dfr$rel)< (1.5)){
                    p2 <- p2 + coord_cartesian(ylim=c(min(dfr$rel),
                                                      max(dfr$rel))) +
                        xlab(paste('Chromosome',chromNum," (Mb)",sep='')) +
                        ylab("Log2 relative expression")
                }
                else{
                    p2 <- p2 + coord_cartesian(ylim=c(-1.5,1.5)) +
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

