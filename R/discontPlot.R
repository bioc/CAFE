mu <- function(data,x,y){
    u <- mean(data[x:y])
    return(u)
}

FindBestPartition <- function(y,gamma){
    N <- length(y)
    tmp1 <- rep(1,(gamma+1))
    result <- NULL
    if((N-gamma)>1){
        tmp2 <- seq(from=2,to=(N-gamma),by=1)
        result <- append(tmp1,tmp2)
    }
    else{
        result <- tmp1
    }
    return(result)
}

segmentationFromPartition <- function(p,y){ 
    N <- length(y)
    r <- N
    p <- p[2:length(p)]
    l <- p[length(p)]
    x <- y 
    while(r > 0){
        if(length(l+1)!= 0){
            for(t in (l+1):r){
                x[t] <- mu(y,(l+1),r)
            }
        }
        if(length(l)!=0){
            r <- l
            l <- p[(r-1)] 
        }
        
        else{
            r <- 0
        }
        
    }
    return(x)
} 

discontSmooth <- function(y,gamma){
    p <- FindBestPartition(y,gamma)
    x <- segmentationFromPartition(p,y)
    return(x)
}

discontPlot <- function(datalist,samples=c(1,2),chromNum=1,gamma=300,
                        idiogram=FALSE, file="default"){
    Loc <- NULL
    val <- NULL
    hg19IdeogramCyto <- NULL
    rel <- NULL
    lab <- NULL
    nrel <- NULL
    if(length(samples)==1){
        wholeSet <- datalist$whole
        if(chromNum=="ALL"){
            returnlist <- list()
            for(i in 1:suppressWarnings(max(as.numeric(wholeSet[[1]]$Chr),
                                            na.rm=TRUE))){
                sample <- wholeSet[[samples]]
                chrom <- sample[sample$Chr==i,]
                chrom <- chrom[order(chrom$Loc),]
                y <- chrom[,4]
                x <- discontSmooth(y,gamma)
                result <- data.frame(val=x,Loc=chrom$Loc)
                ori <- data.frame(val=y,Loc=chrom$Loc)
                nam <- NULL
                if(file=="default"){
                    nam <- paste("Chr",as.character(i),"disc",
                                 as.character(samples),sep="")
                }
                else{
                    nam <- file
                }
                message(paste("Writing plot to ",nam,sep=""))
                message("\n")
                result$Loc <- (result$Loc/1000000)
                ori$Loc <- (ori$Loc/1000000)
                if(idiogram==FALSE){
                    png(nam,1600,1200)
                    pl <- ggplot(result,aes(Loc,val)) + geom_line(size=2) + 
                        geom_point(data=ori,alpha=I(0.5))
                    if(max(result$val)>1.5 | min(result$val)< (-1.5)){
                        pl <- pl + coord_cartesian(ylim=c(min(result$val),
                                                          max(result$val))) + 
                            xlab(paste("Chromosome ",i," (Mb)",sep="")) + 
                            ylab("Log2 relative expression")
                    }
                    else{
                        pl <- pl + coord_cartesian(ylim=c(-0.75,0.75)) + 
                            xlab(paste("Chromosome ",i," (Mb)",sep="")) + 
                            ylab("Log2 relative expression")
                    }
                    
                    print(pl)
                    dev.off()
                    returnlist[[i]] <- pl
                }
                else{
                    data(hg19IdeogramCyto, envir=environment())
                    chr <- paste("chr",as.character(i),sep="")
                    suppressMessages(p <- plotSingleChrom(hg19IdeogramCyto,chr,
                                                          xlabel=FALSE))
                    png(nam,1600,1200)
                    p <- p + theme(strip.background=element_blank(),
                                   strip.text=element_blank()) +
                        theme(rect=element_blank())
                    p2 <- ggplot(result,aes(Loc,val)) +geom_line(size=2) +
                        geom_point(data=ori,alpha=I(0.5)) +
                        theme(legend.position='bottom')
                    if(max(result$val)>1.5 | min(result$val)< (-1.5)){
                        p2 <- p2 + coord_cartesian(ylim=c(min(result$val),
                                                          max(result$val))) +
                            xlab(paste('Chromosome',i," (Mb)",sep='')) +
                            ylab("Log2 relative expression")
                    }
                    else{
                        p2 <- p2 + coord_cartesian(ylim=c(-0.75,0.75)) +
                            xlab(paste('Chromosome',i," (Mb)",sep='')) +
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
                    returnlist[[i]] <- p2
                }
            }
            return(returnlist)
        }
        
        else{
            sample <- wholeSet[[samples]]
            chrom <- sample[sample$Chr==chromNum,]
            chrom <- chrom[order(chrom$Loc),]
            y <- chrom[,4]
            x <- discontSmooth(y,gamma)
            result <- data.frame(val=x,Loc=chrom$Loc)
            ori <- data.frame(val=y,Loc=chrom$Loc)
            nam <- NULL
            if(file=="default"){
                nam <- paste("Chr",as.character(chromNum),"disc",
                             as.character(samples),sep="")
            }
            else{
                nam <- file
            }
            message(paste("Writing plot to ",nam,sep=""))
            message("\n")
            result$Loc <- (result$Loc/1000000)
            ori$Loc <- (ori$Loc/1000000)
            if(idiogram==FALSE){
                png(nam,1600,1200)
                pl <- ggplot(result,aes(Loc,val)) + geom_line(size=2) +
                    geom_point(data=ori,alpha=I(0.5))
                if(max(result$val)>1.5 | min(result$val)< (-1.5)){
                    pl <- pl + coord_cartesian(ylim=c(min(result$val),
                                                      max(result$val))) +
                        xlab(paste("Chromosome ",chromNum," (Mb)",sep="")) +
                        ylab("Log2 relative expression") 
                }
                else{
                    pl <- pl + coord_cartesian(ylim=c(-0.75,0.75)) +
                        xlab(paste("Chromosome ",chromNum," (Mb)",sep="")) +
                        ylab("Log2 relative expression")
                }
                
                print(pl)
                dev.off()
                return(pl)
            }
            else{
                data(hg19IdeogramCyto, envir=environment())
                chr <- paste("chr",as.character(chromNum),sep="")
                suppressMessages(p <- plotSingleChrom(hg19IdeogramCyto, chr,
                                                      xlabel=FALSE))
                png(nam,1600,1200)
                p <- p + theme(strip.background=element_blank(),
                               strip.text=element_blank()) +
                    theme(rect=element_blank())
                p2 <- ggplot(result,aes(Loc,val)) +geom_line(size=2) +
                    geom_point(data=ori,alpha=I(0.5)) +
                    theme(legend.position='bottom')
                if(max(result$val)>1.5 | min(result$val)< (-1.5)){
                    p2 <- p2 + coord_cartesian(ylim=c(min(result$val),
                                                      max(result$val))) +
                        xlab(paste('Chromosome',chromNum," (Mb)",sep='')) +
                        ylab("Log2 relative expression") 
                }
                else{
                    p2 <- p2 + coord_cartesian(ylim=c(-0.75,0.75)) +
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
    else{        
        wholeSet <- datalist$whole
        selSet <- wholeSet[as.numeric(samples)]
        if(chromNum=="ALL"){
            returnlist <- list()
            for(a in 1:suppressWarnings(max(as.numeric(wholeSet[[1]]$Chr),
                                            na.rm=TRUE))){
                bigori <- data.frame(val=NULL,Loc=NULL,samples=NULL)
                bigresult <- data.frame(val=NULL,Loc=NULL,samples=NULL)
                for(i in 1:length(selSet)){
                    sample <- selSet[[i]]
                    chrom <- sample[sample$Chr==a,]
                    chrom <- chrom[order(chrom$Loc),]
                    y <- chrom[,4]
                    x <- discontSmooth(y,gamma)
                    result <- data.frame(val=x,Loc=chrom$Loc,
                                         samples=names(selSet)[i])
                    ori <- data.frame(val=y,Loc=chrom$Loc,
                                      samples=names(selSet)[i])
                    bigori <- rbind(bigori,ori)
                    bigresult <- rbind(bigresult,result)
                }
                nam <- NULL
                if(file=="default"){
                    nam <- paste("Chr",as.character(a),"discMul",sep="")
                }
                else{
                    nam <- file
                }
                message(paste("Writing plot to ",nam,sep=""))
                message("\n")
                bigresult$Loc <- (bigresult$Loc/1000000)
                bigori$Loc <- (bigori$Loc/1000000)
                if(idiogram==FALSE){
                    png(nam,1600,1200)
                    pl <- ggplot(bigresult,aes(Loc,val)) +
                        geom_line(aes(colour=samples),size=2) +
                        geom_point(data=bigori,aes(colour=samples),alpha=I(0.5))
                    if(max(bigresult$val)>1.5 | min(bigresult$val)< (-1.5)){
                        pl <- pl + 
                            coord_cartesian(ylim=c(min(bigresult$val),
                                                   max(bigresult$val))) +
                            xlab(paste("Chromosome ",a," (Mb)",sep="")) +
                            ylab("Log2 relative expression")
                    }
                    else{
                        pl <- pl + coord_cartesian(ylim=c(-0.75,0.75)) +
                            xlab(paste("Chromosome ",a," (Mb)",sep="")) +
                            ylab("Log2 relative expression")
                    }
                    pl <- pl + guides(fill=guide_legend(title="Samples"))
                    print(pl)
                    dev.off()
                    returnlist[[a]] <- pl
                }
                else{
                    data(hg19IdeogramCyto, envir=environment())
                    chr <- paste("chr",as.character(a),sep="")
                    suppressMessages(p <- plotSingleChrom(hg19IdeogramCyto, chr,
                                                          xlabel=FALSE))
                    png(nam,1600,1200)
                    p <- p + theme(strip.background=element_blank(),
                                   strip.text=element_blank()) +
                        theme(rect=element_blank())
                    p2 <- ggplot(bigresult,aes(Loc,val)) +
                        geom_line(aes(colour=samples),size=2) +
                        geom_point(data=bigori,aes(colour=samples),
                                   alpha=I(0.5)) +
                        theme(legend.position='bottom')
                    if(max(bigresult$val)>1.5 | min(bigresult$val)< (-1.5)){
                        p2 <- p2 + coord_cartesian(ylim=c(min(bigresult$val),
                                                          max(bigresult$val))) +
                            xlab(paste('Chromosome',a," (Mb)",sep='')) +
                            ylab("Log2 relative expression")
                    }
                    else{
                        p2 <- p2 + coord_cartesian(ylim=c(-0.75,0.75)) +
                            xlab(paste('Chromosome',a," (Mb)",sep='')) +
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
                    returnlist[[a]] <- p2
                }
            }
            return(returnlist)
        }
        else{
            bigori <- data.frame(val=NULL,Loc=NULL,samples=NULL)
            bigresult <- data.frame(val=NULL,Loc=NULL,samples=NULL)
            for(i in 1:length(selSet)){
                sample <- selSet[[i]]
                chrom <- sample[sample$Chr==chromNum,]
                chrom <- chrom[order(chrom$Loc),]
                y <- chrom[,4]
                x <- discontSmooth(y,gamma)
                result <- data.frame(val=x,Loc=chrom$Loc,
                                     samples=names(selSet)[i])
                ori <- data.frame(val=y,Loc=chrom$Loc,
                                  samples=names(selSet)[i])
                bigori <- rbind(bigori,ori)
                bigresult <- rbind(bigresult,result)
            }
            nam <- NULL
            if(file=="default"){
                nam <- paste("Chr",as.character(chromNum),"discMul",sep="")
            }
            else{
                nam <- file
            }
            message(paste("Writing plot to ",nam,sep=""))
            message("\n")
            bigresult$Loc <- (bigresult$Loc/1000000)
            bigori$Loc <- (bigori$Loc/1000000)
            if(idiogram==FALSE){
                png(nam,1600,1200)
                pl <- ggplot(bigresult,aes(Loc,val)) +
                    geom_line(aes(colour=samples),size=2) +
                    geom_point(data=bigori,aes(colour=samples),alpha=I(0.5))
                if(max(bigresult$val)>1.5 | min(bigresult$val)< (-1.5)){
                    pl <- pl + 
                        coord_cartesian(ylim=c(min(bigresult$val),
                                               max(bigresult$val))) +
                        xlab(paste("Chromosome ",chromNum," (Mb)",sep="")) +
                        ylab("Log2 relative expression") 
                }
                else{
                    pl <- pl + coord_cartesian(ylim=c(-0.75,0.75)) +
                        xlab(paste("Chromosome ",chromNum," (Mb)",sep="")) +
                        ylab("Log2 relative expression")
                }
                pl <- pl + guides(fill=guide_legend(title="Samples"))
                print(pl)
                dev.off()
                return(pl)
            }
            else{
                data(hg19IdeogramCyto, envir=environment())
                chr <- paste("chr",as.character(chromNum),sep="")
                suppressMessages(p <- plotSingleChrom(hg19IdeogramCyto,
                                                      chr,xlabel=FALSE))
                png(nam,1600,1200)
                p <- p + theme(strip.background=element_blank(),
                               strip.text=element_blank()) +
                    theme(rect=element_blank())
                p2 <- ggplot(bigresult,aes(Loc,val)) +
                    geom_line(aes(colour=samples),size=2) +
                    geom_point(data=bigori,aes(colour=samples),alpha=I(0.5)) +
                    theme(legend.position='bottom')
                if(max(bigresult$val)>1.5 | min(bigresult$val)< (-1.5)){
                    p2 <- p2 + coord_cartesian(ylim=c(min(bigresult$val),
                                                      max(bigresult$val))) +
                        xlab(paste('Chromosome',chromNum," (Mb)",sep='')) +
                        ylab("Log2 relative expression") 
                }
                else{
                    p2 <- p2 + coord_cartesian(ylim=c(-0.75,0.75)) +
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
