makelevels <- function(CAFE_sample){
    # the problem: we have a vector of chromosomes,
    # containing numbers and letters (mostly numbers)
    # we need to sort these to become '1,2,3,...10,11...,A,B...' 
    # however, R sort will do '1,10,11..19,2,20..,A,B..'
    # one could use mixedsort from the gtools package, 
    # but that results in an extra package dependency which we don't want
    
    # input of this function should be ONE sample from the Datalist$whole list. 
    
    Chr <- CAFE_sample$Chr
    
    all <- unique(Chr)
    numerics <- suppressWarnings(sort(unique(na.omit(as.numeric(Chr)))))
    chars <- sort(all[which(!(all %in% numerics))])
    levels <- c(as.character(numerics),chars)
    return(levels)
}

facetPlot <- function(datalist,samples=c(1,2),slid=FALSE,combine=FALSE,k=1,
                      file="default"){
    Loc <- NULL
    val <- NULL
    rel <- NULL
    lab <- NULL
    nrel <- NULL
    if(slid==FALSE){
        if(length(samples)==1){
            wholeSet <- datalist$whole
            daf <- wholeSet[[samples]]
            b <- makelevels(wholeSet[[1]])
            daf <- daf[order(daf$Chr),]
            daf <- within(daf,Chr <- factor(Chr,levels=b))
            daf <- daf[order(daf$Chr),]
            daf$Loc <- (daf$Loc)/1000000
            nam <- NULL
            if(file=="default"){
                nam <- paste("Facet","-",as.character(samples),sep="")
            }
            
            else{
                nam <- file
            }
            message(paste("Writing plot to ",nam,sep=""))
            names(daf)[4] <- "rel"
            pl <- ggplot(daf,aes(Loc,rel)) + geom_point() + 
                facet_grid(.~Chr, scales="free_x")
            if(max(daf$rel)>1.5 | min(daf$rel)< (-1.5)){
                pl <- pl + coord_cartesian(ylim = c(min(daf$rel),
                                                    max(daf$rel))) +
                    ylab("Log2 relative expression") +
                    xlab("Chromosomal location (Mb)")
            }
            else{
                pl <- pl + coord_cartesian(ylim =c(-0.75,0.75)) +
                    ylab("Log2 relative expression") +
                    xlab("Chromosomal location (Mb)")
            }
            png(nam,5000,600)
            print(pl) 
            dev.off()
            return(pl)
        }
        
        if(length(samples)>1){
            wholeSet <- datalist$whole
            selSet <- wholeSet[as.numeric(samples)]
            
            Locs <- NULL
            rels <- NULL
            labls <- NULL
            Chrs <- NULL
            
            for(v in 1:length(selSet)){
                daf <- selSet[[v]]
                daf <- daf[order(daf$Chr),]
                len <- length(daf$Chr)
                Locs <- append(Locs,daf$Loc)
                rels <- append(rels,daf[,4])
                Chrs <- append(Chrs,daf$Chr)
                tmp <- NULL
                tmp[1:len] <- names(selSet)[v]
                labls <- append(labls,tmp)
                nreltmp <- NULL
            }
            dfr <- data.frame(Loc=Locs,rel=rels,Chr=Chrs,
                              samples=labls,stringsAsFactors=FALSE)
            b <- makelevels(wholeSet[[1]])
            dfr <- dfr[order(dfr$Chr),]
            dfr <- within(dfr,Chr <- factor(Chr,levels=b))
            dfr <- dfr[order(dfr$Chr),]
            dfr$Loc <- (dfr$Loc/1000000)
            nam <- NULL
            if(file=="default"){
                nam <- paste("Facet","Mul",sep="")
            }
            
            else{
                nam <- file
            }
            message(paste("Writing plot to ",nam,sep=""))
            pl <- ggplot(dfr,aes(Loc,rel,colour=samples)) +
                geom_point(alpha=I(1/5)) + facet_grid(.~Chr, scales="free_x")
            if(max(dfr$rel)>1.5 | min(dfr$rel)< (-1.5)){
                pl <- pl + coord_cartesian(ylim = c(min(dfr$rel),
                                                    max(dfr$rel))) +
                    ylab("Log2 relative expression") +
                    xlab("Chromosomal location (Mb)")
            }
            else{
                pl <- pl + coord_cartesian(ylim =c(-0.75,0.75)) +
                    ylab("Log2 relative expression") +
                    xlab("Chromosomal location (Mb)")
            }
            pl <- pl + guides(fill=guide_legend(title="Samples"))
            png(nam,5000,600)
            print(pl)
            dev.off()
            return(pl)
        }    
    }
    
    if(slid==TRUE){
        if(length(samples)==1){
            if(combine==FALSE){
                wholeSet <- datalist$whole
                daf <- wholeSet[[samples]]
                Locs <- NULL
                rels <- NULL
                daf <- daf[order(daf$Chr,daf$Loc),]
                Locs <- daf$Loc
                rels <- daf[,4]        
                nrels <- slidSmooth(rels,k)
                raw <- data.frame(Chr=as.character(daf$Chr),Loc=Locs,rel=rels,
                                  nrel=nrels,stringsAsFactors=FALSE)
                raw$Loc <- (raw$Loc/1000000)
                raw <- na.omit(raw)
                b <- makelevels(wholeSet[[1]])
                raw <- raw[order(raw$Chr),]
                raw <- within(raw,Chr <- factor(Chr,levels=b))
                raw <- raw[order(raw$Chr),]
                nam <- NULL
                if(file=="default"){
                    nam <- paste("FacetSlid",as.character(samples),sep="")
                }
                else{
                    nam <- file
                }
                message(paste("Writing plot to ",nam,sep=""))
                pl <- ggplot(raw,aes(Loc,nrel)) + geom_line() +
                    facet_grid(.~Chr,scales="free_x")
                if(max(raw$nrel)>1.5 | min(raw$nrel)< (-1.5)){
                    pl <- pl + coord_cartesian(ylim = c(min(raw$nrel),
                                                        max(raw$nrel))) +
                        ylab("Log2 relative expression") +
                        xlab("Chromosomal location (Mb)")
                }
                else{
                    pl <- pl + coord_cartesian(ylim =c(-0.75,0.75)) +
                        ylab("Log2 relative expression") +
                        xlab("Chromosomal location (Mb)")
                }
                pl <- pl + guides(fill=FALSE)
                png(nam,5000,600)
                print(pl)
                dev.off()
                return(pl)
            }
            
            if(combine==TRUE){
                wholeSet <- datalist$whole
                daf <- wholeSet[[samples]]
                Locs <- NULL
                rels <- NULL
                daf <- daf[order(daf$Chr,daf$Loc),]
                Locs <- daf$Loc
                rels <- daf[,4]
                nrels <- slidSmooth(rels,k)
                b <- makelevels(wholeSet[[1]])
                raw <- data.frame(Chr=as.character(daf$Chr),Loc=Locs,rel=rels,
                                  nrel=nrels,stringsAsFactors=FALSE)
                raw$Loc <- (raw$Loc/1000000)
                raw <- na.omit(raw)        
                raw <- raw[order(raw$Chr),]
                raw <- within(raw,Chr <- factor(Chr,levels=b))
                raw <- raw[order(raw$Chr),]
                nam <- NULL
                if(file=="default"){
                    nam <- paste("FacetSlid",as.character(samples),sep="")
                }
                else{
                    nam <- file
                }
                message(paste("Writing plot to ",nam,sep=""))
                message("\n")
                png(nam,5000,600)
                pl <- ggplot(raw,aes(Loc,rel)) +
                    facet_grid(. ~ Chr, scales="free_x")
                pl <- pl + geom_point(alpha=I(1/5))
                pl <- pl + geom_line(aes(y=nrel,colour="green")) #TODO somehow give this other colours
                if(max(raw$nrel,na.rm=TRUE)>1.5 | min(raw$nrel,na.rm=TRUE)< (-1.5)){
                    pl <- pl + 
                        coord_cartesian(ylim = c(min(raw$nrel,na.rm=TRUE),
                                                 max(raw$nrel,na.rm=TRUE))) +
                        ylab("Log2 relative expression") +
                        xlab("Chromosomal location (Mb)")
                }
                else{
                    pl <- pl + coord_cartesian(ylim =c(-0.75,0.75)) +
                        ylab("Log2 relative expression") +
                        xlab("Chromosomal location (Mb)")
                }
                pl <- pl + theme(legend.position="none")
                print(pl)
                dev.off()
                return(pl)
            }      
        }
        
        if(length(samples)>1){
            if(combine==FALSE){
                wholeSet <- datalist$whole
                selSet <- wholeSet[as.numeric(samples)]
                Locs <- NULL
                rels <- NULL
                labls <- NULL
                Chrs <- NULL
                nrels <- NULL
                labls2 <- NULL
                for(v in 1:length(selSet)){
                    daf <- selSet[[v]]
                    daf  <- daf[order(daf$Chr,daf$Loc),]
                    Locs <- append(Locs,daf$Loc)
                    rels <- append(rels,daf[,4])
                    Chrs <- append(Chrs,daf$Chr)
                    len <- length(daf$Chr)
                    tmp <- NULL
                    tmp[1:len] <- names(selSet)[v]
                    labls <- append(labls,tmp)
                    labls2 <- append(labls2,tmp)          
                }
                nrels <- slidSmooth(rels,k)
                
                dfr <- data.frame(Loc=Locs,nrel=nrels,Chr=Chrs,
                                  samples=labls,stringsAsFactors=FALSE)
                dfr <- na.omit(dfr)
                b <- makelevels(wholeSet[[1]])
                dfr <- dfr[order(dfr$Chr),]
                dfr <- within(dfr,Chr <- factor(Chr,levels=b))
                dfr <- dfr[order(dfr$Chr),]
                dfr$Loc <- (dfr$Loc/1000000)
                nam <- NULL
                if(file=="default"){
                    nam <- paste("FacetSlidMul",sep="")
                }
                else{
                    nam <- file
                }
                message(paste("Writing plot to ",nam,sep=""))
                message("\n")
                pl <- ggplot(dfr,aes(Loc,nrel,colour=samples)) + geom_line() +
                    facet_grid(.~Chr,scales="free_x")
                if(max(dfr$nrel)>1.5 | min(dfr$nrel)< (-1.5)){
                    pl <- pl + coord_cartesian(ylim = c(min(dfr$nrel), 
                                                        max(dfr$nrel))) +
                        ylab("Log2 relative expression") +
                        xlab("Chromosomal location (Mb)")
                }
                else{
                    pl <- pl + coord_cartesian(ylim =c(-0.75,0.75)) +
                        ylab("Log2 relative expression") +
                        xlab("Chromosomal location (Mb)")
                }
                pl <- pl + guides(fill=guide_legend(title="Samples"))
                png(nam,5000,600)          
                print(pl)
                dev.off()
                return(pl)
            }
            
            if(combine==TRUE){
                wholeSet <- datalist$whole
                selSet <- wholeSet[as.numeric(samples)]
                Locs <- NULL
                rels <- NULL
                labls <- NULL
                Chrs <- NULL
                nrels <- NULL
                labls2 <- NULL
                for(v in 1:length(selSet)){
                    daf <- selSet[[v]]
                    daf  <- daf[order(daf$Chr,daf$Loc),]
                    Locs <- append(Locs,daf$Loc)
                    rels <- append(rels,daf[,4])
                    Chrs <- append(Chrs,daf$Chr)
                    len <- length(daf$Chr)
                    tmp <- NULL
                    tmp[1:len] <- names(selSet)[v]
                    labls <- append(labls,tmp)
                    labls2 <- append(labls2,tmp)          
                }
                nrels <- slidSmooth(rels,k)
                
                raw <- data.frame(Chr=Chrs,Loc=Locs,rel=rels,nrel=nrels,
                                  samples=labls,stringsAsFactors=FALSE)
                sam <- names(wholeSet)[samples]
                b <- makelevels(wholeSet[[1]])
                raw <- raw[order(raw$Chr),]
                raw <- within(raw,Chr <- factor(Chr,levels=b))
                raw <- raw[order(raw$Chr),]
                raw$Loc <- (raw$Loc/1000000)
                nam <- NULL
                if(file=="default"){
                    nam <- paste("FacetSlidMul",sep="")
                }
                else{
                    nam <- file
                }
                message(paste("Writing plot to ",nam,sep=""))
                message("\n")
                png(nam,5000,600)
                pl <- ggplot(raw,aes(Loc,rel,colour=samples)) +
                    facet_grid(. ~ Chr, scales="free_x") 
                pl <- pl + geom_point(alpha=I(1/20))
                pl <- pl + geom_line(aes(y=nrel,colour=samples),size=I(1)) 
                if(max(raw$nrel,na.rm=TRUE)>1.5 | min(raw$nrel,na.rm=TRUE)< (-1.5)){
                    pl <- pl + coord_cartesian(ylim = c(min(raw$nrel,
                                                            na.rm=TRUE),
                                                        max(raw$nrel,
                                                            na.rm=TRUE))) +
                        ylab("Log2 relative expression") +
                        xlab("Chromosomal location (Mb)")
                }
                else{
                    pl <- pl + coord_cartesian(ylim =c(-0.75,0.75)) +
                        ylab("Log2 relative expression") +
                        xlab("Chromosomal location (Mb)")
                }
                pl <- pl + guides(fill=guide_legend(title="Samples"))
                print(pl)
                dev.off()
                return(pl)
            }      
        }    
    }
}