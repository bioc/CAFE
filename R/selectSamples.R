#Enable user to select which samples to compare against which samples
#possibly in semi-GUI format
#only relevant for Fisher test

cliSubset <- function(datalist,alternative){
    
    message("Which samples do you want to compare (against the rest)?\n")
    nam <- names(datalist$whole)
    message("Would you first like to see which samples you have?\n Y/N:\n")
    y <- readline()
    if(y == "Y" || y == "y"){
        print(as.character(nam))	
    }
    
    message("Specify samples with their number (e.g. '1' for sample 1), separated by comma's\n")
    x <- readline()
    sep <- ","
    x <- strsplit(x,sep)
    x <- x[[1]]
    x <- as.numeric(x)
    subs_o <- NULL
    if(alternative=="greater"){
        subs_o <- datalist$over[x]
    }
    else if(alternative=="less"){
        subs_o <- datalist$under[x]
    }
    #Fisher test requires ENTIRE datalist to compare to
    subs <- list(samples=subs_o,whole=datalist$whole[x])
    return(subs)
}


guiSubset <- function(datalist,alternative){
    #depends on tcl/tk widgets
    #do something
    nam <- names(datalist$whole)
    tt <- tktoplevel()
    scr <- tkscrollbar(tt, repeatinterval=5,
                       command=function(...)tkyview(tl,...))
    tl <- tklistbox(tt,height=10,selectmode="multiple",
                    yscrollcommand=function(...)tkset(scr,...),
                    background="white")
    tkgrid(tklabel(tt,text="select one or multiple samples"))
    tkgrid(tl,scr)
    tkgrid.configure(scr,rowspan=4,sticky="nsw")
    for(i in 1:length(nam)){
        tkinsert(tl,"end",nam[i])
    }
    result <- NULL
    OnOK <- function(){
        choice <- tkcurselection(tl)
        choice <- as.numeric(choice)+1
        tkdestroy(tt)
        subs_o <- NULL
        if(alternative=="greater"){
            subs_o <- datalist$over[choice]
        }
        else if(alternative=="less"){
            subs_o <- datalist$under[choice]
        }  
        subs <- list(samples=subs_o,whole=datalist$whole[choice])
        assign("guiSelectedSet",subs,envir=as.environment(globalenv()))
        return(subs)
    }
    
    #OK.but <- tkbutton(tt,text="OK",command=result<-OnOK)
    OK.but <- tkbutton(tt,text="OK",command=OnOK)
    tkgrid(OK.but)
    tkfocus(tt)
    
    #return(result)
}

