test_CAFE_slidPlot <- function(){
  data(CAFE_data)
  a <- CAFE_data
  single = c(1)
  mul = c(1,2)
  
  p1 <- slidPlot(a,single,chromNum=1,combine=FALSE,idiogram=FALSE)
  p2 <- slidPlot(a,single,chromNum="ALL",combine=FALSE,idiogram=FALSE)
  p3 <- slidPlot(a,single,chromNum=1,combine=TRUE,idiogram=FALSE)
  p4 <- slidPlot(a,single,chromNum=1,combine=TRUE,idiogram=TRUE)
  p5 <- slidPlot(a,single,chromNum=1,combine=FALSE,idiogram=TRUE)
  p6 <- slidPlot(a,single,chromNum="ALL",combine=TRUE,idiogram=FALSE)
  p7 <- slidPlot(a,single,chromNum="ALL",combine=TRUE,idiogram=TRUE)
  p8 <- slidPlot(a,single,chromNum="ALL",combine=FALSE,idiogram=TRUE)
  p9 <- slidPlot(a,mul,chromNum=1,combine=FALSE,idiogram=FALSE)
  p10 <- slidPlot(a,mul,chromNum="ALL",combine=FALSE,idiogram=FALSE)
  p11 <- slidPlot(a,mul,chromNum=1,combine=TRUE,idiogram=FALSE)
  p12 <- slidPlot(a,mul,chromNum=1,combine=TRUE,idiogram=TRUE)
  p13 <- slidPlot(a,mul,chromNum=1,combine=FALSE,idiogram=TRUE)
  p14 <- slidPlot(a,mul,chromNum="ALL",combine=TRUE,idiogram=FALSE)
  p15 <- slidPlot(a,mul,chromNum="ALL",combine=TRUE,idiogram=TRUE)
  p16 <- slidPlot(a,mul,chromNum="ALL",combine=FALSE,idiogram=TRUE)
  
  checkTrue(class(p1)[1]=="gg")
  checkTrue(class(p3)[1]=="gg")
  checkTrue(class(p9)[1]=="gg")
  checkTrue(class(p11)[1]=="gg")
  
  checkTrue(length(p2)==22)
  checkTrue(length(p6)==22)
  checkTrue(length(p7)==22)
  checkTrue(length(p8)==22)
  checkTrue(length(p10)==22)
  checkTrue(length(p14)==22)
  checkTrue(length(p15)==22)
  checkTrue(length(p16)==22)
  
  for(i in 1:22){
    checkTrue(class(p2[[i]][1])=="list")
    checkTrue(class(p6[[i]][1])=="list")
    checkTrue(class(p7[[i]][1])=="list")
    checkTrue(class(p8[[i]][1])=="list")
    checkTrue(class(p10[[i]][1])=="list")
    checkTrue(class(p14[[i]][1])=="list")
    checkTrue(class(p15[[i]][1])=="list")
    checkTrue(class(p16[[i]][1])=="list")
  }
}