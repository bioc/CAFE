test_CAFE_discontPlot <- function(){
  data(CAFE_data)
  a <- CAFE_data
  single = c(1)
  mul = c(1,2)
  
  p1 <- discontPlot(a,single,chromNum=1,idiogram=FALSE)
  p2 <- discontPlot(a,single,chromNum="ALL",idiogram=FALSE)
  p3 <- discontPlot(a,mul,chromNum="ALL",idiogram=TRUE)
  p4 <- discontPlot(a,mul,chromNum=1,idiogram=TRUE)
  
  checkTrue(class(p1)[1]=="gg")
  
  checkTrue(length(p2)==22)
  checkTrue(length(p3)==22)
  
  for (i in 1:22){
    checkTrue(class(p2[[i]][1])=="list")
    checkTrue(class(p3[[i]][1])=="list")
  }
}