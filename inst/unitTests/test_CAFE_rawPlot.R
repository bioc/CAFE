test_CAFE_rawPlot <- function(){
  data(CAFE_data)
  a <- CAFE_data
  single = c(1)
  mul = c(1,2)
  
  p1 <- rawPlot(a,single,chromNum=1,idiogram=FALSE)
  p2 <- rawPlot(a,single,chromNum="ALL",idiogram=FALSE)
  p3 <- rawPlot(a,mul,chromNum="ALL",idiogram=TRUE)
  p4 <- rawPlot(a,mul,chromNum=1,idiogram=TRUE)
  
  checkTrue(class(p1)[1]=="gg")
  checkTrue(class(p4)[1]=="gg")
  checkTrue(length(p3)==22)
  checkTrue(length(p2)==22)
  
  for(i in 1:22){
    #checkTrue(class(p3[[i]])[1]=="gg")
    checkTrue(class(p2[[i]])[1]=="gg")
  }
}