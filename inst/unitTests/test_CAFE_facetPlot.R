test_CAFE_facetPlot <- function(){
  data(CAFE_data)
  a <- CAFE_data
  single = c(1)
  mul = c(1,2)
  
  p1 <- facetPlot(a,single,slid=FALSE,combine=FALSE)
  p2 <- facetPlot(a,single,slid=TRUE,combine=FALSE)
  p3 <- facetPlot(a,single,slid=FALSE,combine=TRUE)
  p4 <- facetPlot(a,single,slid=TRUE,combine=TRUE)  
  p5 <- facetPlot(a,mul,slid=FALSE,combine=FALSE)
  p6 <- facetPlot(a,mul,slid=TRUE,combine=FALSE)
  p7 <- facetPlot(a,mul,slid=FALSE,combine=TRUE)
  p8 <- facetPlot(a,mul,slid=TRUE,combine=TRUE)
  
  checkTrue(class(p1)[1]=="gg")
  checkTrue(class(p2)[1]=="gg")
  checkTrue(class(p3)[1]=="gg")
  checkTrue(class(p4)[1]=="gg")
  checkTrue(class(p5)[1]=="gg")
  checkTrue(class(p6)[1]=="gg")
  checkTrue(class(p7)[1]=="gg")
  checkTrue(class(p8)[1]=="gg")
}