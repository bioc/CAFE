test_CAFE_chromosomeStats <- function(){
  data(CAFE_data)
  f <- chromosomeStats(CAFE_data,chromNum="ALL",samples=c(1,3))
  f2 <- chromosomeStats(CAFE_data,chromNum="ALL",samples=c(1,3),enrichment="less")
  
  checkTrue(length(f)==length(f2))
  checkTrue(length(f)==22)
  checkTrue(length(f2)==22)
  
  for(i in 1:22){
    checkTrue(f[i]<=1)
    checkTrue(f2[i]<=1)
  }
}