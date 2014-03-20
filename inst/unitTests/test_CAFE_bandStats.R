test_CAFE_bandStats <- function(){
  data(CAFE_data)
  f <- bandStats(CAFE_data,chromNum="ALL",samples=c(1,3))
  f2 <- bandStats(CAFE_data,chromNum="ALL",samples=c(1,3),enrichment="less")
  
  checkTrue(length(f)==length(f2))
  checkTrue(length(f[,1])==length(unique(CAFE_data[[1]][[1]]$Band)))
  checkTrue(length(f2[,1])==length(unique(CAFE_data[[1]][[1]]$Band)))
  
  for(i in 1:length(unique(CAFE_data[[1]][[1]]$Band))){
    checkTrue(f[i]<=1)
    checkTrue(f2[i]<=1)
  }
}