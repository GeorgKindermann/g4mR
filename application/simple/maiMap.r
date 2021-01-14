library(raster)
#source("../../g4m.r")
source("https://raw.githubusercontent.com/GeorgKindermann/g4mR/main/g4m.r")

u2f <- function(url) {
   tf <- tempfile(fileext = ".tif")
   download.file(url, tf, quiet=TRUE)
   tf
}

D <- list()
D$t <- simplify2array(lapply(sprintf("%02d", 1:12), function(i) {
    as.matrix(raster(u2f(paste0("https://github.com/GeorgKindermann/g4m/raw/master/application/simple/data/tavg", i, "30.tif"))))
}))
D$p <- simplify2array(lapply(sprintf("%02d", 1:12), function(i) {
    as.matrix(raster(u2f(paste0("https://github.com/GeorgKindermann/g4m/raw/master/application/simple/data/prec", i, "30.tif"))))
}))
D$r <- simplify2array(lapply(sprintf("%02d", 1:12), function(i) {
    as.matrix(raster(u2f(paste0("https://github.com/GeorgKindermann/g4m/raw/master/application/simple/data/srad", i, "30.tif"))))
}))
D$whc <- as.matrix(raster(u2f("https://github.com/GeorgKindermann/g4m/raw/master/application/simple/data/soilAcwG4m30.tif")))
D$nn <- as.matrix(raster(u2f("https://github.com/GeorgKindermann/g4m/raw/master/application/simple/data/demGtopo30.tif")))
D$st <- as.matrix(raster(u2f("https://github.com/GeorgKindermann/g4m/raw/master/application/simple/data/soilFAO90G4m30.tif")))
D$swr <- as.matrix(raster(u2f("https://github.com/GeorgKindermann/g4m/raw/master/application/simple/data/soilSwr30.tif")))

soilRosetta <- findInterval(1:206, c(1,9,18,24,31,38,46,51,57,61,71,78,81,86,92,97,105,114,121,125,131,137,143,148,155,162,170,177,182,999))
soilRosetta[c(63,65,12,196,103,153,106,161,17,62,151)] <- 30:40
soilRosetta <- match(soilRosetta, c(24,6,30,20,31,1,32,25,33,2,21,34,35,16,14,36,17,37,10,38,39,40))
soilRosetta[is.na(soilRosetta)] <- 0L
D$st[is.na(D$st)] <- 1

i <- D$whc < 1
D$whc[i] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,6,1,1,1,1,5,1,1,1,1,1,1,1,1,1,1,6,3,3,3,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,1,1,1,1,5,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,5,5,6,6,5,1,1,1,1,1,1,1,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,1,1,1,1,1,1,1,1,1,5,1,6,6,3,1,1,3,1,1,3,3,6,1,1,1,1,3,1,1,1,1,1,3,3,1,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)[D$st[i]]
D$whc[] <- c(0,150,125,100,75,50,15)[findInterval(D$whc, 1:6)+1]
D$r <- D$r * 1000./(24.*60.*60.)


mai <- apply(expand.grid(lapply(dim(D$nn), seq_len)), 1, function(i) {
    max(g4mMai(D$t[i[1],i[2],], D$p[i[1],i[2],], D$r[i[1],i[2],]
, D$whc[i[1],i[2]], D$nn[i[1],i[2]], soilRosetta[D$st[i[1],i[2]]], 0, 0.038))
})

plot(raster(matrix(mai, dim(D$nn))))


