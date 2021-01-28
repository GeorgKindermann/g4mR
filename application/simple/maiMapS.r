source("https://raw.githubusercontent.com/GeorgKindermann/g4mR/main/g4m.r")

u2f <- function(url) {
   tf <- tempfile(fileext = ".tif")
   download.file(url, tf, quiet=TRUE)
   tf
}

D <- readRDS(u2f("https://user.iiasa.ac.at/~kinder/g4m/simpleData/g4mSimpleData.RData"))

dT <- 0 #delta Temperature
mP <- 1 #Precipitation multiplicator

mai <- apply(expand.grid(lapply(dim(D$nn), seq_len)), 1, function(i) {
    max(g4mMai(D$t[i[1],i[2],] + dT, D$p[i[1],i[2],] * mP, D$r[i[1],i[2],]
, D$whc[i[1],i[2]], D$nn[i[1],i[2]], D$st[i[1],i[2]], 0, 0.038))
})
mai <- matrix(mai, dim(D$nn))

image(t(mai[nrow(mai):1,]))
