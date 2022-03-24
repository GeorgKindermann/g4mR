if(!require(Rcpp)) {print("Trying to install Rcpp")
    install.packages("Rcpp")
    if(!require(Rcpp)) stop("Could not install Rcpp\nPlease install Rcpp")
}

g4m <- new.env()

#sourceCpp("g4mInRAsOs.cc", env=g4m)
tf <- tempfile(fileext = ".cc")
download.file("https://raw.githubusercontent.com/GeorgKindermann/g4mR/main/application/rAsOS/g4mInRAsOs.cc", tf, quiet=TRUE)
sourceCpp(tf, env=g4m)
unlink(tf)

tfIn <- tempfile()
download.file("https://user.iiasa.ac.at/~kinder/g4mr/simpleData/g4mSiteData.txt", tfIn, quiet=TRUE)
tfOut <- tempfile(fileext = ".txt")

g4m$g4mMai(tfIn, tfOut)

unlink(tfIn)
x <- read.table(tfOut)
unlink(tfOut)

names(x) <- c("lon", "lat", "mai")
x$x <- as.integer((x$lon+180.25)*2)
x$y <- as.integer((x$lat+90.25)*2)
x$y <- x$y - min(x$y) + 1
. <- matrix(NA, max(x$x), max(x$y))
.[cbind(x$x, x$y)] <- x$mai
image(.)
