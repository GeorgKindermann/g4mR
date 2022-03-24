source("https://raw.githubusercontent.com/GeorgKindermann/g4mR/main/g4m.r")
#source("g4m.r")

tf <- tempfile()
download.file("https://user.iiasa.ac.at/~kinder/g4mr/simpleData/g4mSiteData.csv.xz", tf, quiet=TRUE)
x <- read.csv(tf)
unlink(tf)

#lon ..... Longitude - not needed in g4m but needed to display
#lat ..... Latitude  - not needed in g4m but needed to display
#nn ...... Altitude
#P ....... Phosphorous 
#whc ..... Water holding capacity
#N ....... Nitrogen
#pH ...... pH-Value
#S ....... Salinity
#p1-p12 .. Precipitation
#t1-t12 .. Temperature
#r1-r12 .. Radiation

#Integer Position on map
x$x <- as.integer((x$lon+180.25)*2)
x$y <- as.integer((x$lat+90.25)*2)
x$y <- x$y - min(x$y) + 1

dT <- 0 #delta Temperature
mP <- 1 #Precipitation multiplicator

#Calculate MAI
mai <- apply(x, 1, function(.) {
  g4mMai(t=.[paste0("t", 1:12)] + dT, p=.[paste0("p", 1:12)] * mP, r=.[paste0("r", 1:12)], whc=.["whc"], nn=.["nn"], co2=380, N=.["N"], P=.["P"], S=.["S"], pH=.["pH"])
})

#Display map
. <- matrix(NA, max(x$x), max(x$y))
.[cbind(x$x, x$y)] <- apply(mai, 2, max)
#image(.)

n <- ceiling(max(c(.), na.rm=TRUE))
COL <- hcl.colors(n+1, "YlOrRd", rev = TRUE)
image(seq_along(.[,1]), seq_along(.[1,]), z=., col = COL, axes=FALSE, asp=1, xlab="", ylab="")
legend("bottomleft", legend=0:n, fill = COL, xpd = NA, bty="n", title="mai\n[tC/ha/year]")
