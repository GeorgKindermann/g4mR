if(!require(Rcpp)) {print("Trying to install Rcpp")
    install.packages("Rcpp")
    if(!require(Rcpp)) stop("Could not install Rcpp\nPlease install Rcpp")
}

g4m <- new.env()
sourceCpp("mai.cc", env=g4m)
g4m$maiInit()

g4mMai <- function(t,p,r,whc,nn,st,swr,co2) {
#Returns the MeanAnnual increment for:
# 1..Nadel-Evergreen-Tropical, 2..Laub-Evergreen-Tropical, 3..Nadel-Deciduous-Tropical, 4..Laub-Deciduous-Tropical, 5..Nadel-Evergreen-Subtropical, 6..Laub-Evergreen-Subtropical, 7..Nadel-Deciduous-Subtropical, 8..Laub-Deciduous-Subtropical, 9..Nadel-Evergreen-Temperate, 10..Laub-Evergreen-Temperate, 11..Nadel-Deciduous-Temperate, 12..Laub-Deciduous-Temperate, 13..Nadel-Evergreen-Boreal, 14..Laub-Evergreen-Boreal, 15..Nadel-Deciduous-Boreal, 16..Laub-Deciduous-Boreal
#Input:
#t .. Temperature for each month of this year [12]
#p .. Precipitation for each month of this year [12]
#r .. Radiation for each month of this year [12]
#whc .. Water holding capacity
#nn .. altitude
#st .. soil type
#swr .. Soild water regime
#co2 .. Co2 concentration
    if(!missing(t)) g4m$maiSetTemperature(t)
    if(!missing(p)) g4m$maiSetPrecipitation(p)
    if(!missing(r)) g4m$maiSetRadiation(r)
    if(!missing(whc)) g4m$maiSetWhc(whc)
    if(!missing(nn)) g4m$maiSetAltitude(nn)
    if(!missing(st)) g4m$maiSetSoilType(st)
    if(!missing(swr)) g4m$maiSetSwr(swr)
    if(!missing(co2)) g4m$maiSetCo2(co2)
    g4m$maiGet()
}
