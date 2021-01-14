#source("g4m.r")
source("https://raw.githubusercontent.com/GeorgKindermann/g4mR/main/g4m.r")

D <- list(t = c(-5.5, -5.1, -1., 4, 9.4, 14.1, 16.1, 15.1, 11., 6., 5., -3.5)
, p = c(52,43,45,47,57,71,76,75,66,63,64,60)
, r = rep(180, 12), whc = 150, nn = 400, st = 3L, swr = 0, co2 = 0.038)
do.call(g4mMai, D)
