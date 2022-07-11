library(qtl2)

r1 <- read_cross2(file.choose())
map <- insert_pseudomarkers(r1$gmap, step=1)
pr <- calc_genoprob(r1,map,error_prob=0.002)
out <- scan1(pr, r1$pheno)
plot(out,map)


r2 <- read_cross2(file.choose())
map <- insert_pseudomarkers(r2$gmap, step=1)
pr <- calc_genoprob(r2,map,error_prob=0.002)
out <- scan1(pr, r2$pheno)
plot(out,map)
