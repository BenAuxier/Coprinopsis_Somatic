library(qtl2)
library(ggplot2)
library(cowplot)

r1 <- read_cross2("/Users/user/Downloads/Coprinopsis_Somatic/som_round1.zip")
map1 <- insert_pseudomarkers(r1$gmap, step=1)
pr1 <- calc_genoprob(r1,map1,error_prob=0.002)
out1 <- scan1(pr1, r1$pheno)

first <- ggplot() + theme_classic()+
  geom_line(aes(y=out1[1:11],x=c(80,115,150,175,200,225,250,275,300,350,400))) +
  labs(x="Position (kb)",y="LOD")

r2 <- read_cross2("/Users/user/Downloads/Coprinopsis_Somatic/som_round2.zip")
map2 <- insert_pseudomarkers(r2$gmap, step=1)
pr2 <- calc_genoprob(r2,map2,error_prob=0.2)
out2 <- scan1(pr2, r2$pheno)

second <- ggplot() + theme_classic()+
  geom_line(aes(y=out[1:5],x=c(200,225,250,275,300))) +
  labs(x="Position (kb)",y="LOD")

svg("fine_mapping.svg",width=4,height=3)
plot_grid(first,second,nrow=2,labels="AUTO")
dev.off()
