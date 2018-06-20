library(rethinking)
library(gtools)
library(tidyr)
setwd("~/Dropbox/Capuchin Monkeys, Coiba National Park/DataCSVs")

#d <- read.csv(file="~/Dropbox/Capuchin Monkeys, Coiba National Park/DataCSVs/stonetoolajpa.csv" , header=TRUE, sep=";")

d <- read.csv(file="jicaron1.csv" , header=TRUE, sep=",")
d <- d %>% drop_na(WEIGHT)
d <- d[d$WEIGHT>100,]

min(d$WEIGHT)
max(d$WEIGHT)
median(d$WEIGHT)
mean(d$WEIGHT)
sd(d$WEIGHT)

m <- map2stan(
    alist(
        WEIGHT ~ dgamma2(mu,scale),
        log(mu) ~ a ,
        a ~ dnorm(0,2),
        scale ~ dexp(0.5)
    ),

       data=d, cores=2 , warmup=2000 , iter=4000 , WAIC=TRUE, constraints=list(scale="lower=0"), 
)

post <- extract.samples(m)

#check prior against model predictions
dens(post$scale , col="red" , xlim=c(0,100))
dens(rexp(n=10000, rate=0.5) , add=TRUE)

dens(post$a , col="red" , xlim=c(0,20))
dens(rnorm(n=10000, mean=0 , sd=2) , add=TRUE)

#begin model plot
png("stone_weights.png",res=300,height=7,width=7, units = 'in')
par(mar=c(4.5,0.5,0.5,0.5))

plot(density(d$WEIGHT, adjust=1.8), xlim=c(0,2000)  , col="white" , ylim=c(-0.0002,0.0027) , main="" , xlab="mean stone tool weight (g)" , cex.lab=1.7, yaxt='n' , ylab="")

for ( i in 1:100 ) {
    curve(dgamma(x , shape=post$a[i] , scale=post$scale[i])  , add=TRUE ,  col=col.alpha("black",0.1))
    #curve(dgamma(x , shape=(post$a[i] + post$b_adult[i]) , scale=post$scale[i] ) , add=TRUE ,  col=col.alpha("blue",0.1) , lty=1)
}

curve(dgamma(x , shape=mean(post$a) , scale=mean(post$scale))  , add=TRUE , col="black" , lw=2)
#curve(dgamma(x , shape=mean((post$a + post$b_adult)) , scale=median(post$scale) ) , add=TRUE , col="blue" , lty=2 , lw=5)


points( d$WEIGHT , rep(-0.0001, nrow(d)) , pch=19 , col="black" , cex=0.5)
dev.off()
write.csv( precis(m)@output , file="tool_weigt_glm_output.csv")
