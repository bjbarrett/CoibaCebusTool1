require(RColorBrewer)
require(lubridate)
require(rethinking)
require(tidyverse)
col.pal <- brewer.pal(3, "Dark2")
d1 <- read.csv("~/Dropbox/Capuchin Monkeys, Coiba National Park/DataCSVs/cnp_camtrap_mono_tool_day_jul_2017.csv" , header=TRUE , sep=",")
d2 <- read.csv("~/Dropbox/Capuchin Monkeys, Coiba National Park/DataCSVs/cnp_camtrap_mono_tool_day_dec_2017.csv" , header=TRUE , sep=",")
d3 <- read.csv("~/Dropbox/Capuchin Monkeys, Coiba National Park/DataCSVs/cnp_camtrap_mono_tool_day_mar_2018.csv" , header=TRUE , sep=",")

d <- rbind(d1,d2,d3)
d$camid_index <- as.integer(d$camid)
d$date_index <- as.integer(as.factor(d[,1]))
d$deployment_cam <- paste(d$camid, d$deployment, sep='') 
d$deployment_cam_index <- as.integer(as.factor(d$deployment_cam))
d$island_index <- as.integer(d$island)



pdf("sampling_effort_mar2017tomar2018.pdf" , width=10 , height=6)

plot(d$camid_index ~ d$date_index , xlab="Date (YYYYMMDD)" , ylab="Camera Station", col="white", xaxt="n" ,  yaxt="n" , ylim=c(0.8,max(d$camid_index+0.2)))

for (i in 1:max(d$deployment_cam_index)) {
	segments( min(d$date_index[d$deployment_cam_index==i  & d$deployed==1]) , d$camid_index[d$deployment_cam_index==i] , 
		max(d$date_index[d$deployment_cam_index==i & d$deployed==1]) ,  d$camid_index[d$deployment_cam_index==i] , col=col.pal[min(d$island_index[d$deployment_cam_index==i])])
}

for (i in 1:max(d$deployment_cam_index)) {
	segments( min(d$date_index[d$deployment_cam_index==i & d$camwork==1]) , d$camid_index[d$deployment_cam_index==i] , 
		max(d$date_index[d$deployment_cam_index==i & d$camwork==1]) ,  d$camid_index[d$deployment_cam_index==i] , lw=3 ,  col=col.pal[min(d$island_index[d$deployment_cam_index==i])])
}

#points(d$date_index[d$toolday==1 & d$deployed==1], d$camid_index[d$toolday==1 & d$deployed==1] , pch=18 , cex=0.9 )
points(d$date_index[d$toolday==1 ], d$camid_index[d$toolday==1 ] , pch=18 , cex=0.9 )

xlabs <- c(20170401,20170501,20170601,20170701,20170801,20170901,20171001,20171101,20171201,20180101,20180201,20180301) 
xloc <- sort(unique(d$date_index[d[,1] %in% xlabs]))
axis(1, at=xloc , padj=0 , cex.axis=0.7, labels=FALSE)
axis(1, at=1:max(d$date_index) , padj=0 , cex.axis=0.7, labels=FALSE , tck=-0.01)
text(x = xloc, par("usr")[3] - 0.9, labels = xlabs, srt = 45, pos = 1, xpd = TRUE, cex=0.6)
axis(2, at=1:max(d$camid_index) , padj=0 , cex.axis=0.7, labels=FALSE , tck=-0.01)
text(y = 1:max(d$camid_index), par("usr")[1], labels = sort(unique(d$camid)), srt = 30, pos = 2, xpd = TRUE, cex=0.5)

legend("topleft", c("Coiba" , "Jicaron" , "Rancheria"), pch=15, col=col.pal, box.col=NA, cex=1 )

dev.off()

#set up GLMM dataset
d1 <- d[d$deployment_cam %in% unique(d$deployment_cam[d$toolday==1]),] #cameras where only tooluse was observed
d1$date <- ymd(d1$yyyymmdd) #makes lubridate work well
d1$month <- month(d1$date)
d1$year <- year(d1$date)



d2 <- droplevels(d1 %>% group_by(month, year, camid , deployment_cam, stillcam ) %>% summarize(monodays = sum(monoday), camworkdays=sum(camwork) , tooldays=sum(toolday) , year_month=min(date) ))
d2$camid_index <- as.integer(d2$camid)
d2$camid_index <- as.integer(d2$camid)

d2$toolRPM <- d2$tooldays/d2$camworkdays
d2$monoRPM <- d2$monodays/d2$camworkdays
d2$toolmonoRPM <- d2$tooldays/d2$monodays

d2$year_month<- d2$month + (d2$year*100) 
d2$year_month_index <- as.integer(as.factor(d2$year_month))
d2$deployment_cam_index <- as.integer(as.factor(d2$deployment_cam))

col.pal2 <- brewer.pal(max(d2$camid_index), "Dark2")
d2 <- as.data.frame(d2[d2$camworkdays!=0,])#drop these so model runs


m <- map2stan(
    alist(
        tooldays ~ dbinom( camworkdays , p ),
        logit(p) <- a_ym[year_month_index] + a_c[camid_index] + a_dc[deployment_cam_index] ,
        a_ym[year_month_index] ~ dnorm( a , sigma_ym ),
        a_c[camid_index] ~ dnorm( a , sigma_c ),
        a_dc[deployment_cam_index] ~ dnorm( a , sigma_dc ),

        a ~ dnorm(0,5),
        c(sigma_ym,sigma_c,sigma_dc) ~ dcauchy(0,2)
    ) ,
    data=d2 , warmup=500 , iter=3000 , chains=3 )





write.csv (precis(m, depth=2)@output , "m_output.csv")
post <- extract.samples(m)
post2 <- post$a_ym
for (i in 1:max(d2$year_month_index)){
	post2[,i] <- logistic(post$a + post$a_ym[,i])
}

pdf("camtrapmonth.pdf" , width=10 , height=7)
plot(toolRPM~year_month_index, data=d2 , col=col.pal2[camid_index] , pch=16 , xlab="Month" , ylab="monthly proportion of days where tool use was observed" , xaxt="n")
points( 1:13 , apply(post2, 2, mean) , col="gray44" , pch=18 , cex=2 )
for (i in 1:max(d2$year_month_index)){
	lines ( c(i,i) , HPDI(post2[,i]) , col="gray44")
}
axis(side=1, at=1:13, labels=c("Mar 2017", "April 2017" , "May 2017" , "June 2017" , "July 2017" , "Aug 2017" , "Sept 2017" , "Oct 2017" , "Nov 2017" , "Dec 2017" , "Jan 2018" , "Feb 2018" , "Mar 2018" ) , las =1 , cex.axis=0.75)


post <- extract.samples(m)
post3 <- post$a_c
for (i in 1:max(d2$camid_index)){
	post3[,i] <- logistic(post$a + post$a_c[,i])

}

pdf("between_camera_prob_observing_tool.pdf" , width=8 , height=5)
dens(post3[,1] , col=col.pal2[1] , ylim=c(0,12) , xlim=c(0,1) , lty=1,  xlab="monthly proportion of days where tool use was observed" , ylab="posterior density")
for (i in 1:max(d2$camid_index)){
	dens(post3[,i] , col=col.pal2[i] , add=TRUE , lty=1)
}
points(d2$toolRPM , rep(-0.15 , nrow(d2) )  , col=col.pal2[d2$camid_index])
dev.off()


######using as proportion of days montkeys were observed

d3 <- as.data.frame(d2[d2$monodays!=0,])

m1 <- map2stan(
    alist(
        tooldays ~ dbinom( monodays , p ),
        logit(p) <- a_ym[year_month_index] + a_c[camid_index] + a_dc[deployment_cam_index] ,
        a_ym[year_month_index] ~ dnorm( a , sigma_ym ),
        a_c[camid_index] ~ dnorm( a , sigma_c ),
        a_dc[deployment_cam_index] ~ dnorm( a , sigma_dc ),

        a ~ dnorm(0,5),
        c(sigma_ym,sigma_c,sigma_dc) ~ dcauchy(0,2)
    ) ,
    data=d3 , warmup=1000 , iter=3500 , chains=3  )

write.csv (precis(m1, depth=2)@output , "m1_output.csv")

precis(m1, depth=2)
post <- extract.samples(m1)
post2 <- post$a_ym
for (i in 1:max(d3$year_month_index)){
	post2[,i] <- logistic(post$a + post$a_ym[,i])
}

pdf("camtrapmonthmono.pdf" , width=10 , height=7)
plot(toolmonoRPM~year_month_index, data=d3 , col=col.pal2[camid_index] , pch=16 , xlab="Month" , ylab="monthly proportion of days where tool use was observed" , xaxt="n", cex=1.25 , cex.lab=1.25)
points( 1:13 , apply(post2, 2, mean) , col="black" , pch=18 , cex=2 )
for (i in 1:max(d3$year_month_index)){
	lines ( c(i,i) , HPDI(post2[,i]) , col="black")
}
axis(side=1, at=1:13, labels=c("Mar 2017", "April 2017" , "May 2017" , "June 2017" , "July 2017" , "Aug 2017" , "Sept 2017" , "Oct 2017" , "Nov 2017" , "Dec 2017" , "Jan 2018" , "Feb 2018" , "Mar 2018" ) , las =1 , cex.axis=0.75)
dev.off()

post <- extract.samples(m1)
post3 <- post$a_c
for (i in 1:max(d3$camid_index)){
	post3[,i] <- logistic(post$a + post$a_c[,i])

}

pdf("between_camera_prob_observing_toolmono.pdf" , width=8 , height=5)
dens(post3[,1] , col=col.pal2[1] , ylim=c(0,8) , xlim=c(0,1) , lty=1,  xlab="monthly proportion of days where tool use was observed" , ylab="posterior density")
for (i in 1:max(d3$camid_index)){
	dens(post3[,i] , col=col.pal2[i] , add=TRUE , lty=1)
}
points(d3$toolmonoRPM , rep(-0.15 , nrow(d3) )  , col=col.pal2[d3$camid_index])
dev.off()

##GET SAMPLING EFFORT PER ISLAND
###get sampling effort

d$date <- ymd(d$yyyymmdd) #makes lubridate work well
samp <- rep(0,max(d$deployment_cam_index))
for (i in 1:max(d$deployment_cam_index)) {
    MAX <-  max(d$date[d$deployment_cam_index==i & d$camwork==1])
    MIN <-  min(d$date[d$deployment_cam_index==i & d$camwork==1])
    samp[i] <- as.integer(MAX-MIN)
}
unique(d$deployment_cam_index[d$island=="Coiba"])
sum(samp[unique(d$deployment_cam_index[d$island=="Coiba"])] , na.rm = TRUE)
sum(samp[unique(d$deployment_cam_index[d$island=="Jicaron"])] , na.rm = TRUE)
sum(samp[unique(d$deployment_cam_index[d$island=="Rancheria"])] , na.rm = TRUE)

#proportion of days where monkeys were observed where tool use was observed
length(unique(d$date[d$toolday==1]))/length(unique(d$date[d$monoday==1]))

length(unique(d$date_index[d$toolday==1 &d$island=="Jicaron"]))/length(unique(d$date_index[d$monoday==1 &d$island=="Jicaron"]))
