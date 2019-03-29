### map the emirates in different colours

## install map packages
library(maps)  
library(mapdata)
library(maptools)


## focus map limits on middle east
map("world")
map.axes()

## looks like middle east is long (xlim) c.40-60, lat (ylim) 10-40

x.lim <-c(50,60)
y.lim <-c(20,30)

admap <- map("world",resolution = 1,xlim=x.lim,ylim=y.lim)
admap <- map("world",resolution = 1,xlim=x.lim,ylim=y.lim)
map.axes()

ad.countries <- admap$names[grep("Emirates",admap$names)]
ad.countries <- ad.countries[c(1:4,6)]

cols <- rainbow(length(ad.countries))

for(i in ad.countries){
  map("world",i, add = T, col = cols[which(ad.countries==i)], fill = T)
}
