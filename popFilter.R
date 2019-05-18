library(readxl)
library(data.table)
library(magrittr)
population<-read_excel("ukmidyearestimates2017finalversion.xls",sheet = "MYE 5",skip = 4) %>%
  as.data.frame() %>%
  as.data.table()
population<-population[!is.na(Geography1) &
                         !(Geography1 %in% c("London Borough","Non-metropolitan District")),
                       .(region=Name,area_type=Geography1,
              population=`Estimated Population mid-2017`,density=`2017 people per sq. km`)]

getNextCity<-function(topBound,band,densitySens=1){
  temp<-population[population<topBound & population>=(topBound-band),]
  temp[,scoor:=population+densitySens*mean(population)*density/mean(density)]
  temp[scoor==max(scoor)]
}
temp<-getNextCity(700500,700500/7,0.5)
temp2<-temp
temp<-getNextCity(temp[1,population]-700500/(2*7),700500/7,0.5)
temp2<-rbind(temp2,temp)
temp<-getNextCity(temp[1,population]-700500/(2*7),700500/7,0.5)
temp2<-rbind(temp2,temp)
temp<-getNextCity(temp[1,population]-700500/(2*7),700500/7,0.5)
temp2<-rbind(temp2,temp)
temp<-getNextCity(temp[1,population]-700500/(2*7),700500/7,0.5)
temp2<-rbind(temp2,temp)
temp<-getNextCity(temp[1,population]-700500/(2*7),700500/7,0.5)
temp2<-rbind(temp2,temp)
temp<-getNextCity(temp[1,population]-700500/(2*7),700500/7,0.5)
temp2<-rbind(temp2,temp)
temp<-getNextCity(temp[1,population]-700500/(2*7),700500/7,0.5)
temp2<-rbind(temp2,temp)
temp2$area_type<-sapply(strsplit(temp2$area_type, " "),function(i){
  paste0(sapply(i,function(j){substr(j, 1, 1)}),collapse = "")
})
temp2<-temp2[,.(area=paste(region,area_type),population)]
save(temp2,file = "filteredmidyearestimates2017.rData")

