#### Tick and possible reservoirs modeling ##########
options(java.parameters = "-Xmx8000m")
library("raster")
library("maptools")
library("spThin")
library("dismo")
library("rJava")
library("ENMeval")
library("rgeos")
library("rgdal")
library("humboldt")
library("tcltk")
library("ENMTools")
library("sp")
library("terra")
library("prettymapr")

#set working directory 
wd <- "."
setwd(wd)

# load in locality data for ticks
all_ticks <- read.csv("Data/Locs/allticks.csv", header=T)
infect_ticks <- read.csv("Data/Locs/infectedticks.csv", header=T)

#load in all locality data
#locs <- list.files(path="data/Re_ Research collaboration/", pattern='csv', full.names=TRUE)
#for (i in 1:length(locs)) assign(locs[i], read.csv(locs[i]))

#load evironemntal variables (path="data/ENV/current/", pattern='tif', full.names=TRUE)
env <- list.files(path="Data/ENV/current/", pattern='tif', full.names=TRUE)
world <- stack(env)
projection(world) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
world
plot(world, 1)

#load in California shape file 
Cali <- readOGR(dsn = "Data/ENV/California_State_Boundary-shp/", 
                layer = "f067d2f7-7950-4e16-beba-8972d811599c2020329-1-18infjv.25og")
plot(Cali)

#look at tick localities 
plot(Cali)
points(all_ticks$longitude, all_ticks$latitude, col="black")
points(infect_ticks$longitude, infect_ticks$latitude, col="red")

## function to create a minimum convex polygon 
# simpleMCP written by Jamie Kass spring 2014 
# Rob Boria vetted code, spring 2014, during his second chapter analyses by comparing code to MCP created in ArcGIS
simpleMCP <- function (xy) {
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2])
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  p <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1)))
}

#clip layers by Cali
layers <- crop(world, Cali)
Current <- mask(layers, Cali)

plot(Current[[1]])
points(all_ticks$longitude, all_ticks$latitude, col="black")

# Generate a buffered MCP (here 1.0 degrees)
tick_mcp <- simpleMCP(all_ticks[2:3])
plot(tick_mcp, add=TRUE, col= "red")

tick_buff <- gBuffer(tick_mcp, width=1.0)
plot(tick_buff, add=T)

tick_layers <- crop(Current, tick_buff)
tick_current <- mask(tick_layers, tick_buff)


#########################################################################
###########Code for uninfected tick modeling starts here ###########################
#########################################################################

#all ticks model evaluation 
tune.args <- list(fc = c("L", "LQ", "H", "LQH", 'LQHT'), rm = 1:6)
model_tick <- ENMevaluate(all_ticks[2:3], tick_current, algorithm = "maxent.jar",
                       tune.args = tune.args, partitions = "block", overlap=F, 
                       updateProgress = F, numCores = 4)
## look at the results 
tick_result <- model_tick@results

write.csv(tick_result, "Output/Eval/All_tick_results.csv", quote=FALSE, row.names=FALSE)

##Generating final models for all ticks 
## Using the top performing model (for now)  ##
#projected to SR; #PM_ best settings were Hinge; RM = 5
args=c("noaddsamplestobackground","noautofeature", "noproduct", "nolinear", 
       "noquadratic", "nothreshold","betamultiplier=5.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#Run model and precit to california
tick_model <- maxent(tick_current, all_ticks[2:3], args = args, 
                   path=paste(getwd(),'/Output/ENM/all_ticks/', sep=''))
tick_predict <- predict(Current, tick_model, args=pred.args) 
plot(tick_predict)
points(all_ticks$longitude, all_ticks$latitude, pch=16, col="black", cex=.35)

#plot with shapefile

plot(tick_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(Cali, add=T)
points(all_ticks$longitude, all_ticks$latitude, pch=16, col="black", cex=.35)

#########################################################################
###########Code for infected tick modeling starts here ##################
#########################################################################

#all ticks model evaluation 
tune.args <- list(fc = c("L", "LQ", "H", "LQH", 'LQHT'), rm = 1:6)
model_infect <- ENMevaluate(infect_ticks[2:3], tick_current, algorithm = "maxent.jar",
                          tune.args = tune.args, partitions = "block", overlap=F, 
                          updateProgress = F, numCores = 4)
## look at the results 
infect_result <- model_infect@results

write.csv(infect_result, "Output/Eval/Infected_tick_results.csv", quote=FALSE, row.names=FALSE)

##Generating final models for all ticks 
## Using the top performing model (for now)  ##
#projected to SR; #PM_ best settings were Hinge; RM = 3
args=c("noaddsamplestobackground","noautofeature", "noproduct", "nolinear",
       "noquadratic", "nothreshold","betamultiplier=3.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#Run model and precit to california
infect_model <- maxent(tick_current, infect_ticks[2:3], args = args, 
                     path=paste(getwd(),'/Output/ENM/infected_ticks/', sep=''))
infect_predict <- predict(Current, infect_model, args=pred.args) 
plot(infect_predict)
points(infect_ticks$longitude, infect_ticks$latitude, pch=16, col="black", cex=.35)

#plot with shapefile

plot(infect_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(Cali, add=T)
points(infect_ticks$longitude, infect_ticks$latitude, pch=16, col="black", cex=.35)

pdf("Figures/tick_ENMs.pdf")
par(mfrow=c(1,2))
plot(tick_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(Cali, add=T)
addnortharrow(pos = "topright", padin = c(0.15, 0.15), scale = .25, lwd = 1, 
              border = "black", cols = c("white", "black"), text.col = "black")
addscalebar(plotunit = "latlon", plotepsg = "EPSG:4326", widthhint = 0.25, 
            unitcategory = "metric", style = "bar", bar.cols = c("black", "white"),
            lwd = 1, linecol = "black", tick.cex = 0.5, labelpadin = 0.08, 
            label.cex = 0.6, label.col = "black", pos = "bottomleft")
points(all_ticks$longitude, all_ticks$latitude, pch=16, col="black", cex=.35)

plot(infect_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(Cali, add=T)
points(infect_ticks$longitude, infect_ticks$latitude, pch=16, col="black", cex=.35)
addnortharrow(pos = "topright", padin = c(0.15, 0.15), scale = .25, lwd = 1, 
              border = "black", cols = c("white", "black"), text.col = "black")

dev.off()




#########################################################################
#########Code for black-tailed jackrabbits modeling starts here #########
#########################################################################

#load in Lepus californicus
lepus <- read.csv("Data/Locs/lecathinned.csv", header=T)

plot(Cali)
points(lepus$longitude, lepus$latitude, pch=16, col="black", cex=.35)
# Generate a buffered MCP (here 1.0 degrees)
lep_mcp <- simpleMCP(lepus[2:3])
plot(lep_mcp, add=TRUE, col= "red")

lep_buff <- gBuffer(lep_mcp, width=1.0)
plot(lep_buff, add=T)

#clip layers by Cali
lep_layers <- crop(Current, lep_buff)
lep_current <- mask(lep_layers, lep_buff)

plot(lep_current[[1]])
points(lepus$longitude, lepus$latitude, col="black")

#all ticks model evaluation 
tune.args <- list(fc = c("L", "LQ", "H", "LQH", 'LQHT'), rm = 1:6)
model_lepus <- ENMevaluate(lepus[2:3], lep_current, algorithm = "maxent.jar",
                            tune.args = tune.args, partitions = "block", overlap=F, 
                            updateProgress = F, numCores = 4)
## look at the results 
lepus_result <- model_lepus@results

write.csv(lepus_result, "Output/Eval/lepus_results.csv", quote=FALSE, row.names=FALSE)

##Generating final models for all ticks 
## Using the top performing model (for now)  ##
#projected to SR; #PM_ best settings were Hinge; RM = 3
args=c("noaddsamplestobackground","noautofeature", "noproduct", "nolinear",
       "noquadratic", "nothreshold","betamultiplier=3.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#Run model and precit to california
lepus_model <- maxent(lep_current, lepus[2:3], args = args, 
                       path=paste(getwd(),'/Output/ENM/Lepus/', sep=''))
lepus_predict <- predict(Current, lepus_model, args=pred.args) 
plot(lepus_predict)
points(lepus$longitude, lepus$latitude, pch=16, col="black", cex=.35)

#plot with shapefile
plot(lepus_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(Cali, add=T)
points(lepus$longitude, lepus$latitude, pch=16, col="black", cex=.35)


#########################################################################
#########Code for brush rabbits modeling starts here #########
#########################################################################

#load in Sylvilagus bachmani
bach <- read.csv("Data/Locs/sybathinned.csv", header=T)

plot(Cali)
points(bach$longitude, bach$latitude, pch=16, col="black", cex=.35)
# Generate a buffered MCP (here 1.0 degrees)
bach_mcp <- simpleMCP(bach[2:3])
plot(bach_mcp, add=TRUE, col= "red")

bach_buff <- gBuffer(bach_mcp, width=1.0)
plot(bach_buff, add=T)

#clip layers by Cali
bach_layers <- crop(Current, bach_buff)
bach_current <- mask(bach_layers, bach_buff)

plot(bach_current[[1]])
points(bach$longitude, bach$latitude, col="black")

#all ticks model evaluation 
tune.args <- list(fc = c("L", "LQ", "H", "LQH", 'LQHT'), rm = 1:6)
model_bach <- ENMevaluate(bach[2:3], bach_current, algorithm = "maxent.jar",
                           tune.args = tune.args, partitions = "block", overlap=F, 
                           updateProgress = F, numCores = 4)
## look at the results 
bach_result <- model_bach@results

write.csv(bach_result, "Output/Eval/Sbachmani_results.csv", quote=FALSE, row.names=FALSE)

##Generating final models for all ticks 
## Using the top performing model (for now)  ##
#projected to SR; #PM_ best settings were LQ; RM = 1
args=c("noaddsamplestobackground","noautofeature", "noproduct", "nohinge", "threshold",
       "betamultiplier=1.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#Run model and precit to california
bach_model <- maxent(bach_current, bach[2:3], args = args, 
                      path=paste(getwd(),'/Output/ENM/Sbachmani/', sep=''))
bach_predict <- predict(Current, bach_model, args=pred.args) 
plot(bach_predict)
points(bach$longitude, bach$latitude, pch=16, col="black", cex=.35)

#plot with shapefile
plot(bach_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(Cali, add=T)
points(bach$longitude, bach$latitude, pch=16, col="black", cex=.35)

#########################################################################
#########Code for desert cottontails  modeling starts here #########
#########################################################################

#load in Sylvilagus audubonii
audub <- read.csv("Data/Locs/syauthinned.csv", header=T)

plot(Cali)
points(audub$longitude, audub$latitude, pch=16, col="black", cex=.35)
# Generate a buffered MCP (here 1.0 degrees)
audub_mcp <- simpleMCP(audub[2:3])
plot(audub_mcp, add=TRUE, col= "red")

audub_buff <- gBuffer(audub_mcp, width=1.0)
plot(audub_buff, add=T)

#clip layers by Cali
audub_layers <- crop(Current, audub_buff)
audub_current <- mask(audub_layers, audub_buff)

plot(audub_current[[1]])
points(audub$longitude, audub$latitude, col="black")

#all ticks model evaluation 
tune.args <- list(fc = c("L", "LQ", "H", "LQH", 'LQHT'), rm = 1:6)
model_audub <- ENMevaluate(audub[2:3], audub_current, algorithm = "maxent.jar",
                          tune.args = tune.args, partitions = "block", overlap=F, 
                          updateProgress = F, numCores = 4)
## look at the results 
audub_result <- model_audub@results

write.csv(audub_result, "Output/Eval/Saudubonii_results.csv", quote=FALSE, row.names=FALSE)

##Generating final models for all ticks 
## Using the top performing model (for now)  ##
#projected to SR; #PM_ best settings were L; RM = 1
args=c("noaddsamplestobackground","noautofeature", "noproduct","nohinge", "noquadratic",
       "nothreshold", "betamultiplier=1.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#Run model and precit to california
audub_model <- maxent(audub_current, audub[2:3], args = args, 
                     path=paste(getwd(),'/Output/ENM/Saudubonii/', sep=''))
audub_predict <- predict(Current, audub_model, args=pred.args) 
plot(audub_predict)
points(audub$longitude, audub$latitude, pch=16, col="black", cex=.35)

#plot with shapefile
plot(audub_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(Cali, add=T)
points(audub$longitude, audub$latitude, pch=16, col="black", cex=.35)

#########################################################################
#########Code for California ground squirrels  modeling starts here #########
#########################################################################

#load in Otospermophilus beecheyi
beech <- read.csv("Data/Locs/otbethinned.csv", header=T)

plot(Cali)
points(beech$longitude, beech$latitude, pch=16, col="black", cex=.35)
# Generate a buffered MCP (here 1.0 degrees)
beech_mcp <- simpleMCP(beech[2:3])
plot(beech_mcp, add=TRUE, col= "red")

beech_buff <- gBuffer(beech_mcp, width=1.0)
plot(beech_buff, add=T)

#clip layers by Cali
beech_layers <- crop(Current, beech_buff)
beech_current <- mask(beech_layers, beech_buff)

plot(beech_current[[1]])
points(beech$longitude, beech$latitude, col="black")

#all ticks model evaluation 
tune.args <- list(fc = c("L", "LQ", "H", "LQH", 'LQHT'), rm = 1:6)
model_beech <- ENMevaluate(beech[2:3], beech_current, algorithm = "maxent.jar",
                           tune.args = tune.args, partitions = "block", overlap=F, 
                           updateProgress = F, numCores = 4)
## look at the results 
beech_result <- model_beech@results

write.csv(audub_result, "Output/Eval/Obeecheyi_results.csv", quote=FALSE, row.names=FALSE)

##Generating final models for all ticks 
## Using the top performing model (for now)  ##
#projected to SR; #PM_ best settings were L; RM = 1
args=c("noaddsamplestobackground","noautofeature", "noproduct","nohinge", "noquadratic",
       "nothreshold", "betamultiplier=1.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#Run model and precit to california
beech_model <- maxent(beech_current, beech[2:3], args = args, 
                      path=paste(getwd(),'/Output/ENM/Obeecheyi/', sep=''))
beech_predict <- predict(Current, beech_model, args=pred.args) 
plot(beech_predict)
points(beech$longitude, beech$latitude, pch=16, col="black", cex=.35)

#plot with shapefile
plot(beech_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(Cali, add=T)
points(beech$longitude, beech$latitude, pch=16, col="black", cex=.35)



#########################################################################
#########Code for bushy tailed woodrat  modeling starts here #########
#########################################################################

#load in Neotoma fuscipes 
neo <- read.csv("Data/Locs/Boria_et_al/Neo_fus_sorted.csv", header=T)

n <- read.csv("Data/Locs/nefuthinned.csv", header=T)

plot(Cali)
points(neo$DEC_LONG, neo$DEC_LAT, pch=16, col="black", cex=.35)

#spatial filter in spThin
nf_thin <- thin(loc.data = neo, 
                 lat.col = "DEC_LAT", long.col = "DEC_LONG", 
                 spec.col = "species", 
                 thin.par = 5, reps = 100, 
                 locs.thinned.list.return = TRUE, 
                 write.files = TRUE, 
                 max.files = 5,
                 out.dir = "Data/Locs/Boria_et_al/", out.base = "NF_f", 
                 write.log.file = TRUE,
                 log.file = "Data/Locs/Boria_et_al/NF_thin_log.txt")

#load in spatial thin dataset
nf_F <- read.csv("Data/Locs/Boria_et_al/NF_f_thin1.csv")

plot(Cali)
points(nf_F$DEC_LONG, nf_F$DEC_LAT, pch=16, col="black", cex=.35)

# Generate a buffered MCP (here 1.0 degrees)
neo_mcp <- simpleMCP(nf_F[2:3])
plot(neo_mcp, add=TRUE, col= "red")

neo_buff <- gBuffer(neo_mcp, width=1.0)
plot(neo_buff, add=T)

#clip layers by Cali
neo_layers <- crop(Current, neo_buff)
neo_current <- mask(neo_layers, neo_buff)

plot(neo_current[[1]])
points(nf_F$DEC_LONG, nf_F$DEC_LAT, col="black")

#all ticks model evaluation 
tune.args <- list(fc = c("L", "LQ", "H", "LQH", 'LQHT'), rm = 1:6)
model_neo <- ENMevaluate(nf_F[2:3], neo_current, algorithm = "maxent.jar",
                           tune.args = tune.args, partitions = "block", overlap=F, 
                           updateProgress = F, numCores = 4)
## look at the results 
neo_result <- model_neo@results

write.csv(audub_result, "Output/Eval/Nfuscipes_results.csv", quote=FALSE, row.names=FALSE)

##Generating final models for all ticks 
## Using the top performing model (for now)  ##
#projected to SR; #PM_ best settings were L; RM = 1
args=c("noaddsamplestobackground","noautofeature", "noproduct","nohinge", "noquadratic",
       "nothreshold", "betamultiplier=1.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#Run model and precit to california
neo_model <- maxent(neo_current, nf_F[2:3], args = args, 
                      path=paste(getwd(),'/Output/ENM/Nfuscipes/', sep=''))
neo_predict <- predict(Current, neo_model, args=pred.args) 
plot(neo_predict)
points(nf_F$DEC_LONG, nf_F$DEC_LAT, pch=16, col="black", cex=.35)

#plot with shapefile
plot(neo_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(Cali, add=T)
points(nf_F$DEC_LONG, nf_F$DEC_LAT, pch=16, col="black", cex=.35)

#########################################################################
#########Code for Douglas squirrels  modeling starts here #########
#########################################################################

#load in tamiasciurus douglasii
doug <- read.csv("Data/Locs/tadothinned.csv", header=T)

plot(Cali)
points(doug$longitude, doug$latitude, pch=16, col="black", cex=.35)

#removing a Baja locality 
doug <- doug[c(1,3:32),]

# Generate a buffered MCP (here 1.0 degrees)
doug_mcp <- simpleMCP(doug[2:3])
plot(doug_mcp, add=TRUE, col= "red")

doug_buff <- gBuffer(doug_mcp, width=1.0)
plot(doug_buff, add=T)

#clip layers by Cali
doug_layers <- crop(Current, doug_buff)
doug_current <- mask(doug_layers, doug_buff)

plot(doug_current[[1]])
points(doug$longitude, doug$latitude, col="black")

#all ticks model evaluation 
tune.args <- list(fc = c("L", "LQ", "H", "LQH", 'LQHT'), rm = 1:6)
model_doug <- ENMevaluate(doug[2:3], doug_current, algorithm = "maxent.jar",
                           tune.args = tune.args, partitions = "block", overlap=F, 
                           updateProgress = F, numCores = 4)
## look at the results 
doug_result <- model_doug@results

write.csv(doug_result, "Output/Eval/Tdouglasii_results.csv", quote=FALSE, row.names=FALSE)

##Generating final models for all ticks 
## Using the top performing model (for now)  ##
#projected to SR; #PM_ best settings were LQH; RM = 2
args=c("noaddsamplestobackground","noautofeature", "noproduct",
       "nothreshold", "betamultiplier=2.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#Run model and precit to california
doug_model <- maxent(doug_current, doug[2:3], args = args, 
                      path=paste(getwd(),'/Output/ENM/Tdouglasii/', sep=''))
doug_predict <- predict(Current, doug_model, args=pred.args) 
plot(doug_predict)
points(doug$longitude, doug$latitude, pch=16, col="black", cex=.35)

#plot with shapefile
plot(doug_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(Cali, add=T)
points(doug$longitude, doug$latitude, pch=16, col="black", cex=.35)


#########################################################################
#########Code for deer mouse modeling starts here #########
#########################################################################

#load in Peromyscus maniculatus
man <- read.csv("Data/Locs/pemathinned.csv", header=T)

plot(Cali)
points(man$longitude, man$latitude, pch=16, col="black", cex=.35)

# Generate a buffered MCP (here 1.0 degrees)
man_mcp <- simpleMCP(man[2:3])
plot(man_mcp, add=TRUE, col= "red")

man_buff <- gBuffer(man_mcp, width=1.0)
plot(man_buff, add=T)

#clip layers by Cali
man_layers <- crop(Current, man_buff)
man_current <- mask(man_layers, man_buff)

plot(man_current[[1]])
points(man$longitude, man$latitude, col="black")

#all ticks model evaluation 
tune.args <- list(fc = c("L", "LQ", "H", "LQH", 'LQHT'), rm = 1:6)
model_man <- ENMevaluate(man[2:3], man_current, algorithm = "maxent.jar",
                          tune.args = tune.args, partitions = "block", overlap=F, 
                          updateProgress = F, numCores = 4)
## look at the results 
man_result <- model_man@results

write.csv(man_result, "Output/Eval/Pmaniculatus_results.csv", quote=FALSE, row.names=FALSE)

##Generating final models for all ticks 
## Using the top performing model (for now)  ##
#projected to SR; #PM_ best settings were LQH; RM = 1
args=c("noaddsamplestobackground","noautofeature", "noproduct",
       "nothreshold", "betamultiplier=1.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#Run model and precit to california
man_model <- maxent(man_current, man[2:3], args = args, 
                     path=paste(getwd(),'/Output/ENM/Pmaniculatus/', sep=''))
man_predict <- predict(Current, man_model, args=pred.args) 
plot(man_predict)
points(man$longitude, man$latitude, pch=16, col="black", cex=.35)

#plot with shapefile
plot(man_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(Cali, add=T)
points(man$longitude, man$latitude, pch=16, col="black", cex=.35)

#########################################################################
#########Code for brush mouse modeling starts here #########
#########################################################################

#load in Peromyscus boylii
boy <- read.csv("Data/Locs/pebothinned.csv", header=T)

plot(Cali)
points(boy$longitude, boy$latitude, pch=16, col="black", cex=.35)

# Generate a buffered MCP (here 1.0 degrees)
boy_mcp <- simpleMCP(boy[2:3])
plot(boy_mcp, add=TRUE, col= "red")

boy_buff <- gBuffer(boy_mcp, width=1.0)
plot(boy_buff, add=T)

#clip layers by Cali
boy_layers <- crop(Current, boy_buff)
boy_current <- mask(boy_layers, boy_buff)

plot(boy_current[[1]])
points(boy$longitude, boy$latitude, col="black")

#all ticks model evaluation 
tune.args <- list(fc = c("L", "LQ", "H", "LQH", 'LQHT'), rm = 1:6)
model_boy <- ENMevaluate(boy[2:3], boy_current, algorithm = "maxent.jar",
                         tune.args = tune.args, partitions = "block", overlap=F, 
                         updateProgress = F, numCores = 4)
## look at the results 
boy_result <- model_boy@results

write.csv(boy_result, "Output/Eval/Pboylii_results.csv", quote=FALSE, row.names=FALSE)

##Generating final models for all ticks 
## Using the top performing model (for now)  ##
#projected to SR; #PM_ best settings were L; RM = 2
args=c("noaddsamplestobackground","noautofeature","nohinge","noquadratic", "noproduct",
       "nothreshold", "betamultiplier=2.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#Run model and precit to california
boy_model <- maxent(boy_current, boy[2:3], args = args, 
                    path=paste(getwd(),'/Output/ENM/Pboylii/', sep=''))
boy_predict <- predict(Current, boy_model, args=pred.args) 
plot(boy_predict)
points(boy$longitude, boy$latitude, pch=16, col="black", cex=.35)

#plot with shapefile
plot(boy_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(Cali, add=T)
points(boy$longitude, boy$latitude, pch=16, col="black", cex=.35)

#########################################################################
#########Code for California deermouse modeling starts here #########
#########################################################################

#load in Peromyscus californicus
califor <- read.csv("Data/Locs/pecathinned.csv", header=T)

plot(Cali)
points(califor$longitude, califor$latitude, pch=16, col="black", cex=.35)

# Generate a buffered MCP (here 1.0 degrees)
califor_mcp <- simpleMCP(califor[2:3])
plot(califor_mcp, add=TRUE, col= "red")

califor_buff <- gBuffer(califor_mcp, width=1.0)
plot(califor_buff, add=T)

#clip layers by Cali
califor_layers <- crop(Current, califor_buff)
califor_current <- mask(califor_layers, califor_buff)

plot(califor_current[[1]])
points(califor$longitude, califor$latitude, col="black")

#all ticks model evaluation 
tune.args <- list(fc = c("L", "LQ", "H", "LQH", 'LQHT'), rm = 1:6)
model_califor <- ENMevaluate(califor[2:3], califor_current, algorithm = "maxent.jar",
                         tune.args = tune.args, partitions = "block", overlap=F, 
                         updateProgress = F, numCores = 4)
## look at the results 
califor_result <- model_califor@results

write.csv(califor_result, "Output/Eval/Pcalifornicus_results.csv", quote=FALSE, row.names=FALSE)

##Generating final models for all ticks 
## Using the top performing model (for now)  ##
#projected to SR; #PM_ best settings were L; RM = 1
args=c("noaddsamplestobackground","noautofeature","nohinge","noquadratic", "noproduct",
       "nothreshold", "betamultiplier=1.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#Run model and precit to california
califor_model <- maxent(califor_current, califor[2:3], args = args, 
                    path=paste(getwd(),'/Output/ENM/Pcalifornicus/', sep=''))
califor_predict <- predict(Current, califor_model, args=pred.args) 
plot(califor_predict)
points(califor$longitude, califor$latitude, pch=16, col="black", cex=.35)

#plot with shapefile
plot(califor_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(Cali, add=T)
points(califor$longitude, califor$latitude, pch=16, col="black", cex=.35)

#########################################################################
#########Code for western fence lizard modeling starts here #########
#########################################################################

#load in sceloporus occidentalis
occid <- read.csv("Data/Locs/scocthinned.csv", header=T)

plot(Cali)
points(occid$longitude, occid$latitude, pch=16, col="black", cex=.35)

# Generate a buffered MCP (here 1.0 degrees)
occid_mcp <- simpleMCP(occid[2:3])
plot(occid_mcp, add=TRUE, col= "red")

occid_buff <- gBuffer(occid_mcp, width=1.0)
plot(occid_buff, add=T)

#clip layers by Cali
occid_layers <- crop(Current, occid_buff)
occid_current <- mask(occid_layers, occid_buff)

plot(occid_current[[1]])
points(occid$longitude, occid$latitude, col="black")

#all ticks model evaluation 
tune.args <- list(fc = c("L", "LQ", "H", "LQH", 'LQHT'), rm = 1:6)
model_occid <- ENMevaluate(occid[2:3], occid_current, algorithm = "maxent.jar",
                          tune.args = tune.args, partitions = "block", overlap=F, 
                          updateProgress = F, numCores = 4)
## look at the results 
occid_result <- model_occid@results

write.csv(occid_result, "Output/Eval/Soccidentalis_results.csv", quote=FALSE, row.names=FALSE)

##Generating final models for all ticks 
## Using the top performing model (for now)  ##
#projected to SR; #PM_ best settings were H; RM = 6
args=c("noaddsamplestobackground","noautofeature", "noproduct", "nolinear", "noquadratic", 
       "betamultiplier=6.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#Run model and precit to california
occid_model <- maxent(occid_current, occid[2:3], args = args, 
                     path=paste(getwd(),'/Output/ENM/Soccidentalis/', sep=''))
occid_predict <- predict(Current, occid_model, args=pred.args) 
plot(occid_predict)
points(occid$longitude, occid$latitude, pch=16, col="black", cex=.35)

#plot with shapefile
plot(occid_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(Cali, add=T)
points(occid$longitude, occid$latitude, pch=16, col="black", cex=.35)



####

pdf("Figures/ENMs_hosts.pdf")
par(mfrow=c(4,3))
plot(lepus_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE, main = expression(italic("Lepus californicus")))
plot(Cali, add=T)
addnortharrow(pos = "topright", padin = c(0.15, 0.15), scale = .25, lwd = 1, 
              border = "black", cols = c("white", "black"), text.col = "black")
addscalebar(plotunit = "latlon", plotepsg = "EPSG:4326", widthhint = 0.25, 
            unitcategory = "metric", style = "bar", bar.cols = c("black", "white"),
            lwd = 1, linecol = "black", tick.cex = 0.5, labelpadin = 0.08, 
            label.cex = 0.6, label.col = "black", pos = "bottomleft")
#points(lepus$longitude, lepus$latitude, pch=16, col="black", cex=.35)

plot(bach_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE, main = expression(italic("Sylvilagus bachmani"))) 
plot(Cali, add=T)
#points(bach$longitude, bach$latitude, pch=16, col="black", cex=.35)

plot(audub_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE, main = expression(italic("Sylvilagus audubonii"))) 
plot(Cali, add=T)
#points(audub$longitude, audub$latitude, pch=16, col="black", cex=.35)

plot(beech_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE, main = expression(italic("Otospermophilus beecheyi"))) 
plot(Cali, add=T)
#points(beech$longitude, beech$latitude, pch=16, col="black", cex=.35)

plot(neo_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE, main = expression(italic("Neotoma fuscipes"))) 
plot(Cali, add=T)
#points(nf_F$DEC_LONG, nf_F$DEC_LAT, pch=16, col="black", cex=.35)

plot(doug_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE, main = expression(italic("Tamiasciurus douglasii"))) 
plot(Cali, add=T)
#points(doug$longitude, doug$latitude, pch=16, col="black", cex=.35)

plot(man_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE, main = expression(italic("Peromyscus maniculatus"))) 
plot(Cali, add=T)
#points(man$longitude, man$latitude, pch=16, col="black", cex=.35)

plot(boy_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE, main = expression(italic("Peromyscus boylii"))) 
plot(Cali, add=T)
#points(boy$longitude, boy$latitude, pch=16, col="black", cex=.35)

plot(califor_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE, main = expression(italic("Peromyscus californicus"))) 
plot(Cali, add=T)
#points(califor$longitude, califor$latitude, pch=16, col="black", cex=.35)

plot(occid_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE, main = expression(italic("Sceloporus occidentalis"))) 
plot(Cali, add=T)
#points(occid$longitude, occid$latitude, pch=16, col="black", cex=.35)
dev.off()

################################################################################
##############################Niche Comparisons#################################
################################################################################

#comparing g-space between infected ticks and Lepus californicus, 
#Sylvilagus audubonii, and Otospermophilus beecheyi

#convert to Spatraster 
current_x <- rast(Current)

#Create emtools.species object for infected ticks
infected_ticks <- enmtools.species()
infected_ticks 

infected_ticks  <- enmtools.species(species.name = "occidentalis", 
            presence.points = infect_ticks[,2:3])

infected_ticks$range <- background.raster.buffer(infected_ticks$presence.points[,2:3], 
                                                 100000, mask = current_x)

infected_ticks$background.points <- background.points.buffer(points = infected_ticks$presence.points,
                                                        radius = 100000, n = 1000, mask = current_x[[1]])

infected_ticks <- check.species(infected_ticks)
infected_ticks

interactive.plot(infected_ticks)

#Create emtools.species object for Lepus
lepus_cal <- enmtools.species()
lepus_cal 

lepus_cal  <- enmtools.species(species.name = "californicus", 
                                    presence.points = lepus[,2:3])

lepus_cal$range <- background.raster.buffer(lepus_cal$presence.points, 
                                                 100000, mask = current_x)

lepus_cal$background.points <- background.points.buffer(points = lepus_cal$presence.points,
                               radius = 100000, n = 1000, mask = current_x[[1]])

lepus_cal <- check.species(lepus_cal)
lepus_cal

interactive.plot(lepus_cal)

#compare two niches in geographic space
infect.lep.sym_id <- identity.test(species.1 = infected_ticks, species.2 = lepus_cal, 
                                   env = current_x, type = "mx", nreps = 100)
infect.lep.sym_id

infect.lep.sym <- background.test(species.1 = infected_ticks, species.2 = lepus_cal,
                   env = current_x, type = "mx", nreps = 100, test.type = "asymmetric", ...=c("noaddsamplestobackground","noautofeature", "noproduct", "nohinge", 
                   "nothreshold","betamultiplier=3.0", "responsecurves=true"))

infect.lep.sym


#Create emtools.species object for Sylvilagus audubonii
audub_cal <- enmtools.species()
audub_cal 

audub_cal  <- enmtools.species(species.name = "audubonii", 
                               presence.points = audub[,2:3])

audub_cal$range <- background.raster.buffer(audub_cal$presence.points, 
                                            100000, mask = current_x)

audub_cal$background.points <- background.points.buffer(points = audub_cal$presence.points,
                                                        radius = 100000, n = 1000, mask = current_x[[1]])

audub_cal <- check.species(audub_cal)
audub_cal

interactive.plot(audub_cal)

#compare two niches in geographic space
infect.aud.sym <- background.test(species.1 = infected_ticks, species.2 = audub_cal,
                                  env = current_x, type = "mx", nreps = 100, test.type = "asymmetric", ...=c("noaddsamplestobackground","noautofeature", "noproduct", "nohinge", 
                                                                                                             "nothreshold","betamultiplier=3.0", "responsecurves=true"))

infect.aud.sym


#Create emtools.species object for Otospermophilus beecheyi
beech_cal <- enmtools.species()
beech_cal 

beech_cal  <- enmtools.species(species.name = "beecheyi", 
                               presence.points = beech[,2:3])

beech_cal$range <- background.raster.buffer(beech_cal$presence.points, 
                                            100000, mask = current_x)

beech_cal$background.points <- background.points.buffer(points = beech_cal$presence.points,
                                                        radius = 100000, n = 1000, mask = current_x[[1]])

beech_cal <- check.species(beech_cal)
beech_cal

interactive.plot(beech_cal)

#compare two niches in geographic space
infect.beech.sym <- background.test(species.1 = infected_ticks, species.2 = beech_cal,
                                  env = current_x, type = "mx", nreps = 100, test.type = "asymmetric", ...=c("noaddsamplestobackground","noautofeature", "noproduct", "nohinge", 
                                                                                                             "nothreshold","betamultiplier=3.0", "responsecurves=true"))

infect.beech.sym$empirical.species.1.model





#####Niche models for California and Oregon 

#load in California shape file 
OR<- readOGR(dsn = "Data/ENV/OR/", 
                  layer = "tl_2022_41_sldl")
OR <- aggregate(OR, dissolve = TRUE)
projection(OR) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
plot(OR)
#combine OR and CA
west <- bind(Cali, OR)
plot(west)


#clip layers by Caliw + OR
layers_pro <- crop(world, west)
Current_pro <- mask(layers_pro, west)

plot(Current_pro[[1]])
points(lepus[2:3]$longitude, lepus[2:3]$latitude, col="black")

lep <- read.csv("Data/Locs/OR_locs/LC_f_thin1.csv")
# Generate a buffered MCP (here 1.0 degrees)
leca_mcp <- simpleMCP(lep[2:3])
plot(leca_mcp, add=TRUE, col= "red")

leca_buff <- gBuffer(leca_mcp, width=1.0)
plot(leca_buff, add=T)

leca_layers <- crop(Current_pro, leca_buff)
leca_current <- mask(leca_layers, leca_buff)

##Lepus californicus,
#projected to SR; #PM_ best settings were Hinge; RM = 3
args=c("noaddsamplestobackground","noautofeature", "noproduct", "nolinear",
       "noquadratic", "nothreshold","betamultiplier=3.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#Run model and precit to california
lep_model <- maxent(leca_current, lep[2:3], args = args, 
                      path=paste(getwd(),'/Output/ENM/Lepus/', sep=''))
lep_predict <- predict(Current_pro, lep_model, args=pred.args) 
plot(lep_predict)
points(lep$longitude, lep$latitude, pch=16, col="black", cex=.35)

#plot with shapefile
plot(lep_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(west, add=T)
points(lep$longitude, lep$latitude, pch=16, col="black", cex=.35)


lepus_model <- maxent(lep_current, lepus[2:3], args = args, 
                      path=paste(getwd(),'/Output/ENM/Lepus/', sep=''))
lepus_predict <- predict(Current_pro, lepus_model, args=pred.args) 
plot(lepus_predict)
points(lepus$longitude, lepus$latitude, pch=16, col="black", cex=.35)





##Otospermophilus beecheyi
# Generate a buffered MCP (here 1.0 degrees)
bee <- read.csv("Data/Locs/OR_locs/OB_f_thin1.csv")
plot(west)
otbe_mcp <- simpleMCP(bee[2:3])
plot(otbe_mcp, add=TRUE, col= "red")

otbe_buff <- gBuffer(otbe_mcp, width=1.0)
plot(otbe_buff, add=T)

otbe_layers <- crop(Current_pro, otbe_buff)
otbe_current <- mask(otbe_layers, otbe_buff)

#projected to SR; #PM_ best settings were L; RM = 1
args=c("noaddsamplestobackground","noautofeature", "noproduct","nohinge", "noquadratic",
       "nothreshold", "betamultiplier=1.0", "responsecurves=true") 
pred.args <- c("outputformat=Cloglog", "doclamp=TRUE")

#Run model and precit to california
bee_model <- maxent(otbe_current, bee[2:3], args = args, 
                      path=paste(getwd(),'/Output/ENM/otbe/', sep=''))
bee_predict <- predict(Current_pro, bee_model, args=pred.args) 
plot(bee_predict)
points(bee$longitude, bee$latitude, pch=16, col="black", cex=.35)

#plot with shapefile
plot(bee_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(west, add=T)
points(bee$longitude, bee$latitude, pch=16, col="black", cex=.35)

or_ticks <- read.csv("Data/Locs/OR_tick_locs.csv")

pdf("Figures/ENM_CA_OR.pdf")
par(mfrow=c(1,2))
plot(lep_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(west, add=T)
points(lep$longitude, lep$latitude, pch=16, col="black", cex=.35)
points(or_ticks[2:3], pch=17, col="red", cex=.35)

plot(bee_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(west, add=T)
points(bee$longitude, bee$latitude, pch=16, col="black", cex=.35)
points(or_ticks$longitude, or_ticks$latitude, pch=17, col="black", cex=.35)

dev.off()




otbe_predict <- predict(Current_pro, lepus_model, args=pred.args) 
plot(otbe_predict)
points(beeech$longitude, beech$latitude, pch=16, col="black", cex=.35)

plot(otbe_predict, col=rev(terrain.colors(10)), breaks= seq(0, 1, by = .1), axes = 0,
     frame.plot=0, box = FALSE) 
plot(west, add=T)
points(beech$longitude, beech$latitude, pch=16, col="black", cex=.35)




