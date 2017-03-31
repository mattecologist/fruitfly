#####
# CLIMEX output

##need to update to "sf" once it works cleanly with raster
library(sp)
library(rgdal)
library(raster)
library (maptools)
library (rasterVis)
library (ggplot2)


# Read the csv table output from CLIMEX
spp1975 <- read.csv ("./data/models/ECB_Europe.csv")
sppext <- c(min(spp1975$Longitude), max(spp1975$Longitude), min(spp1975$Latitude), max(spp1975$Latitude))

##threshold for EI values used in binary presence/absence maps.
EI_thresh <- 10

#raster for calling extent for output
# This is the CliMond data that 
r1 <- raster ("./data/CM10_1975H_Bio_V1.2/CM10_1975H_Bio01_V1.2.txt")

# Crop to extent
r1 <- crop(r1, extent(sppext))

#fishnet for table joins
fishnet <- readOGR("./data/CM10_Fishnet_V1.2.shp")
fishnet <- fishnet[which(fishnet$Location %in% spp1975$Location),]

#fishnet <- crop(fishnet, extent(sppext))
#aus <- readOGR("./data/Aus_admin.shp")

spdf1975 <- merge(fishnet, spp1975, by.x="Location")

## Ecoclimatic index raster
#system.time(spp1975R_EI <- rasterize(spdf1975, r1, field="EI", fun='last'))

###################################################

## This code taken from: http://gis.stackexchange.com/questions/213225/processing-vector-to-raster-faster-with-r?newreg=eefcb4ac453d4fd781a30c2f9dd63828

# Load 'parallel' package for support Parallel computation in R
library('parallel')

# Calculate the number of cores
no_cores <- detectCores() - 1

features <- 1:nrow(spdf1975[,])

# Split features in n parts
n <- 2
parts <- split(features, cut(features, n))

# Initiate cluster (after loading all the necessary object to R environment: BRA_adm2, parts, r.raster, n)
cl <- makeCluster(no_cores, type = "FORK")
print(cl)

# Parallelize rasterize function
system.time(rParts <- parLapply(cl = cl, X = 1:n, fun = function(x) rasterize(spdf1975[parts[[x]],], r1, 'EI')))

# Finish
stopCluster(cl)

# Merge all raster parts
spp1975R_EI <- do.call(merge, rParts)

# Plot raster
plot(spp1975R_EI)

## Growth Index Raster
cl <- makeCluster(no_cores, type = "FORK")
system.time(rParts <- parLapply(cl = cl, X = 1:n, fun = function(x) rasterize(spdf1975[parts[[x]],], r1, 'GI')))
stopCluster(cl)
spp1975R_GI <- do.call(merge, rParts)



## Growth Index raster
spp1975R_GI <- rasterize(spdf1975, r1, field="GI", fun='last')

#######################################################################################
### Plotting
#######################################################################################

gplot(spp1975R_EI) +
  geom_raster(aes(fill = value)) +
  scale_x_continuous(expand = c(0,0), limits= c(sppext[1], sppext[2])) +
  scale_y_continuous(expand = c(0,0), limits= c(sppext[3], sppext[4]))+
  labs(x = "Latitude", y="Longitude") +
  scale_fill_continuous(low="white", high="dark red")+
  coord_equal()+
  ggtitle("Ecoclimatic Index - European Corn Borer")
ggsave(file="./maps/ECB_Europe_EI.png", last_plot())

gplot(spp1975R_GI) +
  geom_raster(aes(fill = value)) +
  scale_x_continuous(expand = c(0,0), limits= c(sppext[1], sppext[2])) +
  scale_y_continuous(expand = c(0,0), limits= c(sppext[3], sppext[4]))+
  labs(x = "Latitude", y="Longitude") +
  scale_fill_continuous(low="white", high="dark green")+
  coord_equal()+
  ggtitle("Growth Index - European Corn Borer")
ggsave(file="./maps/ECB_Europe_GI.png", last_plot())

