#####
# CLIMEX output

library(sp)
library(rgdal)
library(raster)
library (maptools)
library (rasterVis)
library (grid)
library (ggplot2)
library (dismo)


##threshold for EI values used in binary presence/absence maps.
EI_thresh <- 10

setwd("/media/matt/OS/climex_out/")

## I used a naming convention of "species_1975" and "species_2050" here
## this will just list the unique species in the output directory 
## this list will be called for all analysis

# files <- data.frame(list.files(pattern='.asc', full.names=TRUE))
# species <- as.data.frame(sapply(files,gsub,pattern=".asc",replacement=""))
# species <- as.data.frame(sapply(species,gsub,pattern="_1975",replacement=""))
# species <- as.data.frame(sapply(species,gsub,pattern="_2030",replacement=""))
# species <- as.data.frame(sapply(species,gsub,pattern="_2050",replacement=""))
# species <- as.data.frame(sapply(species,gsub,pattern="_2070",replacement=""))
# species <- as.data.frame(sapply(species,gsub,pattern="_CS",replacement=""))
# species <- as.data.frame(sapply(species,gsub,pattern="_MR",replacement=""))
# species <- as.data.frame(sapply(species,gsub,pattern="./",replacement=""))
# 
# species <- as.data.frame(sapply(species,gsub,pattern="bs.*",replacement=""))
# 
# species <- unique(species)

species <- c (#"a_fraterculus", 
  "a_ludens", 
  "a_obliqua",
  "b_correcta",
  "b_cucurbitae",
 # "b_cucumis",
  "b_dorsalis",
  #"b_jarvisi",
  "b_latifrons",
  #"b_musae",
  #"b_neohumeralis",
"b_tryoni",
"b_zonata",
"c_capitata",
"c_rosa",
"r_indifferens",
"r_pomonella")

#data frame for storing values from tests
species_df <-data.frame(species)
colnames (species_df) <- "Species"

#List for looping
species <- as.list(sapply(species, levels))

#raster for calling extent for output
r1 <- raster ("./fishnet/bio1.asc.txt")

#fishnet for table joins
fishnet <- readShapePoly("./fishnet/CM30_Fishnet_V1.2.shp")

############## FUNCTIONS ################

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#This next bit is most time intensive! Joins the climex outut to a SPDF and then writes out raster.
#only run once!
for (spp in species) { cat(spp,'\n')

  #spp1975 <- read.csv (paste0("./",spp,"_1975.csv"), skip=2)
  spp2030 <- read.csv (paste0("./",spp,"_2050_MR.csv"), skip=2)

  
  #spdf1975 <- merge(fishnet, spp1975, by.x="Location")
  #spp1975R <- rasterize(spdf1975, r1, field="EI", fun='last')
  
  spdf2030 <- merge(fishnet, spp2030, by.x="Location")
  spp2030R <- rasterize(spdf2030, r1, field="EI", fun='last')
  
  #writeRaster(spp1975R, filename=(paste0("./",spp,"_1975.asc")), overwrite=TRUE)
  writeRaster(spp2030R, filename=(paste0("./",spp,"_2050_MR.asc")), overwrite=TRUE)
  
}



#### AVERAGE THE CLIMATE CHANGE PREDICTIONS

for (spp in species) { cat(spp,'\n')
  
CS2030 <-raster (paste0(spp, "_2030_CS.asc"))
CS2050 <-raster (paste0(spp, "_2050_CS.asc"))
CS2070 <-raster (paste0(spp, "_2070_CS.asc"))

MR2030 <-raster (paste0(spp, "_2030_MR.asc"))
MR2050 <-raster (paste0(spp, "_2050_MR.asc"))
MR2070 <-raster (paste0(spp, "_2070_MR.asc"))

R2030 <- CS2030 + MR2030
R2030 <- R2030/2

R2050 <- CS2050 + MR2050
R2050 <- R2050/2

R2070 <- CS2070 + MR2070
R2070 <- R2070/2

writeRaster(R2030, filename=(paste0("./",spp,"_2030.asc")), overwrite=TRUE)
writeRaster(R2050, filename=(paste0("./",spp,"_2050.asc")), overwrite=TRUE)
writeRaster(R2070, filename=(paste0("./",spp,"_2070.asc")), overwrite=TRUE)
}

### SUMMARY MAP PER TIME SLICE

rc1975 <- raster (nrows=280, ncol=720, xmn=-180, xmx=180, ymn=-56, ymx=84)
rc1975 <- 0

rc2030 <- rc1975
rc2050 <- rc1975
rc2070 <- rc1975

for (spp in species) {cat(spp,'\n')
  
  R1975 <- raster(paste0(spp, "_1975.asc"))
  R2030 <- raster(paste0(spp, "_2030.asc"))
  R2050 <- raster(paste0(spp, "_2050.asc"))
  R2070 <- raster(paste0(spp, "_2070.asc"))
  
  rcR.1 <- reclassify (R1975, c(0, EI_thresh, 0, EI_thresh, Inf, 1))
  rcR.2 <- reclassify (R2030, c(0, EI_thresh, 0, EI_thresh, Inf, 1))
  rcR.3 <- reclassify (R2050, c(0, EI_thresh, 0, EI_thresh, Inf, 1))
  rcR.4 <- reclassify (R2070, c(0, EI_thresh, 0, EI_thresh, Inf, 1))
  rc1975 <- rc1975 + rcR.1
  rc2030 <- rc2030 + rcR.2
  rc2050 <- rc2050 + rcR.3
  rc2070 <- rc2070 + rcR.4
}


#grey 90 and red were original 

rc1975[rc1975==0] <- NA
rc2030[rc2030==0] <- NA
rc2050[rc2050==0] <- NA
rc2070[rc2070==0] <- NA


library (maptools)

fruit <- readShapePoly ("/media/matt/OS/climex_out/fruit/FruitProd.shp")
fruit <- fruit[fruit@data$NAME_1!="Antarctica",]

plot75 <- gplot(rc1975) + geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'grey90', high = 'dark blue', na.value='white', limits = c(0, 12)) +
  ggtitle("Baseline")+
  mytheme +
 geom_path(data=fruit, aes(x=long, y=lat, group=group))+
  guides(fill = guide_colorbar(barwidth = 26, barheight = 2, ticks=FALSE, title="Species No.", nbin=7)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+

  coord_equal()



plot2030 <- gplot(rc2030) + geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'grey90', high = 'dark blue', na.value="white") +
  ggtitle("A2 SRES 2030")+
  mytheme +
  geom_path(data=fruit, aes(x=long, y=lat, group=group))+
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_equal()

plot2050 <- gplot(rc2050) + geom_tile(aes(fill = value)) +
  #facet_wrap(~ variable) +
  scale_fill_gradient(low = 'grey90', high = 'dark blue', na.value="white") +
  ggtitle("A2 SRES 2050")+
  mytheme +
  geom_path(data=fruit, aes(x=long, y=lat, group=group))+
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_equal()

plot2070 <- gplot(rc2070) + geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'grey90', high = 'dark blue', na.value="white") +
  ggtitle("A2 SRES 2070")+
  mytheme +
  geom_path(data=fruit, aes(x=long, y=lat, group=group))+
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_equal()

pdf (file = "./revisions/figure1.pdf", paper="special", width=15, height=8)
grid_arrange_shared_legend(plot75, plot2030, plot2050, plot2070)
dev.off()


# vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(2, 2)))
# print(plot75, vp = vplayout(1, 1))
# print(plot2030, vp = vplayout(1, 2))
# print(plot2050, vp = vplayout(2, 1))
# print(plot2070, vp = vplayout(2, 2))

#loading rasters for nicheOverlap tests in dismo

overlap <- species_df
overlap$baseline <- 1
overlap$D2030 <- NA
overlap$D2050 <- NA
overlap$D2070 <- NA

for (spp in species) {cat(spp,'\n')
    
  R1975 <- raster(paste0(spp, "_1975.asc"))
  R2030 <- raster(paste0(spp, "_2030.asc"))
  R2050 <- raster(paste0(spp, "_2050.asc"))
  R2070 <- raster(paste0(spp, "_2070.asc"))
  
  D2030 <- nicheOverlap(R1975, R2030, stat='D')
  D2050 <- nicheOverlap(R1975, R2050, stat='D')
  D2070 <- nicheOverlap(R1975, R2070, stat='D')
  
  ii = which(overlap$Species==spp)
  overlap$D2030[ii] <- D2030
  overlap$D2050[ii] <- D2050
  overlap$D2070[ii] <- D2070
}

library(reshape)

library (RColorBrewer)

matt <- colorRampPalette(brewer.pal(9,"Paired"))(12)

overlap$col <- matt
overlap$shp <- c(1,2,1,2,3,4,5,6,1,2,1,2)

overlap$Species <- c("A. ludens",
                    "A. obliqua",
                    "B. correcta",
                    "Z. cucurbitae",
                    "B. dorsalis",
                    "B. latifrons",
                    "B. tryoni",
                    "B. zonata",
                    "C. capitata",
                    "C. rosa",
                    "R. indifferens",
                    "R. pomonella")

overlap$genus <- "NA"
overlap$genus[1:2] <- "Anastrepha"
overlap$genus[3] <- "Bactrocera"
overlap$genus[4] <- "Zeugodacus"
overlap$genus[5:8] <- "Bactrocera"
overlap$genus[9:10] <- "Ceratitis"
overlap$genus[11:12] <- "Rhagoletis"

mdata <- melt(overlap, id=c("Species", "genus", "col", "shp"))

mytheme.2 <-  theme(#axis.line = element_blank(),
                  #axis.text = element_blank(),
                  #axis.line = element_line (colour="black"),
                  axis.text = element_text (colour="black", size=10),
                  legend.text = element_text (colour="black", size=10),
                  #legend.text = element_blank(),
                  legend.text.align=0,
                  legend.title = element_blank(),
                  legend.title.align = 0,
                  legend.direction = "vertical",
                  legend.position = "right",
                  #axis.title.y = element_blank(),
                  axis.title.y = element_text(face="bold", colour="#000000", size=10),
                  #axis.title.y = element_blank(),
                  axis.title.x = element_blank())



q <- ggplot (mdata, aes(x=variable, y=value), group=Species) +
  geom_point(aes(colour=Species), size=3, alpha=0.8) +
  geom_line(aes(group=Species, colour=Species), show_guide = FALSE)+
  mytheme.2 +
  scale_y_continuous(name="Schoener's D") +
  scale_x_discrete("variable",labels=c("Baseline","2030","2050","2070")) +
  facet_wrap(~ genus, ncol=2)
q

pdf (file = "./revisions/figure3b.pdf", paper="special", width=7, height=7)
q
dev.off()

q <- q + scale_fill_brewer(palette="Paired")

matt <- colorRampPalette(brewer.pal(9,"Paired"))(13)

q <- ggplot (mdata, aes(x=variable, y=value, fill=Species)) +
  geom_point(size=5, alpha=0.8, colour=matt) +
  geom_line(aes(fill=Species), show_guide = FALSE)+
  mytheme.2 +
  ggtitle("Schoener's D")

p <- ggplot(df, aes(x=cond, y=yval, fill=cond)) + geom_bar(stat="identity") +
  geom_point(aes(x=cond,y=yval2, size=10,color=factor(yval2)))


q <- q + scale_fill_manual(values = matt)

pdf (file = "figure3b.pdf", paper="special", width=7, height=5)
q
dev.off()

##examining changes in latitude

for (spp in species) {

spp1975 <- read.csv (paste0("./",spp,"_1975.csv"), skip=2)
spp2030 <- read.csv (paste0("./",spp,"_2030_CS.csv"), skip=2)
spp2050 <- read.csv (paste0("./",spp,"_2050_CS.csv"), skip=2)
spp2070 <- read.csv (paste0("./",spp,"_2070_CS.csv"), skip=2)
#remove zeros

spp1975 <- spp1975[spp1975$EI >EI_thresh,]
spp2030 <- spp2050[spp2030$EI >EI_thresh,]
spp2050 <- spp2050[spp2050$EI >EI_thresh,]
spp2070 <- spp2050[spp2070$EI >EI_thresh,]

spp1975$ID <- "1975"
spp2030$ID <- "2030"
spp2050$ID <- "2050"
spp2070$ID <- "2070"

sppAll<- rbind(spp1975, spp2030, spp2050, spp2070)


m <- ggplot(sppAll, aes(x = Latitude, fill=ID))+
  geom_density(alpha=0.4)+
  xlim(-90, 90)+
  guides(fill=FALSE, colour=FALSE)+
  ggtitle(spp)

assign(paste0(spp,"_density",sep=""), m)

rm(m)

}

plot_list <- Filter(function(x) is(x, "ggplot"), mget(ls()))
multiplot (plotlist=plot_list, cols=3)


