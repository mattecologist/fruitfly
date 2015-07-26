#FRUIT!

library (maptools)
library (raster)

#masks the "extract" function, may need to be unloaded
detach("package:tidyr", unload=TRUE)

setwd("~/Documents/climex_out")


##flies PA is missing 2 of the species.

species <- c (#"a_fraterculus", 
  "a_ludens", 
  "a_obliqua",
  "b_correcta",
  "b_cucurbitae", 
  "b_dorsalis",
  "b_latifrons",
  "b_tryoni",
  "b_zonata",
  "c_capitata",
  "c_rosa",
  #"r_indifferens",
  "r_pomonella")


fliesPA <- read.csv ("fliesPA.csv")

import <- readShapePoly ("~/Documents/climex_out/fruit/ImportIndx.shp")
export <- readShapePoly ("~/Documents/climex_out/fruit/ExportIndx.shp")
fruit <- readShapePoly ("~/Documents/climex_out/fruit/FruitProd.shp")

# go home Antarctica, you're drunk 
fruit <- fruit[fruit@data$NAME_1!="Antarctica",]


species_out <- data.frame (country=NA, import=NA, export=NA, production=NA, species=NA, PA=NA)

####CAUTION###

#THIS LOOP TAKES A WHILE
# Assumes you ahve the raster outputs for each species e.g. spp_1975.asc

for (spp in species) { cat(spp,'\n')
  ##This creates a point for each grid cell that has an EI > 10
  r1 <- raster (paste0(spp, "_1975.asc", sep=""))
  r2 <- reclassify (r1, c(-Inf, 10, NA, 10, Inf, 1))
  r3 <- rasterToPoints(r2)
  
  ## now extracts the fruit variable at each of those points (for each country the species is present in)
  r4 <- extract (import, r3[,1:2])
  r4.1 <- extract (export, r3[,1:2])
  r4.2 <- extract (fruit, r3[,1:2])
  
  #data frame of the relevant indicies
  r5 <- data.frame(r4$NAME)
  r5$import  <- (r4$ImportIndx)
  r5$export  <- (r4.1$ExportIndx)
  r5$fruit <- (r4.2$Production)
  
  colnames (r5) <- c("country", "import", "export", "production")
  
  # seems to select only unique countries (variables are averaged across countries, so this is fine to reduce
  # to a single variable)
  r6 <- unique(r5)
  
  #replace to spp name when headings fixed...
  spdf <- data.frame (fliesPA$cntry)
  spdf[,1] <- spp 
  spdf[,2] <- (fliesPA$cntry)
  spdf[,3] <- (fliesPA[,spp])
  colnames (spdf) <- c("species","country", "PA")
  
  
  ###
  # This is now a table of the countries that hold *any* suitable climate space, whether the species has been recorded,
  # import, export and fruit production values. **it does exclude some countries for which a presence is recorded**
  ###
  
  
  test <- merge.data.frame(r6, spdf, by.x="country", all.x=TRUE)
  
  test$species[is.na(test$species)] <- spp
  
  species_out <- rbind(species_out, test)
  
}

###### WRITE OUT AT THIS POINT #########

write.csv (species_out, file="species_out.csv")


##BRING BACK IN - starting point to avoid that loop from before.
species_out <- read.csv ("species_out.csv")
species_out$X <- NULL

# Leaving the NAs as NA for now - might need to match up some country names here..
## Need to determine how many important countries get excluded because names are different

species_out <- na.omit(species_out)
import_df <- species_out[species_out$PA==0,]
export_df <- species_out[species_out$PA==1,] 

import_df$export <- NULL
export_df$import <- NULL

import_df$imprtr <- 0
export_df$exprtr <- 0

#this is long winded, but works fine 

#### indices defined here as "log(production) * export" (or import)
import_df$imprtr[import_df$PA==0] <- log(import_df$production[import_df$PA==0])*import_df$import[import_df$PA==0]
export_df$exprtr[export_df$PA==1] <- log(export_df$production[export_df$PA==1])*export_df$export[export_df$PA==1]

require (tidyr)

#EXPORT DATA
#reshape the data.frame so each species is column
export_df <- spread (export_df, species, exprtr)
# OK to us change NA here as we are looking to sum across species
export_df[is.na(export_df)] <- 0

#sum up the export * production values
export_df$sumexprt <- NULL
export_df$sumexprt <- rowSums(export_df[5:15])
#count how many species per country
export_df$nospec <- rowSums(export_df[,5:15] > 0)


## Same for import
import_df <- spread (import_df, species, imprtr)
import_df[is.na(import_df)] <- 0
import_df$sumimprt <- NULL
import_df$sumimprt <- rowSums(import_df[5:15])
import_df$nospec <- rowSums(import_df[,5:15] > 0)


###### Mapping 

## determine the X highest countries at export risk (sources)
imprt_sort <- arrange(import_df, -sumimprt)
exprt_sort <- arrange(export_df, -sumexprt)

## here X = 20
imprt_sort <- imprt_sort[1:20,]
exprt_sort <- exprt_sort[1:20,]

colnames (imprt_sort)[1] <- "NAME"
imprtmap <- merge(fruit, imprt_sort, by.x="NAME")
imprtmap$sumimprt[is.na(imprtmap$sumimprt)] <- 0
imprtmap <- imprtmap[imprtmap$sumimprt > 0,]

colnames (exprt_sort)[1] <- "NAME"
exprtmap <- merge(fruit, exprt_sort, by.x="NAME")
exprtmap$sumexprt[is.na(exprtmap$sumexprt)] <- 0
exprtmap <- exprtmap[exprtmap$sumexprt > 0,]

imprtmap <- fortify (imprtmap)
exprtmap <- fortify (exprtmap)
base <- fortify (fruit)

require(ggplot2)

#ggplot this map up
map <- ggplot() +  
  geom_polygon(data=base, aes(x=long, y=lat, group=group), fill="grey40", 
               colour="grey90", alpha=1)+
  labs(x="", y="", title="")+ #labels
  theme(axis.line = element_line (colour="black"),
        axis.text = element_text (colour="black", size=20),
        legend.text = element_text (colour="black", size=15),
        legend.text.align=0,
        legend.title = element_blank(),
        legend.title.align = 0,
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.title.y = element_text(face="bold", colour="#000000", size=15),
        axis.title.x = element_text(face="bold", colour="#000000", size=15))+
  geom_polygon(data=exprtmap, aes(x=long, y=lat, group=group), fill="red", 
               colour="grey90", alpha=0.5)+
  geom_polygon(data=imprtmap, aes(x=long, y=lat, group=group), fill="blue", 
               colour="grey90", alpha=0.5)+
  
  guides(col = guide_legend(nrow = 5))+
  coord_equal(ratio=1) # square plot to avoid the distortion


pdf (file = "trademap.pdf", paper="special", width=15, height=8)
map
dev.off()
