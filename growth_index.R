library (raster)
library (ggplot2)
library (rasterVis)

host <- raster ("./maps/ECB_global_GI.asc")

parasitoid <- raster ("./maps/TO_global_GI.asc")


## rules
#parasitoid2 <- reclassify (parasitoid, c(-Inf, 10, NA))
#host2 <- reclassify (host, c(-Inf, 10, NA))
difference <- parasitoid2 - host2



# Relative suitability for P. xylostella (RSPx) = [(GIPx - GIDs)/GIPx] where GIPx = annual growth index for P. xylostella and
# GIDs = annual growth index for D. semiclausum. White = locations not suitable for P. xylostella (GIPx = 0). Yellow = locations where
# conditions are equally suitable for P. xylostella and D. semiclausum (GIPx = GIDs). Blue = locations where suitability for D. semiclausum is
# greater than that for P. xylostella (GIDs > GIPx). In all other locations (dark red through pink) conditions are more suitable for P. xylostella
# than for D. semiclausum (GIPx > GIDs): dark red = locations which favour P. xylostella over D. semiclausum most (RSPx > 0.45); 
# red (RSPx = 0.31–0.45); light red (RSPx = 0.16–0.30); pink (RSPx = 0.01–  0.15).


## Reclassify host; 

#GI (parasitoid) - GI (host)

# where new map = 0 yellow
# where new map is negative = blue
#

difference <- parasitoid - host

host2 <- reclassify (host, c(-Inf, 10, NA, 10, Inf, 0))
parasitoid2 <- reclassify (parasitoid, c(-Inf, 10, NA, 10, Inf, 0))

difference <- host2+ difference

outmap <- reclassify (difference, c(-100, -0.01, 0, -0.009, 0.001, 1, 0, 5, 2, 5, 10, 3, 10, 40, 4))

outmap <- ratify (outmap)

rat <- levels(outmap)[[1]]
rat$class <- c('Negative', 'None', "Positive1", "Positive2", "Positive3")
levels(outmap) <- rat

levelplot (outmap, col.regions=c("blue", "yellow", "pink", "red", "dark red"))


