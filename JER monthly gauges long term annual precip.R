##############
## Spatial interpolation of Jornada precipitation
## R Work Group 11/30
## Written by Darren James and Heather Savoy
## If you want the data to run this script, ask Darren
##############



library(tidyverse)
library(raster)
library(rgdal)
library(akima)
library(gstat)

#### Keep everything in UTM Zone 13N WGS84
# "+proj=utm +zone=13 +datum=WGS84 +units=m +no_defs"

# Import boundary shapefiles
JER.boundary <- readOGR(dsn="data", layer="jer_boundary") 

CRDDC.boundary <- readOGR(dsn="data", layer="cdrrc_boundary") 

projection(JER.boundary) == projection(CRDDC.boundary)

# Import DEM
DEM.import <- raster("data/JER_CDRRC_DEM.tif") 

projection(DEM.import) == projection(JER.boundary)

plot(DEM.import)
lines(JER.boundary)
lines(CRDDC.boundary)

# Create hillshade from DEM.import
slope <- terrain(DEM.import, opt="slope", unit='radians')
aspect <- terrain(DEM.import, opt="aspect", unit='radians')
hillshade <- hillShade(slope, aspect, angle=55, direction=300)

plot(hillshade, col = gray.colors(30, start = 0, end = 1), legend = FALSE)
lines(JER.boundary)
lines(CRDDC.boundary)

# Import long-term annual means
long.term.annual.means <- read_csv("data/LongTerm annual means from monthly gauges.csv") %>%
  mutate(annual_precip = format(round(precip, digits=1), nsmall = 1))

range(long.term.annual.means$precip)

# Create shapefile of the long-term annual means
long.term.annual.means.shp <- SpatialPointsDataFrame(coords = long.term.annual.means[,c("Xcoord", "Ycoord")],
                                                     proj4string = CRS("+proj=utm +zone=13 +datum=WGS84 +units=m +no_defs +ellps=WGS84"),
                                                     data = long.term.annual.means)


# Plot of long-term annual means overlaid on the hillshade
# jpeg(filename = "figures/Long_term_AnnualPrecip.jpg", width = 13, height = 13, units = "in", res = 600)
old.par <- par(mar = c(0,0,0,0))
plot(hillshade, col = gray.colors(30, start = 0, end = 1), legend = FALSE, axes = FALSE, box = FALSE)
plot(JER.boundary, add = TRUE, lwd = 1)
plot(CRDDC.boundary, add = TRUE, lwd = 1)
raster::text(x = long.term.annual.means.shp$Xcoord, y = long.term.annual.means.shp$Ycoord, labels=long.term.annual.means.shp$annual_precip, cex=1.4, font=1, col = "darkblue")
scalebar(6437.38, xy=c(336000, 3590500), type='bar', divs=6, label = c("0", "2 miles", "4 miles"))
# dev.off()
par(old.par)

# Set the extent of the interpolation area
view.xmin <- 312294
view.xmax <- 361500
view.ymin <- 3589459
view.ymax <- 3633904

# Bilinear interpolation of the long-term means
fld.linear <- with(long.term.annual.means, akima::interp(x = Xcoord, y = Ycoord, z = precip,
                                                  xo  = seq(view.xmin, view.xmax, 100), 
                                                  yo = seq(view.ymin, view.ymax, 100), 
                                                  linear = TRUE, 
                                                  extrap = TRUE)) 



# Create raster of the bilinear interpolation and plot it
# Notice the interpolation is limited to the minimum convex polygon of the gauge locations
raster.linear <- raster(fld.linear)

plot(raster.linear, col = blues9)
plot(JER.boundary, add = TRUE, lwd = 1)
plot(CRDDC.boundary, add = TRUE, lwd = 1)
raster::text(x = long.term.annual.means.shp$Xcoord, 
             y = long.term.annual.means.shp$Ycoord, 
             labels=long.term.annual.means.shp$annual_precip, 
             cex=1, font=1, col = "black")


# Bilinear spline interpolation of the long-term means
fld.spline <- with(long.term.annual.means, akima::interp(x = Xcoord, y = Ycoord, z = precip,
                                                  xo  = seq(view.xmin, view.xmax, 100), 
                                                  yo = seq(view.ymin, view.ymax, 100), 
                                                  linear = FALSE, 
                                                  extrap = TRUE)) 

range(fld.spline$z) # Yikes!! Some areas got a lot of rain?!

# Create raster of the bilinear interpolation and plot it
# Notice the interpolation is overwhelmed by artifacts near the Dona Ana Mountains
raster.spline <- raster(fld.spline)
plot(raster.spline, col = blues9)
plot(JER.boundary, add = TRUE, lwd = 1)
plot(CRDDC.boundary, add = TRUE, lwd = 1)
raster::text(x = long.term.annual.means.shp$Xcoord, 
             y = long.term.annual.means.shp$Ycoord, 
             labels=long.term.annual.means.shp$annual_precip, 
             cex=1, font=1, col = "black")


# Plot the bilinear spline with a different color scale
cuts=c(6,7,8,9,10,11,12, 13, 14, 21000) #set breaks

plot(raster.spline, breaks=cuts, col = blues9) #plot with defined breaks
plot(JER.boundary, add = TRUE, lwd = 1)
plot(CRDDC.boundary, add = TRUE, lwd = 1)
text(x = long.term.annual.means.shp$Xcoord, 
     y = long.term.annual.means.shp$Ycoord, 
     labels=long.term.annual.means.shp$annual_precip,  
     cex=1, font=1, col = "black")


#### Geostatistical interpolation (universal kriging) of the long-term means
# See the gstat_examples.R script for an intro to kriging
# First, add the DEM values to the long term spatial object
long.term.annual.means.shp$JER_CDRRC_DEM <- extract(DEM.import,
         long.term.annual.means.shp,
         method = "bilinear")

# Generate an empirical variogram of residuals
vgm_emp <- variogram(precip~JER_CDRRC_DEM, 
                     long.term.annual.means.shp,
                     alpha = 45, # you can modify the coordinate direction
                     tol.hor = 15) # and filter to pairs of points in the directions
plot(vgm_emp)
# Fit variogram model to empirical variogram
vgm_res <- fit.variogram(vgm_emp, 
                         vgm("Mat",
                             psill=0.2,
                             anis = c(45,0.5)), # heterogeneity can be anisotropic
                         fit.kappa = TRUE)

plot(vgm_emp,vgm_res)

# If you want to gauge model fitting error, you can extract
# the sum of square errors. Not so helpful on its own, but you can compare
# among candidate models
attr(vgm_res, "SSErr")

# Make a gridded spatial version of the DEM
DEM.grid <- DEM.import %>%
  crop(extent(c(view.xmin, view.xmax,view.ymin, view.ymax)))  %>% # to match earlier examples
  as('SpatialGridDataFrame')

# Now perform the universal kriging
# This is time consuming! Note that the entext matches earlier examples, but the
# resolution is much finer: 5655x4315 vs 445x493
system.time(univ_krig_interp <- krige(precip~JER_CDRRC_DEM,
                                      long.term.annual.means.shp,
                                      DEM.grid,
                                      model = vgm_res,
                                      nmax = 10)) # can save computational time
# elapsed: 676.966s, or ~11 min

# Create raster of the geostatistical interpolation and plot it
# Pattern is mostly elevation, but you can see influence of local minima as well,
# e.g. the lower 8.6 value surrounded by 9+ values in the northwest corner has lower
# local values.
raster.uk <- raster(univ_krig_interp["var1.pred"])
plot(raster.uk, col = blues9)
plot(JER.boundary, add = TRUE, lwd = 1)
plot(CRDDC.boundary, add = TRUE, lwd = 1)
raster::text(x = long.term.annual.means.shp$Xcoord, 
             y = long.term.annual.means.shp$Ycoord, 
             labels=long.term.annual.means.shp$annual_precip, 
             cex=1, font=1, col = "black")



