library(sp)  #for example dataset
library(gstat) #for spatial interpolation
library(tidyverse) #just for some piping

## For a similar vignette, check out:
## https://cran.r-project.org/web/packages/gstat/vignettes/gstat.pdf


### Read in example data and make spatial objects
# Point data of field measurements
data(meuse)
coordinates(meuse) = ~x+y
# Read in grid over which interpolation will happen
data(meuse.grid)
gridded(meuse.grid) = ~x+y

### Visualize the data available
spplot(meuse, "zinc") 
hist(meuse$zinc)
spplot(meuse.grid) 

#######
### 0. Inverse-distance weighting
#######
idw_interp <- idw(log(zinc)~1, meuse, meuse.grid)
spplot(idw_interp["var1.pred"], main = "IDW predictions") 

#######
### 1. Fitting the variogram
#######
# Calculate an empirical variogram
vgm_emp <- variogram(log(zinc)~1, meuse)
plot(vgm_emp)
# Fit variogram model to empirical variogram, with or without 
# initial values of parameters
vgm_obs <- fit.variogram(vgm_emp, vgm("Sph"))
plot(variogramLine(vgm_obs,maxdist = max(vgm_emp$dist)), 
     type = "l")
points(vgm_emp[,c(2,3)])


#######
### 2. Kriging
#######

### ordinary kriging: stationary but unknown mean
ord_krig_interp <- krige(log(zinc)~1, meuse, meuse.grid, model = vgm_obs)
spplot(ord_krig_interp["var1.pred"], main = "ordinary kriging predictions")
spplot(ord_krig_interp["var1.var"],  main = "ordinary kriging variance")

### simple kriging: stationary with KNOWN mean
# only difference in code is specifying the beta parameter
simp_krig_interp <- krige(log(zinc)~1, meuse, meuse.grid, model = vgm_obs, beta = 5.9)
spplot(simp_krig_interp["var1.pred"], main = "simple kriging predictions")
spplot(simp_krig_interp["var1.var"],  main = "simple kriging variance")

### universal kriging: kriging residuals from a linear model
# a.k.a external drift kriging
# e.g. kriging PPT while also including PPT ~ elevation relationship
# The variogram should be based on residuals instead of observations
# residual variogram:
vgm_res <- variogram(log(zinc)~dist, meuse) %>%
  fit.variogram(vgm("Sph"))
# and then just include covariates in formula;
univ_krig_interp <- krige(log(zinc)~dist, meuse, meuse.grid, model = vgm_res)
spplot(univ_krig_interp["var1.pred"], main = "universal kriging predictions")


### Compare different interpolated fields
all_interp <- idw_interp["var1.pred"] %>%
  cbind(ord_krig_interp["var1.pred"]) %>%
  cbind(simp_krig_interp["var1.pred"]) %>%
  cbind(univ_krig_interp["var1.pred"]) 

names(all_interp) <- c("Inverse_Distance_Weighting",
                       "Ordinary_Kriging",
                       "Simple_Kriging",
                       "Universal_Kriging")

all_interp %>%
  spplot()
