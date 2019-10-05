# ====================================
#
# On the beaten Path - Salvi, Williamson, Draper
# Replication Material
#
# ====================================

# ----- Loading Packages -----
pacman::p_load(ggplot2, reshape2, tidyverse, readr, ggmap, dplyr, raster, maptools, spatstat, sp, 
               cshapes, rgdal, geosphere, sf, rgeos, scales, sp, rgeos, dplyr, mwa )

# --------------------

#dev.off()

# ========================================
# Setting Up Observed and Simulated Data
# ========================================

## Important: set your working directory here

getwd()
#setwd()

set.seed(112)

# ---------- DRC shape area ----------

worldmap <- cshp(as.Date("1997-1-1"), useGW = FALSE)
worldmap@proj4string <- CRS("+proj=longlat +ellps=WGS84 ")
drc_map <- worldmap[worldmap$CNTRY_NAME == "Congo, DRC",]

# Check map
plot(drc_map)

# ---------- Create the DRC roads buffer ----------
roads <- readOGR("cod_trs_roads_osm.shp")
roads@proj4string <- CRS("+proj=longlat +ellps=WGS84 ")

# take only highways, primary and secondary roads, clip and buffer
roads_buff <- subset(roads, ntlclass %in% c("motorway", "primary", "primary_link",
                                            "secondary", "secondary_link")) %>%
  gSimplify(tol = 0.01) %>%
  gIntersection(drc_map, byid = TRUE) %>% 
  spTransform(CRS("+proj=utm +zone=32")) %>% # need this for accurate buffer width
  gBuffer(width = 5000) %>%
  spTransform(CRS("+proj=longlat +ellps=WGS84 ")) %>%
  gIntersection(drc_map, byid = TRUE)

# Visualize
plot(drc_map)
plot(roads_buff, add = TRUE, col = "lightgrey")


# ---------- Load observed events ----------

acled_drc <- read.csv("acled_drc.csv")
colnames(acled_drc)
conflict_years <- (1998:2000)
acled_drc <- acled_drc %>% filter(YEAR %in% conflict_years)
sp_point <- cbind(acled_drc$LONGITUDE, acled_drc$LATITUDE) 
colnames(sp_point) <- c("LONG","LAT") 

drc_events <- SpatialPointsDataFrame(coords = sp_point, data = acled_drc, 
                                     proj4string = CRS("+proj=longlat +ellps=WGS84 "))

drc_battles <- subset(drc_events, EVENT_TYPE %in% c("Battle-Government regains territory",
                                                    "Battle-No change of territory",
                                                    "Battle-Non-state actor overtakes territory"))

drc_vac <- subset(drc_events, EVENT_TYPE == "Violence against civilians")

# Visualize observed events
plot(drc_map)
points(drc_battles, col = alpha("blue", 0.4), pch = 1) # Battles
points(drc_vac, col = alpha("red", 0.4), pch = 2) # VAC
plot(roads_buff, add = TRUE, col = "burlywood1")
legend("bottomleft", legend=c("Battles", "VAC"),
       col=c("blue", "red"), pch = c(1, 2),
       cex = 0.9, bty = "n")

# ---------- Take events within buffer ---------- #
road_battles <- raster::intersect(drc_battles, roads_buff)
road_vac <- raster::intersect(drc_vac, roads_buff)

# Percentage of events falling within road buffer
nrow(road_battles@coords) / nrow(drc_battles@coords) # 69% of BATTLES w/i 5km roads
nrow(road_vac@coords) / nrow(drc_vac@coords) # 68% of VAC events fall w/i 5km of roads

# Percentage of country territory covered by 5km road buffer
area(roads_buff) / area(drc_map) # 13.8%

# Visualize
plot(drc_map)
plot(roads_buff, add = TRUE, col = "lightgrey")
points(drc_battles, col = alpha("grey", 0.4), pch = 2)
points(drc_vac, col = alpha("grey", 0.4), pch = 1)
points(road_battles, col = alpha("blue", 0.5), pch = 17)
points(road_vac, col = alpha("red", 0.5), pch = 19)
legend("bottomleft", legend=c("Violence against civilians (in buffer)", "Battle event (in buffer)", 
                              "Violence against civilians (outside buffer)",
                              "Battle event (outside buffer)"),
       col=c("red", "blue", "grey", "grey"), pch = c(19, 17, 1, 2),
       cex = 0.9, bty = "n")

# ---------- Simulate events ----------

# PPP function
to_ppp <- function(events,polygon){
  events  <- as(events,"SpatialPoints")
  events  <- as(events,"ppp")
  events  <- unique(events)
  polygon <- as(polygon,"SpatialPolygons")
  polygon <- as(polygon,"owin")
  ppp  <- ppp(events$x, events$y,window=polygon)
  return(ppp)
}

# Generate PPP from data
drc_ppp <- to_ppp(drc_battles, drc_map) # all drc
road_ppp <- to_ppp(drc_events, roads_buff) # near roads

# Looping a randomized version of point pattern (Road Buffer)
sim.road.df <- data.frame()
for (i in 1:5) {
  sims.road <- rpoispp(intensity(road_ppp), win=Window(road_ppp), ex = road_ppp)
  sims.road <- data.frame(sims.road) 
  sim.road.df <- rbind(sim.road.df, sims.road)
}

# dates of simulated events are sampled from the real ones
#sim.road.df$date <- sample(seq(as.Date('1997/01/01'), as.Date('2018/07/20'), by="day"), nrow(sim.road.df))
sim.road.df$date <-sample(as.Date(road_battles$EVENT_DATE), nrow(sim.road.df), replace = TRUE)
sim.road.df <- sim.road.df %>% arrange(date)


# Comparing # obs to sim obs
nrow(road_battles@coords)
nrow(sim.road.df)

# Looping a randomized version of point pattern (All DRC)
sim.drc.df <- data.frame()
for (i in 1:5) {
  sims.drc <- rpoispp(intensity(drc_ppp), win=Window(drc_ppp), ex = drc_ppp)
  sims.drc <- data.frame(sims.drc) 
  sim.drc.df <- rbind(sim.drc.df, sims.drc)
}

sim.drc.df$date <-sample(as.Date(drc_battles$EVENT_DATE), nrow(sim.drc.df), replace = TRUE)
sim.drc.df <- sim.drc.df %>% arrange(date)



# Comparing obs to sim
nrow(drc_battles@coords)
nrow(sim.drc.df)

# Naming columns
names(sim.road.df) <- c("LONG", "LAT", "TIMESTAMP")
names(sim.drc.df) <- c("LONG", "LAT", "TIMESTAMP")

# Transforming to sp
proj <- CRS("+proj=longlat +datum=WGS84") 

sim.road.sp <- SpatialPointsDataFrame(coords=sim.road.df[,1:2], data=sim.road.df, proj4string=proj) 
sim.drc.sp <- SpatialPointsDataFrame(coords=sim.drc.df[,1:2], data=sim.drc.df, proj4string=proj) 
plot(sim.drc.sp)

# Visualize observed AND simulated battles
plot(drc_map)
#points(sim.drc.sp, col = alpha("green", 0.4), pch = 20) # Sim DRC battles
points(sim.road.sp, col = alpha("blue", 0.4), pch = 20) # Sim road Battles
#points(drc_battles, col = alpha("black", 0.4), pch = 20) # Obs DRC battles
points(road_battles, col = alpha("red", 0.4), pch = 20) # Obs road battles

# ===========================================================================
# Descriptive: Computing Avg Distances (between observed/sim - battles/VAC)
# ===========================================================================

# ----- Coverage sensitivity check ------- # 

# iteratively change buffer width
widths <- seq(1000, 10000, by = 1000)
battle_capture <- vector()
area_capture <- vector()

for (i in 1:length(widths)) {
  
  # create buffer based on iterative width spec
  roads_buff <- roads_sub %>% 
    gBuffer(width = widths[i]) %>%
    spTransform(CRS("+proj=longlat +ellps=WGS84 ")) %>%
    gIntersection(drc_map, byid = TRUE)
  
  # find intersection of events with buffer area
  road_battles <- gIntersection(roads_buff, drc_battles, byid = TRUE)
  
  # percentage of battles falling within road buffer:
  battle_capture[i] <- nrow(road_battles@coords) / nrow(drc_battles@coords) 
  # percentage of country territory covered by buffer
  area_capture[i] <- area(roads_buff) / area(drc_map) 
  
}

buffer_percentages <- data.frame(battle_capture, area_capture, width = 1:10)


ggplot() + 
  geom_line(data = buffer_percentages, aes(width, battle_capture*100), col = "darkgrey") + 
  geom_line(data = buffer_percentages, aes(width, area_capture*100), col = "black") + 
  geom_text(aes(c(5.5, 5.5), c(80,25), 
                label = c("Percentage of battles \n in road buffer", 
                          "Percentage of country \n area in road buffer"))) + 
  labs(x = "Width of buffer (km)", y = "%", 
       title = "Coverage sensitivity of road buffer") + 
  coord_cartesian(ylim = c(0, 100)) + 
  scale_x_continuous(breaks = 1:10)


# ---------- Observed events ----------

# Calculating distance
d.road <- distm(road_battles, drc_vac)
d.drc  <- distm(drc_battles, drc_vac)

# Checking dimensions
dim(d.road)
dim(d.drc)

# Empty vectors to hold avg. distances
mean_dist_obs_road <- vector()
mean_dist_obs_drc <- vector()

# Observed Battles to VAC (5 nearest in roads buffer)
for (i in 1:nrow(d.road)) {
  distances <- sort(c(d.road[i,]))
  nearest_neigh <- distances[1:5]
  mean_dist_obs_road[i] <- mean(nearest_neigh)
}

# Observed Battles to VAC (5 nearest all DRC)
for (i in 1:nrow(d.drc)) {
  distances <- sort(c(d.drc[i,]))
  nearest_neigh <- distances[1:5]
  mean_dist_obs_drc[i] <- mean(nearest_neigh)
}

summary(mean_dist_obs_road)
summary(mean_dist_obs_drc)

# ---------- Simulated events ----------

# Calculating distance
d.sim.road <- distm(sim.road.sp, drc_vac)
d.sim.drc <- distm(sim.drc.sp, drc_vac)

# Checking dimensions
dim(d.sim.road)
dim(d.sim.drc)

# Empty vectors to hold avg. distances
mean_dist_sim_road <- vector()
mean_dist_sim_drc <- vector()

# Simulated Battles to VAC (5 nearest in roads buffer)
for (i in 1:nrow(d.sim.road)) {
  distances <- sort(c(d.sim.road[i,]))
  nearest_neigh <- distances[1:5]
  mean_dist_sim_road[i] <- mean(nearest_neigh)
}

# Simulated Battles to VAC (5 nearest all DRC)
for (i in 1:nrow(d.sim.drc)) {
  distances <- sort(c(d.sim.drc[i,]))
  nearest_neigh <- distances[1:5]
  mean_dist_sim_drc[i] <- mean(nearest_neigh)
}

par(mfrow=c(1,2))
summary(mean_dist_sim_drc[,1])
summary(mean_dist_sim_road[,1])

# ---------- T-tests / Creating Dataframe ----------

# T-test between road buffer and all DRC (obs vs. simulated)
t.test(as.vector(mean_dist_obs_road), as.vector(mean_dist_sim_road))
t.test(as.vector(mean_dist_obs_drc), as.vector(mean_dist_sim_drc))

# Create dataframes, add columns
mean_dist_obs_road <- as.data.frame(mean_dist_obs_road)
mean_dist_sim_road <- as.data.frame(mean_dist_sim_road)
mean_dist_obs_drc <- as.data.frame(mean_dist_obs_drc)
mean_dist_sim_drc <- as.data.frame(mean_dist_sim_drc)

mean_dist_obs_road$cat <- "Observed battles (road)"
mean_dist_sim_road$cat <- "Simulated battles (road)"
mean_dist_obs_drc$cat <- "Observed battles (DRC)"
mean_dist_sim_drc$cat <- "Simulated battles (DRC)"

mean_dist_obs_road$area <- "Road Buffer"
mean_dist_sim_road$area <- "Road Buffer"
mean_dist_obs_drc$area <- "DRC"
mean_dist_sim_drc$area <- "DRC"

colnames(mean_dist_obs_road) <- c("avg_dist", "cat", "area")
colnames(mean_dist_sim_road) <- c("avg_dist", "cat", "area")
colnames(mean_dist_obs_drc) <- c("avg_dist", "cat", "area")
colnames(mean_dist_sim_drc) <- c("avg_dist", "cat", "area")

# Combine into one (for ggplot)
combined.df <-rbind(mean_dist_obs_road,mean_dist_sim_road,mean_dist_obs_drc,mean_dist_sim_drc)
head(combined.df)

# ---------- Boxplot ----------

# Theme settings
main.theme <- theme_bw() + 
  theme(legend.key = element_rect(fill = NA, color = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(family = "Roboto", size = 14))
theme_set(main.theme)

# Changing order of cateogires
combined.df$cat <- ordered(combined.df$cat, levels = c("Observed battles (road)", "Simulated battles (road)","Observed battles (DRC)","Simulated battles (DRC)"))

# ggplot: distance of battles to five nearest VAC events, comparison between observed vs. simulated battles
ggplot() + 
  geom_boxplot(data = combined.df, aes(x=cat, avg_dist/1000)) + 
  labs(x = "", y = "Avg. Distance (in km, logged)", title = "Distance of battles to five nearest VAC events", subtitle = "Comparison between observed vs. simulated battles") +
  geom_text(aes(1.5,200, label = "T-test, p < 0.001")) + 
  geom_text(aes(3.5,200, label = "T-test, p < 0.001"))

ggplot() + 
  geom_boxplot(data = combined.df, aes(x=cat, avg_dist/1000, fill = area, color = area)) + 
  labs(x = "", y = "Avg. Distance (in km, logged)", title = "Distance of battles to five nearest VAC events", subtitle = "Comparison between observed vs. simulated battles") +
  geom_text(aes(1.5,200, label = "T-test, p < 0.001")) + 
  geom_text(aes(3.5,200, label = "T-test, p < 0.001")) + stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) + guides(color=FALSE, fill = FALSE)

# ===========================================================================
# Setup of the Analysis
# ===========================================================================


# MWA

#Load population data
population <- raster("gpw_v4.asc", 
                     proj4string = CRS("+proj=longlat +ellps=WGS84"))

#Build subsets of events and the country polygon
population_drc <-  crop(population, drc_map)

#Combine all events types in one df
# getting rid of some vars we don't need

colnames(sim.road.df)

road_sim_events <- sim.road.df %>%
  dplyr::rename(LATITUDE = LAT, LONGITUDE = LONG, EVENT_DATE = TIMESTAMP) %>%
  dplyr::mutate(type = "control") %>% dplyr::select(LONGITUDE, LATITUDE, EVENT_DATE, type)

road_battle_events <- road_battles@data %>% 
  dplyr::mutate(type = "treatment") %>%
  dplyr::select(LATITUDE, LONGITUDE, EVENT_DATE, type)


vac_events <- drc_vac@data %>%
  dplyr::mutate(type = "dependent") %>%
  dplyr::select(LATITUDE, LONGITUDE, EVENT_DATE, type)

names(road_sim_events)
names(road_battle_events)
names(vac_events)

mwa_events <- rbind(road_sim_events, road_battle_events, vac_events)

# Spatial coordinate for simulated events

sp_point2 <- cbind(road_sim_events$LONGITUDE, road_sim_events$LATITUDE) 
colnames(sp_point2) <- c("LONG","LAT") 

sp_point3 <- cbind(road_battle_events$LONGITUDE, road_battle_events$LATITUDE) 
colnames(sp_point3) <- c("LONG","LAT") 

sp_point4 <- cbind(vac_events$LONGITUDE, vac_events$LATITUDE) 
colnames(sp_point4) <- c("LONG","LAT") 

mwa_events_1 <- SpatialPointsDataFrame(coords = sp_point2, data = road_sim_events, 
                                       proj4string = CRS("+proj=longlat +ellps=WGS84 "))
mwa_events_2 <- SpatialPointsDataFrame(coords = sp_point3, data = road_battle_events, 
                                       proj4string = CRS("+proj=longlat +ellps=WGS84 "))
mwa_events_3 <- SpatialPointsDataFrame(coords = sp_point4, data = vac_events, 
                                       proj4string = CRS("+proj=longlat +ellps=WGS84 "))


# Adding Covariates

## Covariate 1: Distance to capital 
kinshasa <- as.data.frame(rbind(c(15.307045, -4.322447)))
colnames(kinshasa) <- c("LONG","LAT") 
kinshasa <- SpatialPointsDataFrame(coords = kinshasa, data = kinshasa,
                                   proj4string = CRS("+proj=longlat +ellps=WGS84 "))

mwa_events_1$capdist <- distm(mwa_events_1, kinshasa)
mwa_events_2$capdist <- distm(mwa_events_2, kinshasa)
mwa_events_3$capdist <- distm(mwa_events_3, kinshasa)

## Covariate 4: Terrain
# Source: https://topotools.cr.usgs.gov/GMTED_viewer/viewer.htm
#     Mean 30 arc sec
terrain1 <- raster("10S000E_20101117_gmted_mea300.tif")
terrain1@crs <- CRS("+proj=longlat +ellps=WGS84 ")
terrain2 <- raster("10S030E_20101117_gmted_mea300.tif")
terrain2@crs <- CRS("+proj=longlat +ellps=WGS84 ")
terrain3 <- raster("30S000E_20101117_gmted_mea300.tif")
terrain3@crs <- CRS("+proj=longlat +ellps=WGS84 ")
terrain4 <- raster("30S030E_20101117_gmted_mea300.tif")
terrain4@crs <- CRS("+proj=longlat +ellps=WGS84 ")

terrain <- raster::merge(terrain1, terrain2, terrain3, terrain4, extent = drc_map)
terrain_drc <- crop(terrain, drc_map)

# find elevation at each point 

mwa_events_1$terrain <- extract(terrain_drc, mwa_events_1)
mwa_events_2$terrain <- extract(terrain_drc, mwa_events_2)
mwa_events_3$terrain <- extract(terrain_drc, mwa_events_3)

## Covariate 5: Ethnic composition
geoepr <- readOGR("GeoEPR.shp")
drc_map@proj4string <- CRS("+proj=longlat +ellps=WGS84")

# crop to border
drc_eth <- gIntersection(geoepr, drc_map, byid = TRUE, drop_lower_td = TRUE)
drc_eth@proj4string <- CRS("+proj=longlat +ellps=WGS84")



# Count how many eth group polygons each point intersects with
pointsInPolygons <- sp::over(x = mwa_events_1, y = drc_eth, returnList = TRUE)
counting <- lapply(pointsInPolygons, FUN = function(x) length(x))
mwa_events_1$num_eth_grp <- t(do.call("cbind", counting))

pointsInPolygons <- sp::over(x = mwa_events_2, y = drc_eth, returnList = TRUE)
counting <- lapply(pointsInPolygons, FUN = function(x) length(x))
mwa_events_2$num_eth_grp <- t(do.call("cbind", counting))

pointsInPolygons <- sp::over(x = mwa_events_3, y = drc_eth, returnList = TRUE)
counting <- lapply(pointsInPolygons, FUN = function(x) length(x))
mwa_events_3$num_eth_grp <- t(do.call("cbind", counting))


# Creating an empty dataset for MWA

dataset <- as.data.frame(cbind(rep(as.character("NA"),nrow(mwa_events))))
names(dataset) <- "type"
dataset$lon <- 0.0
dataset$lat <- 0.0
dataset$timestamp <- as.Date("1900-01-01")
dataset$population <- 0.0

dataset$num_eth_grp <- 0
dataset$terrain <- 0.0
dataset$capdist <- 0.0

## Copy data into the DF

dataset$type <- c(rep("control",length(mwa_events_1[,1])),rep("treatment",length(mwa_events_2[,1])),rep("dependent",length(mwa_events_3[,1])))

dataset$lat  <- c(mwa_events_1$LATITUDE,mwa_events_2$LATITUDE,mwa_events_3$LATITUDE)
dataset$lon  <- c(mwa_events_1$LONGITUDE,mwa_events_2$LONGITUDE,mwa_events_3$LONGITUDE)
dataset$timestamp  <-c(mwa_events_1$EVENT_DATE,mwa_events_2$EVENT_DATE,mwa_events_3$EVENT_DATE)
dataset$timestamp <- as.Date(dataset$timestamp)

# add covariates
dataset$population <- c(population_drc[mwa_events_1,],population_drc[mwa_events_2,],population_drc[mwa_events_3,])

dataset$num_eth_grp  <- c(mwa_events_1$num_eth_grp,mwa_events_2$num_eth_grp,mwa_events_3$num_eth_grp)

dataset$terrain <- c(mwa_events_1$terrain,mwa_events_2$terrain,mwa_events_3$terrain)

dataset$capdist <- c(mwa_events_1$capdist,mwa_events_2$capdist,mwa_events_3$capdist)

##MWA Analysis
# Specify required parameters:
t_window <- c(5,50,5)
spat_window <- c(5,60,5)
t_unit <- "days" 
TCM <- TRUE
weighted <- FALSE

treatment <- c("type","treatment")
# - column and entries that indicate control events 
control  <- c("type","control")
# - column and entries that indicate dependent events 
dependent <- c("type","dependent")
# - columns for matching
matchColumns <- c("population", "num_eth_grp", "terrain", "capdist")

# Execute method (replicates the result presented on the paper):

dataset1<- dataset[complete.cases(dataset), ] 

library(rJava) # needs to be loaded separatly on some machines
options(java.parameters = "-Xmx1g")


results <- matchedwake(dataset1, t_window, spat_window, treatment, control, dependent, matchColumns, weighted = weighted, t_unit = t_unit, TCM = TCM)

plot(results)
summary(results)



