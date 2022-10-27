# this is R code put together by Tim Meehan to get eBird relative abundances
# for sharp-shinned hawk prey species in northern California during fall 
# migration


# setup ------------------------------------------------------------------------
# set up a bunch of stuff. load libraries, set options, define spatial 
# reference systems, etc.

# libraries
library(RColorBrewer)
library(colorspace)
library(forcats)
library(grid)
library(mapview)
library(doParallel)
library(sf)
library(ggplot2)
library(tidyr)
library(lubridate)
library(stringr)
library(dplyr)
library(errors)
library(ebirdst)
library(raster)

# options
options(scipen=9999999)
options(max.print=99999)
options(errors.notation = "plus-minus")

# pointers
extract <- raster::extract
select <- dplyr::select

# time series plot function 
theme_timeseries <- function (base_size = 11, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = rel(0.9), color="gray10", 
                                     angle = 0),
          axis.text.y = element_text(size = rel(0.9), color="gray10", 
                                     angle = 0),
          strip.background = element_rect(fill = "grey80"),
          legend.key = element_rect(fill = "white", colour = NA),
          plot.title = element_text(size=14, hjust = 0.5,
                                    margin=margin(t=5, b=10)),
          legend.position="right",
          complete = TRUE)
}; theme_set(theme_timeseries())

# working directory
setwd("~/GitHub/ggro_sharpie_prey/")

# ebird key
# set_ebirdst_access_key("6tsm3ctdenmv", overwrite=T)

focal_weeks <- 31:48

# crs's
ebird_crs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs "
epsg102008 <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
# ------------------------------------------------------------------------------



# get focal area map ----------=------------------------------------------------
# focal areas
focal <- read_sf("./data/focal_area.shp") %>%
  mutate(FID=0) %>% dplyr::select(FID) %>% st_transform(ebird_crs)

# plot it all for QAQC
mapview::mapview(st_geometry(st_transform(focal, epsg102008)), border="black", 
     lwd=2)
# ------------------------------------------------------------------------------



# get spp data -----------------------------------------------------------------
species_table <- read.csv(file = "./data/ebird_prey_list.csv") %>%
  rename(prey_spp=common_name)
spp_list <- species_table[,1]

spp_list <- c("American Goldfinch","American Robin","Anna's Hummingbird",
  "Band-tailed Pigeon","Blue-gray Gnatcatcher","Brewer's Blackbird",
  "Brown-headed Cowbird","Brown Creeper","Bushtit","California Scrub-Jay",     
  "California Thrasher","California Towhee","Cassin's Vireo","Cedar Waxwing",
  "Chestnut-backed Chickadee","Common Yellowthroat",      
  "Dark-eyed Junco","Eurasian Collared-Dove","Fox Sparrow",
  "Golden-crowned Kinglet","Grasshopper Sparrow","Green-tailed Towhee",      
  "Hermit Thrush","Hooded Oriole","Horned Lark","House Finch","Hutton's Vireo",
  "Lincoln's Sparrow", "MacGillivray's Warbler","Mountain Bluebird",
  "Mountain Chickadee",       
  "Mourning Dove","Nashville Warbler","Northern Mockingbird",     
  "Olive-sided Flycatcher","Orange-crowned Warbler","Pacific-slope Flycatcher", 
  "Pine Siskin","Purple Finch","Pygmy Nuthatch",         
  "Red-breasted Nuthatch","Red-winged Blackbird","Ruby-crowned Kinglet",     
  "Savannah Sparrow","Song Sparrow","Spotted Towhee",           
  "Steller's Jay","Swainson's Thrush","Townsend's Warbler",       
  "Tricolored Blackbird","Varied Thrush","Vaux's Swift",             
  "Warbling Vireo","Western Bluebird","Western Meadowlark",    
  "Western Tanager","White-breasted Nuthatch","White-crowned Sparrow",    
  "White-throated Sparrow","Wilson's Warbler","Wrentit",                  
  "Yellow-rumped Warbler","Yellow Warbler")

top_choice <- c("Hermit Thrush",  "Spotted Towhee",  "Fox Sparrow",  "Yellow Warbler",  
"Swainson's Thrush", "Red-breasted Nuthatch", "Dark-eyed Junco",
"Chestnut-backed Chickadee", "American Robin", "California Towhee")

bottom_choice <- c("Grasshopper Sparrow",  "Olive-sided Flycatcher",
"Nashville Warbler",  "Vaux's Swift", "White-throated Sparrow",
"Horned Lark", "Hooded Oriole", "Mountain Chickadee", "Green-tailed Towhee", 
"Mountain Bluebird")

spp_list <- c(top_choice, bottom_choice)
# make sure that there are ebird layers for these species
all(spp_list %in% ebirdst_runs$common_name==T)
# ------------------------------------------------------------------------------



# get sharpie occurrence -------------------------------------------------------
sharp_path <- ebirdst_download(species = "Sharp-shinned Hawk", force=F)
sharp_occ <- load_raster("occurrence", path = sharp_path)
sharp_occ <- crop(sharp_occ[[focal_weeks]], focal)
sharp_occ[is.na(sharp_occ)] <- 0
sharp_occ <- mask(sharp_occ, focal)
plot(sharp_occ)
tmpfilter <- sharp_occ[[9]] < 0.05
filtered_image <- mask(sharp_occ[[8]], tmpfilter, maskvalue=1)
filtered_image %>% mapview(na.color=NA, alpha.regions=0.7)
saveRDS(sharp_occ, "sharpie_occurrence_map.RData")
# ------------------------------------------------------------------------------



# this function does most of the work ------------------------------------------
# source this function and then run it in a parallel loop below for multiple
# species.
out_stk <- stack()
out_names <- c()
for(i in 1:length(spp_list)){
  spp1 <- spp_list[i]
  print(paste("starting species", spp1))
  # libs
  library(ebirdst)
  library(raster)
  library(sf)
  # get ebird status and trends abundance layers using the ebirdst package
  sp_path <- ebirdst_download(species = spp1, force=F)
  abd <- load_raster("abundance", path = sp_path)
  # trim everything for focal prey species
  abd <- crop(abd[[focal_weeks]], focal)
  abd[is.na(abd)] <- 0
  abd <- mask(abd, focal)
  abd_early <- mean(abd[[1:9]])
  abd_late <- mean(abd[[10:18]])
  earl_name <- gsub("'", "", tolower(paste0(gsub(" ", "_", spp1), "_abd_early")))
  late_name <- gsub("'", "", tolower(paste0(gsub(" ", "_", spp1), "_abd_late")))
  out_stk <- raster::addLayer(out_stk, abd_early, abd_late)
  out_names <- c(out_names, earl_name, late_name)
}

names(out_stk) <- out_names
mapview(out_stk[[20]], na.color=NA, alpha.regions=0.7)
prey_stk <- raster::projectRaster(out_stk, crs = 4326)
mapview(prey_stk[[20]], na.color=NA, alpha.regions=0.7)
names(prey_stk) <- gsub("[.]", "_", names(prey_stk))
saveRDS(prey_stk, "prey_maps.RData")
# ------------------------------------------------------------------------------





