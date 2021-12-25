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
set_ebirdst_access_key("4kur40litsen")

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

# make sure that there are ebird layers for these species
all(spp_list %in% ebirdst_runs$common_name==T)
# ------------------------------------------------------------------------------



# get sharpie occurrence -------------------------------------------------------
sharp_path <- ebirdst_download(species = "Sharp-shinned Hawk", force=T)
sharp_occ <- load_raster("occurrence", path = sharp_path)
sharp_occ <- crop(sharp_occ[[focal_weeks]], focal)
sharp_occ[is.na(sharp_occ)] <- 0
sharp_occ <- mask(sharp_occ, focal)
plot(sharp_occ)

# get stuff to estimate occurrence proportional error
sharp_abd <- crop(load_raster("abundance", path = sharp_path), focal)
sharp_abd[is.na(sharp_abd)] <- 0
sharp_abd <- mask(sharp_abd, focal)
sharp_lcl <- crop(load_raster("abundance_lower", path = sharp_path), focal)
sharp_lcl[is.na(sharp_lcl)] <- 0
sharp_lcl <- mask(sharp_lcl, focal)
sharp_ucl <- crop(load_raster("abundance_upper", path = sharp_path), focal)
sharp_ucl[is.na(sharp_ucl)] <- 0
sharp_ucl <- mask(sharp_ucl, focal)
plot(sharp_ucl)
# ------------------------------------------------------------------------------



# this function does most of the work ------------------------------------------
# source this function and then run it in a parallel loop below for multiple
# species.
run_species <- function(spp1="American Dipper"){
  print(paste("starting species", spp1))
  
  # libs
  library(ebirdst)
  library(raster)
  library(sf)
  library(errors)
  
  # get ebird status and trends abundance layers using the ebirdst package
  sp_path <- ebirdst_download(species = spp1, force=T)
  abd <- load_raster("abundance", path = sp_path)
  lcl <- load_raster("abundance_lower", path = sp_path)
  ucl <- load_raster("abundance_upper", path = sp_path)
  
  # trim everything for focal prey species
  abd <- crop(abd[[focal_weeks]], focal)
  abd[is.na(abd)] <- 0
  abd <- mask(abd, focal)
  # plot(abd)
  lcl <- crop(lcl[[focal_weeks]], focal)
  lcl[is.na(lcl)] <- 0
  lcl <- mask(lcl, focal)
  # plot(lcl)
  ucl <- crop(ucl[[focal_weeks]], focal)
  ucl[is.na(ucl)] <- 0
  ucl <- mask(ucl, focal)
  # plot(ucl)

  # make empty df and loop through weeks 
  out_df <- c()
  for(w in 1:nlayers(sharp_occ)){
    # get ids
    week_i <- names(sharp_occ)[w]
    
    # pull appropriate weekly maps for sharpies and prey
    sharp_i <- sharp_occ[[w]]
    sharp_abd_i <- sharp_abd[[w]]
    sharp_abd_i_se <- ((sharp_ucl[[w]] - sharp_lcl[[w]]) / 2.56)
    # plot(sharp_i)
    # plot(sharp_abd_i)
    # plot(sharp_abd_i_se)
    lay_i <- abd[[w]]
    # plot(lay_i)
    lay_i_se <- ((ucl[[w]] - lcl[[w]]) / 2.56)
    # plot(lay_i_se)
    
    # get values and assign errors for propagation for prey
    vals_i <- values(lay_i)
    vals_i <- vals_i[!is.na(vals_i)]
    # summary(vals_i)
    se_i <- values(lay_i_se)
    se_i <- se_i[!is.na(se_i)] 
    # summary(se_i)
    errors(vals_i) <- se_i
    
    # get values and assign errors for propagation for sharpies
    s_vals_i <- values(sharp_i)
    s_vals_i <- s_vals_i[!is.na(s_vals_i)]
    # summary(s_vals_i)
    s_abd_vals_i <- values(sharp_abd_i)
    s_abd_vals_i <- s_abd_vals_i[!is.na(s_abd_vals_i)]
    s_abd_vals_i_se <- values(sharp_abd_i_se)
    s_abd_vals_i_se <- s_abd_vals_i_se[!is.na(s_abd_vals_i_se)]
    prop_i <- s_abd_vals_i_se / s_abd_vals_i
    prop_i[is.na(prop_i)] <- 0
    # summary(prop_i)
    errors(s_vals_i) <- prop_i * s_vals_i
    
    # extract mean ai, oi and se
    sharp_wt_mean_prey_ai <- as.numeric(mean(vals_i*s_vals_i))
    sharp_wt_mean_prey_ai_se <- errors(mean(vals_i*s_vals_i))
    mean_prey_ai <- as.numeric(mean(vals_i))
    mean_prey_ai_se <- errors(mean(vals_i))
    mean_sharp_oi <- as.numeric(mean(s_vals_i))
    mean_sharp_oi_se <- errors(mean(s_vals_i))
    
    # build df
    out <- data.frame(week=week_i, prey_spp=spp1,
                      mean_sharp_oi=mean_sharp_oi,
                      mean_sharp_oi_se=mean_sharp_oi_se,
                      mean_prey_ai=mean_prey_ai,
                      mean_prey_ai_se=mean_prey_ai_se,
                      sharp_wt_mean_prey_ai=sharp_wt_mean_prey_ai,
                      sharp_wt_mean_prey_ai_se=sharp_wt_mean_prey_ai_se)
    
    # add row to the output data
    out_df <- rbind(out_df, out)
    print(paste("finished week", w, "of", nlayers(sharp_occ)))
  }
  write.csv(out_df, file=paste0("./output/Prey_species_", spp1, ".csv"), row.names=F)
  return(out_df)
}

# test
# test_df <- run_species()
# ------------------------------------------------------------------------------



# run the work horse function in parallel --------------------------------------
# make cluster
cl <- makeCluster(6)
registerDoParallel(cl)

# run function
loop_out1 <- foreach(i=1:length(spp_list), 
                     .packages=c("sf", "errors", "ebirdst", "raster")) %dopar% {
                                   run_species(spp_list[i])
                                 }
# stop cluster
stopCluster(cl)
# ------------------------------------------------------------------------------



# time series plots and tables -------------------------------------------------
# get data, clean and save
dat1 <- bind_rows(loop_out1) %>% mutate(week=ymd(gsub("w", "", week))) %>%
  left_join(species_table)
mean_abund <- dat1 %>% group_by(prey_spp) %>%
  summarise(seasonal_mean_prey_ai=mean(mean_prey_ai)) %>%
  ungroup()
dat1 <- dat1 %>% left_join(mean_abund) %>% 
  arrange(desc(seasonal_mean_prey_ai)) 
dat1$prey_spp <- factor(as.character(dat1$prey_spp))
dat1$detected_in_diet <- factor(as.character(dat1$detected_in_diet))
dat1$prey_spp <- fct_reorder(dat1$prey_spp, 
                             dat1$seasonal_mean_prey_ai, .desc=T)
write.csv(dat1, "./output/relative_predator_occurrence_prey_abundance_north_coast_CA_fall.csv")

# plotting data
dat2 <- dat1 %>% filter(detected_in_diet=="yes") 
dat2$lt <- rep(rep(1:12, each=18), length.out=nrow(dat2))

ggplot(dat2, 
       aes(x=week, y=mean_prey_ai, color=prey_spp,
           linetype=prey_spp)) + 
  geom_line() +
  theme(legend.text = element_text(size = 7)) + 
  scale_y_continuous(trans="sqrt") +
  labs(x="Week", 
       y="Mean prey abundance index (square root scale)",
       color="",
       linetype="") +
  scale_linetype_manual(values=1:nlevels(dat2$prey_spp))
ggsave("./output/time_series.pdf", height=6, width=12, dpi = 300)

# table
dat1 %>% group_by(prey_spp, detected_in_diet) %>% summarise(mean_ai=mean(mean_prey_ai)) %>%
  arrange(desc(mean_ai))
# ------------------------------------------------------------------------------

 

   
