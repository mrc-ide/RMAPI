### --------------------------------------------------- START --------------------------------------------------- ###

library(devtools)
directory="C:/Users/Keith/Documents/Work/Work_backup"	# Location of folder in which RMAPI file is located, e.g. "C:/Users/YourName/Documents/PairwiseDataProcessing"
setwd(directory)					# Set specified folder as working directory
load_all("RMAPI")					# Load RMAPI files

#### -------------------------------------- LOADING EXISTING PAIRWISE DATA -------------------------------------- ###

# Specify and load data #

latlongdata_file="C:/Users/Keith/Documents/Work/Work_backup/RMAPI/data/old_coord_data.txt"		# File containing latitude/longitude (or other X/Y) data in 2 columns
pairwisedata_file="C:/Users/Keith/Documents/Work/Work_backup/RMAPI/data/old_v_data.txt"			# File containing pairwise data in 1 column

d <- load_data(latlongdata_file,pairwisedata_file)							# Data loaded from files and converted to

# --------------------------------- Create new RMAPI project ----------------------------------- #

p <- rmapi_project()
summary(p)

# ----------------------------------------- Bind data ------------------------------------------ #

p <- bind_data(p, d)
summary(p)

# ----------------------------------------- Set up map ----------------------------------------- #

hex_size=500						# Size of hexes for map, in same units as latitude/longitude data
p <- create_map(p, hex_size)

# -------------------------------------- Run simulations --------------------------------------- #

Nperms=1  						  # Number of permutations to run when establishing statistical significance
eccentricity=0.975					# Eccentricity of ellipses
flag_nullmap=0						# Indicator of whether to use "null" distance-based map or not (0=no, 1=yes)
dist_model=1						# Indicator of distance model to use for null map ()
p <- run_sims(p, Nperms, eccentricity, flag_nullmap, dist_model)

# ------------------------------------- Plot maps ---------------------------------------------- #

if(flag_nullmap == 0)
{
# Basic map
plot_map(p, "map_values1")
points(p$data$long, p$data$lat, pch=20, cex=0.8)
} else 
{
# Basic map with distance-related "background" map, compensated map and permutation results #
par(mfrow=c(2,2))
plot_map(p, "map_values1")
points(p$data$long, p$data$lat, pch=20, cex=0.8)
plot_map(p, "map_values2")
points(p$data$long, p$data$lat, pch=20, cex=0.8)
plot_map(p, "map_values3")
points(p$data$long, p$data$lat, pch=20, cex=0.8)
plot_map(p, "empirical_p")
points(p$data$long, p$data$lat, pch=20, cex=0.8)
}