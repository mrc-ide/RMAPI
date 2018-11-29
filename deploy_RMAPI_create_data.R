### --------------------------------------------------- START --------------------------------------------------- ###

library(devtools)
directory="C:/Users/Keith/Documents/Work/Work_backup"	# Location of folder in which RMAPI file is located, e.g. "C:/Users/YourName/Documents/PairwiseDataProcessing"
setwd(directory)					# Set specified folder as working directory
load_all("RMAPI")					# Load RMAPI files

#### ---------------------------------- CREATING AND SIMULATING DATA WITH BARRIER ------------------------------------------

latlongdata_file="C:/Users/Keith/Documents/Work/Work_backup/RMAPI/data/old_coord_data.txt"			# File containing latitude/longitude (or other X/Y) data in 2 columns
barrier_data_file="C:/Users/Keith/Documents/Work/Work_backup/RMAPI/data/barrier_data01.txt"			# File containing barrier data in 3 columns

coord_data=read.table(coord_data_file)
barrier_data=read.table(barrier_data_file)

long=coord_data["V1"]$V1
lat=coord_data["V2"]$V2
xbarrier=barrier_data["V1"]$V1
ybarrier=barrier_data["V2"]$V2
vbarrier=barrier_data["V3"]$V3
rbarrier=0.25
distance_model=0  # 0=linear dependence of pairwise data on distance
                  # 1=?

# generate data
d <- sim_custom_barrier(long, lat, xbarrier, ybarrier, rbarrier, vbarrier, distance_model)

# --------------------------------- Create new RMAPI project ----------------------------------- #

p <- rmapi_project()
summary(p)

# ----------------------------------------- Bind data ------------------------------------------ #

p <- bind_data(p, d, check_delete_output = FALSE)
summary(p)

# ----------------------------------------- Set up map ----------------------------------------- #

rhex=0.5
p <- create_map(p, hex_size = rhex)


# -------------------------------------- Run simulations --------------------------------------- #

Nperms=1  						  # Number of permutations to run when establishing statistical significance
eccentricity=0.975					# Eccentricity of ellipses
flag_nullmap=0						# Indicator of whether to use "null" distance-based map or not (0=no, 1=yes)
dist_model=distance_model				# Indicator of distance model to use for null map ()
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
