# ---------------------------------------------------------------
# Description: Go from raster climate to population-weighted climate
# at the country-level
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# 1. Directories and packages
# ---------------------------------------------------------------

# Clean up workspace and load or install necessary packages if necessary
rm(list=ls())
want <- c("sf","terra","RColorBrewer","Matrix","tmaptools","rworldmap","parallel")
need <- want[!(want %in% installed.packages()[,"Package"])]
if (length(need)) install.packages(need)
lapply(want, function(i) require(i, character.only=TRUE))
rm(want, need)

# Working directories
dir <- list()
dir$root <- dirname(getwd())
dir$weather   <- paste(dir$root,"/hw2/data/weather",sep="")
dir$pop <- paste(dir$root,"/hw2/data/population",sep="")


# ---------------------------------------------------------------
# 2. Load data
# ---------------------------------------------------------------

# Raster - elevation
e <- rast(list.files(dir$weather, full.names=T)[1])
id <- land <- e
land[] <- ifelse(e[]>0,1,NA)
id[] <- ifelse(e[]>0,1:ncell(e),NA)
plot(e, main="Elevation")
plot(land, main="Land indicator")
plot(id, main="ID = original cell position")
# Note that climate data appears centered in Pacific

# Raster - climate
s <- rast(list.files(dir$weather, full.names=T)[2])

# Raster population
fname <- list.files(dir$pop, full.names=T, recursive = T)
fname <- fname[grepl("[.]tif",fname)] # only select the file name with ".tif" in the name
p <- rast(fname)
# Note that the population data is centered in the Atlantic

# Shapefiles
co <- st_as_sf(countriesLow)

# The difference in raster projection (center in pacific or atlantic)
# will cause issues when you try to combin raster with polygons
# here is an illustration:
plot(land)
plot(vect(co), add=T)
# polygons are only overlaid on the Eastern hemisphere

# Solution? Recenter cliamte data so that it is centered in the Atlantic

rland <- terra::rotate(land)
rid   <- terra::rotate(id)

# Corresponds to the position BEFORE the raster was "recentered"
newid <- rid
newid[] <- ifelse(rland[]==1,1:ncell(e),NA)

# Check it out
matchid <- data.frame(original.id=c(rid[]), new.id=c(newid[]))
matchid <- matchid[!is.na(matchid[,1]),] # drop rows with NAs
head(matchid)
# How to understand this?
# The first row indicates that the cell that was located in the
# ORIGINAL projection (centered in Pacific) in position 9285 is
# now located in position 8925 in the NEW projection (centered)
# in the Atlantic.
# Note this is simply a shift of all IDs
table(apply(matchid, 1, diff)) # the difference is 360 for half of the cells and -360 fo rthe other half

# ---------------------------------------------------------------
# 3. Aggregation of raster data in time (monthly to annual)
# ---------------------------------------------------------------

# Test run the code for a smaller set of years
if (F) {
  t <- s[[1:(12*10)]]
  years <- 1901:(1901+9)#1901:2012
  index <- rep(years, each=12)
  system.time({ # 36s for 5 years, 73s for 10 years, etc... so plan accordingly
    t2 <- tapp(t, index=index, fun=mean)
  })
}

# Create year index for full data
years <- 1901:2012
index <- rep(years, each=12)

mcores <- detectCores()-1

# Perform annual aggregation - I have saved s2 in a disk to save time. You can just run s2 <- readRDS("s2_annual_climate.rds")
system.time({
  s2 <- tapp(s, index, fun=mean, cores=mcores)
})

saveRDS(s2, file = "s2_annual_climate.rds")
s2 <- readRDS("s2_annual_climate.rds")

# ---------------------------------------------------------------
# 4. Aggregation of raster data in space
# ---------------------------------------------------------------

# 4.1 Obtain aggregation weights in matrix form -----

# Recenter the climate raster as extent of s2 and p do not match
s2_recentered <- terra::rotate(s2)

# Resample population raster to match climate raster resolution
p_resampled <- resample(p, s2_recentered, method="bilinear")

# The rest of this part is mostly, if not fully, taken from lecture notes.
# Adjust population density by dividing the resampled population raster by 5 and then match to climate raster resolution again.
dens <- p_resampled/5
dens2 <- resample(dens, s2_recentered, method="bilinear")

# Extract climate data for each country polygon from the recentered climate raster and assign a unique ID to each country
info <- extract(s2_recentered, vect(co), weights=TRUE, cells=TRUE, exact=TRUE)
co$ID <- 1:nrow(co)

# Add country info
co.info <- st_drop_geometry(co[, c("ID", "ISO3", "ADMIN")])
info <- merge(info, co.info)

# Filter out rows where either the ISO3 code or the cell index is missing
info <- info[!is.na(info$ISO3) & !is.na(info$cell), ]

# Compute population weights for aggregation
info2 <- lapply(unique(info$ISO3), function(coname) {
  df <- info[info$ISO3 == coname, ]
  x <- dens2[df$cell]
  x <- x / sum(x, na.rm=T)
  df$popweight <- unlist(x)
  df <- df[!is.na(df$popweight), c("ID", "cell", "popweight")]
  return(df)
})

# Assign country ISO3 codes as names to the list elements
names(info2) <- as.character(unique(info$ISO3))

# Remove countries with no valid data (i.e., empty data frames)
info2 <- info2[-c(which(sapply(info2, nrow)==0))]

# Combine all country data frames into a single data frame
temp <- do.call("rbind", info2)

# Create a sparse matrix for population weights
P <- sparseMatrix(i=temp$ID, j=temp$cell, x=temp$popweight, dims=c(nrow(co), ncell(s2_recentered)))

# Assign country ISO3 codes as row names
rownames(P) <- co$ISO3 

# 4.2 Perform aggregation ----

# Perform matrix multiplication to aggregate the climate data (M) using the population weights (P)
M <- s2[]

# The result (Mco) is a matrix where each row corresponds to a country and each column corresponds to a year
Mco <- P %*% M

# Convert the aggregated climate matrix (Mco) into a data frame
d <- as.data.frame(as.matrix(Mco))

# Add a column for country ISO3 codes using the row names of Mco
d$ISO3 <- rownames(d)

# 4.3 Combine it with a map/shapefile ----

co <- merge(co, d, by="ISO3")

print(co)

# ------------------------------------------------------------------------------
# 5. Data Cleaning and Transformation
# ------------------------------------------------------------------------------

#Delete the rows that have 0 values aka NA, which are small countries that likely lack data like Aruba, Aland..

# Drop the geometry column temporarily
co_df <- st_drop_geometry(co)

# Identify columns that contain temperature data
temp_cols <- grep("^X\\d{4}$", names(co_df), value = TRUE)

# Filter out rows where all temperature values are zero
filter_indices <- rowSums(co_df[temp_cols] != 0, na.rm = TRUE) > 0
co_filtered_df <- co_df[filter_indices, ]

# Reattach the geometry column correctly
co_filtered <- st_as_sf(co_filtered_df, geometry = st_geometry(co)[filter_indices])

# Check the result
print(co_filtered)

# ------------------------------------------------------------------------------
# 6. Unit Conversion: Kelvin to Celsius
# ------------------------------------------------------------------------------
# Convert temperature values from Kelvin to Celsius
co_filtered_df <- st_drop_geometry(co_filtered)
temp_cols <- grep("^X\\d{4}$", names(co_filtered_df), value = TRUE)
co_filtered_df[temp_cols] <- co_filtered_df[temp_cols] - 273.15

# Reattach the geometry column to the updated data frame
co_filtered <- st_as_sf(co_filtered_df, geometry = st_geometry(co_filtered))

# Create a subset of the data with only country names and temperature columns
country_temps <- co_filtered[, c("ADMIN", temp_cols)]

# ------------------------------------------------------------------------------
# 7. Calculate Temperature Anomalies
# ------------------------------------------------------------------------------

# Define baseline columns (1931–1960)
baseline_cols <- grep("^X19(3[1-9]|4[0-9]|5[0-9]|60)$", names(co_filtered), value = TRUE)

# Calculate the baseline temperature for each country (mean of baseline years)
co_filtered$baseline <- rowMeans(st_drop_geometry(co_filtered)[baseline_cols], na.rm = TRUE)

# Compute temperature anomalies by subtracting the baseline from annual temperatures
anomalies <- st_drop_geometry(co_filtered)[temp_cols] - co_filtered$baseline

# Rename anomaly columns to include "_anomaly" suffix
colnames(anomalies) <- paste0(colnames(anomalies), "_anomaly")

# Combine the anomalies with the main data frame
co_filtered <- cbind(co_filtered, anomalies)

# ------------------------------------------------------------------------------
# 8. Calculate Global Average Anomalies
# ------------------------------------------------------------------------------

# Load tidyr now because otherwise it interferes with previous data
install.packages("tidyr")
install.packages("dplyr")
library(tidyr)
library(dplyr)

# Calculate global average anomaly for each year
global_anomalies <- co_filtered %>%
  st_drop_geometry() %>%
  summarise(across(ends_with("_anomaly"), ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "year", values_to = "global_anomaly") %>%
  mutate(year = as.numeric(gsub("X|_anomaly", "", year)))

# ------------------------------------------------------------------------------
# 9. Visualization: Mapping and Line Plot
# ------------------------------------------------------------------------------

##install required packages for mapping

install.packages("ggplot2")
install.packages("magick")
install.packages("av")
install.packages("cowplot")

library(ggplot2)
library(magick)
library(av)
library(cowplot)

# Define the breaks and corresponding colors
color_breaks <- c(-Inf, -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, Inf)
color_palette <- c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0",
                   "#f0f0f0", "#FFC0CB", "#FA8072", "#FF0000", "#b2182b", "#67001f")

world <- st_as_sf(getMap(resolution = "low"))


# Loop through each year to create the plots
for (year in 1901:2012) {
  # Filter the data for the current year
  year_col <- paste0("X", year, "_anomaly")
  co_filtered$current_anomaly <- co_filtered[[year_col]]
  
  # Create the map plot
  map_plot <- ggplot() +
    geom_sf(data = co_filtered, aes(fill = current_anomaly), color = NA) +
    geom_sf(data = world, fill = NA, color = "black", linewidth = 0.2) + # Add country borders
    scale_fill_stepsn(
      colors = color_palette,
      breaks = color_breaks,
      limits = c(-3,3),
      name = "Temperature anomaly (°C)",
      guide = guide_colorsteps(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(18, "cm") # Adjust legend width
      )
    ) +
    theme_minimal() +
    labs(title = paste0(year), # Bold year as title
         fill = NULL) + # Remove legend title
    theme(
      legend.position = "bottom",
      legend.key.width = unit(5, "cm"), # Make legend 5 times bigger
      legend.text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold") # Center and bold title
    )
  
  # Create the line plot for global anomalies
  line_data <- global_anomalies %>% filter(year <= !!year)
  line_plot <- ggplot(line_data, aes(x = year, y = global_anomaly)) +
    geom_line(color = "black", linewidth = 1.2) + # Black and bold line
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Dashed line at y = 0
    geom_point(data = line_data %>% filter(year == !!year), # Add point only for the current year
               color = "black", size = 2) + 
    scale_x_continuous(breaks = seq(1900, 2020, 20), limits = c(1900, 2020)) + # Set x-axis breaks to 20-year increments
    theme_minimal() +
    labs(x = NULL, y = "Global temperature anomaly (°C)") + # Update y-axis title
    theme(
      axis.title.y = element_text(size = 12),
      panel.grid.major = element_blank(), # Remove grid lines
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black"), # Add x-axis line
      axis.line.y = element_line(color = "black") # Add y-axis line
    )
  
  # Combine the two plots into one image
  combined_plot <- cowplot::plot_grid(map_plot, line_plot, ncol = 1, rel_heights = c(2, 1))
  
  # Save the combined plot as a PNG file with a white background
  ggsave(paste0("frames/plot_", year, ".png"), combined_plot, width = 10, height = 8, bg = "white")
}

frame_files <- list.files(path = "frames", pattern = "plot_.*\\.png", full.names = TRUE)
frame_files <- sort(frame_files)

# Create a temporary directory to store resized frames
temp_dir <- tempdir()
resized_files <- file.path(temp_dir, basename(frame_files))

# Resize frames to 800x600 (or smaller if needed)
for (i in seq_along(frame_files)) {
  img <- image_read(frame_files[i])
  img <- image_scale(img, "800x600") # Resize to 800x600 pixels
  image_write(img, resized_files[i])
  rm(img) # Free up memory
  gc() # Trigger garbage collection
}

# Create the video using the resized frames
av::av_encode_video(
  input = resized_files, # Use resized frames
  output = "output.mp4", # Output video file name
  framerate = 5 # Lo7wer frame rate (5 fps) to reduce load
)
