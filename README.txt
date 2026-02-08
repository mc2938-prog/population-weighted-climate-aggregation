Project description: The goal of this project is to analyze and visualize population-weighted climate data at the country level.
The project processes raster climate data, aggregates it into annual averages, computes temperature anomalies, and creates visualizations (maps and line plots) for each year from 1901 to 2012.
It generates a video showing the evolution of temperature anomalies over time.

Section 1 (Initial setup)
_________

- Here, I am downloading the main packages, namely sf, terra, RColorBrewer, Matrix, tmaptools, rworldmap and parallel.
- I am not downloading certain packages here but later because some of them intervene with codes. For example, tidyr blocks "info" variable to be processed.

Section 2 (Loading the data)
_________

- Loading the data that is given, namely climate, elevation and population raster
- Shifting the data so that it is centered on Pacific.

Section 3 (Aggregation of raster data in time)
__________

- Turning the monthly data to annual data.
- It's important to note that performing the annual aggregation, s2, takes a long time. Therefore, I have ran the code once, and
then downloaded it to my disk and the project file for ease.
 
Section 4 (Obtaining aggregation weights in matrix form)
_________

- Recentering climate raster to population so that they match
- Defining population density by dividing it by 5 according to the data
- Extracting climate data for each country polgygon from the climate raster
- Computing population weights for aggregation
- Combining all data frames into one
- Creating a sparse matrix which is used to aggregate climate data by country, weighted by population.
- Combining climate matrix and the shapefile

Section 5 (Data Cleaning and Transformation)
__________

- Here, because of the info processing that I have made, the data that have NA values are entered as 0.
- I have checked which countries had NA values and as expected, none of them were big countries and most of them were island countries so
I have decided to exclude them.

Section 6 (Unit Conversion: Kelvin to Celsius)
_________
- I wanted to confirm the rationality behind each country's temperatures, so I converted the Kelvin data to Celcius.
- The graphs are also in Celcius.

Section 7 (Calculate Temperature Anomalies
_________
- Calculating the temperature anomalies by the 1931-1960 baseline
- Adding the anomalies as a seperate column for ease.

Section 8 (Calculate Global Average Anomalies)
_________
- Downloading tidy and dplyr
- Calculating global average anomaly for each year

Section 9 (Mapping and Line plot)
__________
- Downloading ggplot2, magick, av and cowplot for mapping and making the video.
- Defining the color brackets similar to the one in the video
- Looping through each year to create plots: Inside, I am plotting both the map and the lineplot.
- Saving each year as a frame in the "frames" path.
- Making a video using magick and av








