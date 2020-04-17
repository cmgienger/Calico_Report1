#required packages
library(adehabitatHR) #for home range calculations
library(data.table) #manipulate S3 and S4 data tables
library(ggplot2) #for graphic output
library(ggfortify) #to allow ggplot2 to read spatial data
library(grid) #to add annotations to the output
library(lubridate) #working with date-time formats
library(OpenStreetMap) #for obtaining raster images
library(pbapply) #needed for progress bar
library(plotly) #for interactive xy plot
library(rgdal) #for converting spatial data
library(scales) #determining breaks and labels for axes and legends
library(sp) #for converting spatial data

#load data with fake date and time
data <- read.csv("test_data.csv")

data$DATE <- mdy(data$DATE) #use lubridate to specify incoming date format (mdy); 
#the current format is "mdy", it then converts column to "ymd" format
str(data) #look at DATE format to verify that it is converted to "ymd"

#make a chart showing the tracking duration for each individual
q <- ggplot(data, aes(DATE, LIZARDNUMBER)) +
  geom_count (color = "blue", alpha=0.5) + #set symbol size proportional to # of overlapping obs (same day)
  scale_x_date(date_breaks = "12 month", labels = date_format("%Y"))+
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 14), 
        legend.key = element_rect(fill = NA)) +
  labs(title = "Tracking Duration and Intensity") +
  labs(size = "Number of Obs.")
q

#interactive examination of data points for outliers with Plotly
p <- ggplot() + geom_point(data=data, aes(EASTING,NORTHING, color=LIZARDNUMBER)) +
  labs(x="Easting", y="Northing")
ggplotly(p)

#split into multiple files for individuals
lapply(split(data, data$LIZARDNUMBER), 
       function(x)write.csv(x, file = paste(x$LIZARDNUMBER[1],".csv"), row.names = FALSE))

#create list of individual files created in previous step #edit pattern based on individual ids
files <- list.files(path = ".", pattern = "[MF]+[0-9]", full.names = TRUE)

#raster from openstreetmap +-0.001 added buffer from min max easting and northing
#should be tuned to match the perimeter of the analyses if necessary
utm_points <- cbind(data$EASTING, data$NORTHING)
utm_locations <- SpatialPoints(utm_points, proj4string=CRS("+proj=utm +zone=12 +datum=WGS84"))
latlon_locations <- as.data.frame(spTransform(utm_locations, CRS("+proj=longlat +datum=WGS84")))
colnames(latlon_locations) <- c("x","y")
raster <- openmap(c(max(latlon_locations$y)+0.002, min(latlon_locations$x)-0.002), 
                  c(min(latlon_locations$y)-0.002, max(latlon_locations$x)+0.002), type = "bing")
raster_utm <- openproj(raster, projection = "+proj=utm +zone=12 +ellps=WGS84 +units=m +no_defs")

#plot location on raster map #add + facet_wrap(~LIZARDNUMBER) separate by individual
autoplot(raster_utm, expand = TRUE) + theme_bw() + theme(legend.position="none") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_point(data=data, aes(EASTING,NORTHING, color=LIZARDNUMBER), size = 3, alpha = 0.8) +
  theme(axis.title = element_text(face="bold")) + labs(x="Easting", y="Northing")

#creating spatial data frame for all points
x <- as.data.frame(data$EASTING)
y <- as.data.frame(data$NORTHING)
xy <- c(x,y)
data.proj <- SpatialPointsDataFrame(xy,data, proj4string = CRS("+proj=utm +zone=12 +ellps=WGS84 +units=m +no_defs"))

#creating homerange in adehabitat for all points
xy <- SpatialPoints(data.proj@coords)
mcp.out <- mcp(xy, percent=100, unout="ha")
#plot(xy)
#plot(mcp.out)
mcp_area <- as.data.frame(mcp.out@data$area)
colnames(mcp_area) <- "Hectares"
write.table(mcp_area, "MCP_Area.csv", sep = ",", row.names = TRUE)

#KDE creation for all points
kde<-kernelUD(xy, h="href", kern="bivnorm", grid=1000)
ver <- getverticeshr(kde, 95)
ver$area
kde@h$h
#plot(data.proj@coords)
#plot(ver)
kde_area <- as.data.frame(ver$area)
colnames(kde_area) <- "Hectares"
write.table(kde_area, "KDE_Area.csv", sep = ",", row.names = TRUE)

#mcp plot for all points
mcp.points <- cbind((data.frame(xy)),data$LIZARDNUMBER)
colnames(mcp.points) <- c("x","y", "lizardnumber")
mcp.poly <- fortify(mcp.out, region = "id")
mcp.plot <- ggplot()+
  geom_polygon(data=mcp.poly, aes(x=mcp.poly$long, y=mcp.poly$lat))+
  geom_point(data=mcp.points, aes(x=x, y=y,color = lizardnumber))
mcp.plot

#kde plot for all points
kde.points <- cbind((data.frame(data.proj@coords)),data$LIZARDNUMBER)
colnames(kde.points) <- c("x","y","lizardnumber")
kde.poly <- fortify(ver, region = "id")
kde.plot <- ggplot()+
  geom_polygon(data=kde.poly, aes(x=kde.poly$long, y=kde.poly$lat))+
  geom_point(data=kde.points, aes(x=x, y=y, color = lizardnumber))
kde.plot

#MCP plots and raster mapping for all points
autoplot(raster_utm, expand = TRUE) + theme_bw() + theme(legend.position="none") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_polygon(data=mcp.poly, x=mcp.poly$long, y=mcp.poly$lat, alpha = 0.5) +
  geom_point(data=mcp.points, aes(x=x, y=y, color = lizardnumber)) +
  theme(axis.title = element_text(face="bold")) + labs(x="Easting", y="Northing")

#KDE plots and raster mapping for all points
autoplot(raster_utm, expand = TRUE) + theme_bw() + theme(legend.position="none") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_polygon(data=kde.poly, x=kde.poly$long, y=kde.poly$lat, alpha = 0.5) +
  geom_point(data=kde.points, aes(x=x, y=y, color = lizardnumber)) +
  theme(axis.title = element_text(face="bold")) + labs(x="Easting", y="Northing")

#looping function for mcp
mcp_analysis <- function(filename){
  data <- read.csv(file = filename)
  x <- as.data.frame(data$EASTING)
  y <- as.data.frame(data$NORTHING)
  xy <- c(x,y)
  data.proj <- SpatialPointsDataFrame(xy,data, proj4string = CRS("+proj=utm +zone=12 +ellps=WGS84 +units=m +no_defs"))
  xy <- SpatialPoints(data.proj@coords)
  mcp.out <- mcp(xy, percent=100, unout="ha")
  area <- as.data.frame(round(mcp.out@data$area,4))
  .rowNamesDF(area, make.names=TRUE) <- data$LIZARDNUMBER
  write.table(area,file="MCP_Hectares.csv", 
              append=TRUE,sep=",", col.names=FALSE, row.names=TRUE)
  mcp.points <- cbind((data.frame(xy)),data$LIZARDNUMBER)
  colnames(mcp.points) <- c("x","y", "lizardnumber")
  mcp.poly <- fortify(mcp.out, region = "id")
  units <- grid.text(paste(round(mcp.out@data$area,2)," ha"), x=0.9,  y=0.95,
                     gp=gpar(fontface=4, cex=0.9), draw = FALSE)
  mcp.plot <- ggplot() +
    geom_polygon(data=mcp.poly, aes(x=mcp.poly$long, y=mcp.poly$lat), alpha=0.5) +
    geom_point(data=mcp.points, aes(x=x, y=y)) + theme_bw() + 
    labs(x="Easting (m)", y="Northing (m)", title=mcp.points$lizardnumber) +
    theme(legend.position="none", plot.title = element_text(face = "bold", hjust = 0.5)) + 
    annotation_custom(units)
  mcp.plot
}

#individual mcp run
mcp_analysis("./F24 .csv")

#run all individuals for mcp, takes 12sec
#lapply(files,mcp_analysis)
#pblapply(files, mcp_analysis) #runs with progressbar

#looping function for kde
kde_analysis <- function(filename){
  data <- read.csv(file = filename)
  x <- as.data.frame(data$EASTING)
  y <- as.data.frame(data$NORTHING)
  xy <- c(x,y)
  data.proj <- SpatialPointsDataFrame(xy,data, proj4string = CRS("+proj=utm +zone=12 +ellps=WGS84 +units=m +no_defs"))
  xy <- SpatialPoints(data.proj@coords)
  kde<-kernelUD(xy, h="href", kern="bivnorm", grid=1000)
  ver <- getverticeshr(kde, 95)
  area <- as.data.frame(round(ver$area,4))
  .rowNamesDF(area, make.names=TRUE) <- data$LIZARDNUMBER
  write.table(area,file="KDE_Hectares.csv", 
              append=TRUE,sep=",", col.names=FALSE, row.names=TRUE)
  kde.points <- cbind((data.frame(data.proj@coords)),data$LIZARDNUMBER)
  colnames(kde.points) <- c("x","y","lizardnumber")
  kde.poly <- fortify(ver, region = "id")
  units <- grid.text(paste(round(ver$area,2)," ha"), x=0.9,  y=0.95,
                     gp=gpar(fontface=4, cex=0.9), draw = FALSE)
  kde.plot <- ggplot() +
    geom_polygon(data=kde.poly, aes(x=kde.poly$long, y=kde.poly$lat), alpha = 0.5) +
    geom_point(data=kde.points, aes(x=x, y=y)) + theme_bw() + 
    labs(x="Easting (m)", y="Northing (m)", title=kde.points$lizardnumber) +
    theme(legend.position="none", plot.title = element_text(face = "bold", hjust = 0.5)) + 
    annotation_custom(units)
  kde.plot
}

#individual kde run
kde_analysis("./F24 .csv")

#run all individuals for kde, takes 3min
#lapply(files,kde_analysis)
#pblapply(files, kde_analysis)  #runs with progressbar

#looping function for mcp with raster
mcp_raster <- function(filename){
  data <- read.csv(file = filename)
  x <- as.data.frame(data$EASTING)
  y <- as.data.frame(data$NORTHING)
  xy <- c(x,y)
  data.proj <- SpatialPointsDataFrame(xy,data, proj4string = CRS("+proj=utm +zone=12 +ellps=WGS84 +units=m +no_defs"))
  xy <- SpatialPoints(data.proj@coords)
  mcp.out <- mcp(xy, percent=100, unout="ha")
  mcp.points <- cbind((data.frame(xy)),data$LIZARDNUMBER)
  colnames(mcp.points) <- c("x","y", "lizardnumber")
  mcp.poly <- fortify(mcp.out, region = "id")
  units <- grid.text(paste(round(mcp.out@data$area,2)," ha"), x=0.9,  y=0.95,
                     gp=gpar(fontface=4, col="white", cex=0.9), draw = FALSE)
  mcp.plot <- autoplot(raster_utm, expand = TRUE) + theme_bw() + theme(legend.position="none") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    geom_polygon(data=mcp.poly, aes(x=mcp.poly$long, y=mcp.poly$lat), alpha=0.8) +
    geom_point(data=mcp.points, aes(x=x, y=y)) + 
    labs(x="Easting (m)", y="Northing (m)", title=mcp.points$lizardnumber) +
    theme(legend.position="none", plot.title = element_text(face = "bold", hjust = 0.5)) + 
    annotation_custom(units)
  mcp.plot
}

#individual mcp run with raster
mcp_raster("./F24 .csv")

#run all individuals for mcp with raster, takes 1min
#lapply(files,mcp_raster)
#pblapply(files, mcp_raster) #runs with progressbar

kde_raster <- function(filename){
  data <- read.csv(file = filename)
  x <- as.data.frame(data$EASTING)
  y <- as.data.frame(data$NORTHING)
  xy <- c(x,y)
  data.proj <- SpatialPointsDataFrame(xy,data, proj4string = CRS("+proj=utm +zone=12 +ellps=WGS84 +units=m +no_defs"))
  xy <- SpatialPoints(data.proj@coords)
  kde<-kernelUD(xy, h="href", kern="bivnorm", grid=1000)
  ver <- getverticeshr(kde, 95)
  kde.points <- cbind((data.frame(data.proj@coords)),data$LIZARDNUMBER)
  colnames(kde.points) <- c("x","y","lizardnumber")
  kde.poly <- fortify(ver, region = "id")
  units <- grid.text(paste(round(ver$area,2)," ha"), x=0.9,  y=0.95,
                     gp=gpar(fontface=4, col="white", cex=0.9), draw = FALSE)
  kde.plot <- autoplot(raster_utm, expand = TRUE) + theme_bw() + theme(legend.position="none") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    geom_polygon(data=kde.poly, aes(x=kde.poly$long, y=kde.poly$lat), alpha = 0.8) +
    geom_point(data=kde.points, aes(x=x, y=y)) +
    labs(x="Easting (m)", y="Northing (m)", title=kde.points$lizardnumber) +
    theme(legend.position="none", plot.title = element_text(face = "bold", hjust = 0.5)) + 
    annotation_custom(units)
  kde.plot
}

#individual kde run
kde_raster("./F24 .csv")

#run all individuals for kde, takes 3.5min
#lapply(files,kde_raster)
#pblapply(files, kde_raster) #runs with progressbar

#trajectory analysis
traj_analysis <- function(filename){
  relocs_data <- read.csv(file = filename)
  relocs <- as.ltraj(cbind(relocs_data$EASTING, relocs_data$NORTHING),id=relocs_data$LIZARDNUMBER, typeII = FALSE, date=NULL)
  relocs.df <- ld(relocs)
  relocs_dist <- as.data.frame(sum(sapply(relocs.df$dist, sum, na.rm=TRUE)))
  colnames(relocs_dist) <- "Total Distance"
  name <- relocs.df$id[1]
  row.names(relocs_dist) <- name
  relocs_units <- grid.text(paste(round(relocs_dist,2),"m"), x=0.9, y=0.9, 
                            gp=gpar(fontface=3, col="black", cex=0.9), draw = FALSE)
  reloc.plot <- ggplot() + theme_classic() + geom_path(data=relocs.df, aes(x=x,y=y), linetype = "dashed", colour = "red",
    arrow = arrow(length=unit(.5,"cm"), angle = 20, ends="last", type = "closed")) +
    geom_point(data=relocs.df, aes(x=x, y=y)) + geom_point(data=relocs.df, aes(x=x[1], 
    y=y[1]), size = 3, color = "darkgreen", pch=0) +
    labs(x="Easting (m)", y="Northing (m)", title=relocs.df$id[1]) +
    theme(legend.position="none", plot.title = element_text(face = "bold", hjust = 0.5)) +
    annotation_custom(relocs_units)
  reloc.plot
}

#pblapply(files, traj_analysis)
traj_analysis("F24 .csv")

#distance over time analysis
dist_analysis <- function(filename){
  relocs_data <- read.csv(file = filename)
  relocs <- as.ltraj(cbind(relocs_data$EASTING, relocs_data$NORTHING),id=relocs_data$LIZARDNUMBER, typeII = FALSE, date=NULL)
  relocs.df <- ld(relocs)
  relocs_dist <- as.data.frame(sum(sapply(relocs.df$dist, sum, na.rm=TRUE)))
  colnames(relocs_dist) <- "Total Distance"
  name <- relocs.df$id[1]
  row.names(relocs_dist) <- name
  write.table(relocs_dist,file="reloc_dist.csv", 
              append=TRUE,sep=",", col.names=FALSE, row.names=TRUE)
  dist.plot <- ggplot() + geom_point(data = relocs.df, aes(x=date, y=dist), na.rm = TRUE) + theme_classic() +
    geom_segment(data = relocs.df, aes(x=date, xend=date, y=0, yend=dist,), stat = "identity", na.rm = TRUE) +
    labs(x="Date", y="Distance (m)", title=relocs.df$id[1]) + 
    theme(legend.position="none", plot.title = element_text(face = "bold", hjust = 0.5))
  dist.plot
}

#pblapply(files, dist_analysis)
dist_analysis("F24 .csv")



#remove individual files after analysis, this works so be careful!
#unlink(files, recursive = FALSE)