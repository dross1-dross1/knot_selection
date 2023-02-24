# Setup -----------------------------------------------------------------------------------------------------------

cat("\014")  # ctrl+L
#if (!is.null(dev.list()["RStudioGD"])) { dev.off(dev.list()["RStudioGD"]) }  # remove all plots
try(dev.off(dev.list()["RStudioGD"]), silent=T)
rm(list = ls())  # remove all variables
set.seed(100)

# Imports
library(dplyr)
library(plotly)
library(ggplot2)
library(zoo)
library(ggforce)

# Set WD to directory of current script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Data Loading ----------------------------------------------------------------------------------------------------

df.ozone = read.csv("8hour_44201_2021_ozone.csv")
# df.co = read.csv("8hour_42101_2021_co.csv")

# filter only mainland cases
df.ozone = filter(df.ozone, Longitude > -130)
df.ozone = filter(df.ozone, Latitude > 25)

# convert and combine date and time columns
df.ozone$dt.local = as.POSIXct(paste(df.ozone$Date.Local, df.ozone$Time.Local), format="%Y-%m-%d %H:%M")
df.ozone$dt.gmt = as.POSIXct(paste(df.ozone$Date.GMT, df.ozone$Time.GMT), format="%Y-%m-%d %H:%M")

df.ozone = subset(df.ozone, select=-c(Date.Local, Time.Local, Date.GMT, Time.GMT))

# allow for matching of lat/lon and x/y coordinates
normalize = function(x) { (x-min(x))/(max(x)-min(x)) }
df.ozone$x = df.ozone$Longitude %>% normalize
df.ozone$y = df.ozone$Latitude %>% normalize

# need year seperate
df.ozone$year = df.ozone$dt.local %>% format("%Y") %>% as.numeric

# Functions -------------------------------------------------------------------------------------------------------

sample_df = function(df, sample_size=nrow(df)) { df[sample(1:nrow(df), sample_size), ] }

graph_point = function(df, sample_size=nrow(df)) {
  df = df %>% distinct(Latitude, Longitude, .keep_all=T)
  sample.df = sample_df(df)
  ggplot(sample.df) + geom_point(aes(x=Longitude, y=Latitude))
}

graph_point_site = function(df, id, sample_size=nrow(df)) {
  df = df %>% distinct(Latitude, Longitude, .keep_all=T)
  sample.df = sample_df(df)
  lat = filter(df, site.id == id)$Latitude[1]
  print(lat)
  lon = filter(df, site.id == id)$Longitude[1]
  ggplot(sample.df) +
    geom_point(aes(x=Longitude, y=Latitude)) +
    geom_circle(aes(x0=lon, y0=lat, r=.5), color="red", size=1)
}

graph_point_value = function(df, sample_size=nrow(df)) {
  df = df %>% distinct(Latitude, Longitude, .keep_all=T)
  sample.df = sample_df(df)
  ggplot(sample.df) + geom_point(aes(x=Longitude, y=Latitude, color=Mean.Excluding.Concurred.Flags))
}

# turns a selection of the data into a format where it can be fed into a knot selection algorithm
df_dataify = function(df.data.in) {
  data.frame(
    id=1:nrow(df.data.in),
    real_id=df.data.in$site.id,
    x=df.data.in$x,
    y=df.data.in$y,
    lat=df.data.in$Latitude,
    lon=df.data.in$Longitude,
    signal=df.data.in$Mean.Excluding.Concurred.Flags,
    year=df.data.in$year
  )
}

# Code ------------------------------------------------------------------------------------------------------------

# check for na (after adding the datetime stuff)
nrow(df.ozone[complete.cases(df.ozone), ]) / nrow(df.ozone)  # are there any na values

# collection map
graph_point(df.ozone)  # graph all collection locations

# regarding each site
df.ozone = df.ozone %>% group_by(Latitude, Longitude) %>% mutate(site.id=cur_group_id()) %>% as.data.frame  # create id for each site
df.ozone.sites = df.ozone %>% count(site.id)  # how many times each site shows up in data
df.ozone.sites = df.ozone.sites[order(-df.ozone.sites$n), ]  # sort descending site occurrences

# 553 = most samples from this site
graph_point_site(df.ozone, 553)

# plot site over time with a moving average
df.ozone.sample = filter(df.ozone, site.id == 553)  # df.ozone subset with only collections from this site
df.ozone.sample = df.ozone.sample %>% mutate(ravg_value=rollapplyr(Mean.Excluding.Concurred.Flags, 100, mean, partial=T))  # make rolling average and assign it to "ravg_value"
ggplot(df.ozone.sample) + geom_line(aes(x=dt.gmt, y=ravg_value)) + geom_smooth(aes(x=dt.gmt, y=ravg_value), method="gam") # plot "ravg_value" over time

# snapshot in time
df.ozone.snap = filter(df.ozone, dt.gmt == df.ozone$dt.local[1])
df.ozone.snap$site.idx = 1:nrow(df.ozone.snap)
graph_point_value(df.ozone.snap)

# for simultaneous_knot_selection
# df.data = subset(df.ozone.snap, select=c(id, x, y, health, year))
# df.data = data.frame(
#   id=df.ozone.snap$site.idx,
#   x=df.ozone.snap$Longitude %>% normalize,
#   y=df.ozone.snap$Latitude %>% normalize,
#   health=df.ozone.snap$Mean.Excluding.Concurred.Flags,
#   year=1950,
#   case_cntl=df.ozone.snap$Mean.Excluding.Concurred.Flags > mean(df.ozone.snap$Mean.Excluding.Concurred.Flags)
# )

# 12 Snapshots ----------------------------------------------------------------------------------------------------

# times for all snapshots put into one dataframe
snapshots = data.frame(times=c(
  "2021-01-01 07:00:00",
  "2021-02-01 07:00:00",
  "2021-03-01 07:00:00",
  "2021-04-01 07:00:00",
  "2021-05-01 07:00:00",
  "2021-06-01 07:00:00",
  "2021-07-01 07:00:00",
  "2021-08-01 07:00:00",
  "2021-09-01 07:00:00",
  "2021-10-01 07:00:00",
  "2021-11-01 07:00:00",
  "2021-12-01 07:00:00"
))
snapshots$times = as.POSIXct(snapshots$times, format="%Y-%m-%d %H:%M:%S")

# still creating snapshot data
dfo.snapshots = filter(df.ozone, dt.local %in% snapshots$times)
dfo.snapshots$site.idx = 1:nrow(dfo.snapshots)
dfo.snapshots$month = dfo.snapshots$dt.local %>% format("%m") %>% as.numeric
dfo.snapshots %>% head

# creating both df with all snapshots (df.data.all), and a list with the 12 snapshots split up into a list (df.data.12list)
df.data.all = df_dataify(dfo.snapshots)
df.data.12list = vector(mode="list", length=12)
for (i in 1:nrow(snapshots)) { df.data.12list[[i]] = df_dataify(dfo.snapshots %>% filter(month == i)) }

# current month
ccnt = 1

# prepping for knot selection
df.data = df.data.12list[[ccnt]]
df.algo = data.frame(df.data)

# graphing single month
ggplot(df.data) + geom_point(aes(x=x, y=y, color=signal), size=3)

# Temp for model fitting ------------------------------------------------------------------------------------------

# df.sks.12list = vector(mode="list", length=12)
# df.kav4.12list = vector(mode="list", length=12)

# ggplot() +
#   geom_point(data=df.data.12list[[2]], aes(x=x, y=y, color=signal)) +
#   geom_circle(data=df.sks.12list[[11]], aes(x0=x, y0=y, r=.01), color="red")# +
  # geom_circle(data=df.kav4.12list[[ccnt]], aes(x0=x, y0=y, r=.05), color="blue")

# Saving Data -----------------------------------------------------------------------------------------------------

# saveRDS(df.data.all, file="df_data_all.RDS")
# saveRDS(df.data.12list, file="df_data_12list.RDS")
# saveRDS(df.sks.12list, file="df_sks_12list.RDS")
# saveRDS(df.kav4.12list, file="df_kav4_12list.RDS")

# saveRDS(df.sks.12list, file="df_sks_12list_NEW.RDS")
