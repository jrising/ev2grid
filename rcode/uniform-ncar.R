setwd("~/research/ev2grid/ev2grid/rcode")

library(dplyr)
library(ggplot2)

df <- read.csv("../results/bytime-xsoc.csv")
df$hour <- sapply(as.POSIXlt(df$datetime, format="%Y-%m-%dT%H:%M"), function(dt) dt$hour) - 1
df$today <- substring(df$datetime, 1, 10) == "2023-07-17"
df$plugged <- df$hour < 8 | df$hour > 16

gp <- ggplot(subset(df, today), aes(group=hour)) +
    geom_blank()

for (hh in 0:23) {
    gp <- gp +
        geom_density(data=subset(df, today & plugged & hour == hh & soc_1 > 0),
                     aes(y=soc_plugged, x=..density.. / 10),
                     position=position_nudge(x=hh),
                     fill=scales::hue_pal()(24)[hh + 1],
                     alpha=0.6, adjust=1)
}

gp
