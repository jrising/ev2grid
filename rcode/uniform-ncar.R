setwd("~/research/ev2grid/ev2grid/rcode")

library(dplyr)
library(ggplot2)

for (suffix in c('', '-sol2')) {
    df <- read.csv(paste0("../results/bytime-xsoc", suffix, ".csv"))
    df$hour <- sapply(as.POSIXlt(df$datetime, format="%Y-%m-%dT%H:%M"), function(dt) dt$hour) - 1
    df$today <- substring(df$datetime, 1, 10) == "2023-07-17"
    df$plugged <- df$hour < 8 | df$hour > 16

    gp <- ggplot(subset(df, today), aes(group=hour)) +
        coord_cartesian(xlim=c(1, 26)) +
        geom_blank()

    for (hh in 0:23) {
        gp <- gp +
            geom_density(data=subset(df, today & plugged & hour == hh & soc_1 >= 0),
                         aes(y=soc_plugged, x=..density.. / 10),
                         position=position_nudge(x=hh + 1),
                         fill=scales::hue_pal()(24)[hh + 1],
                         alpha=0.6, adjust=1)
    }

    gp + theme_bw() + scale_y_continuous("Plugged-in SOC (%)", labels=scales::percent) + xlab("Hour of the day")
    ggsave(paste0("../figures/uniform-ncar", suffix, ".pdf"), width=6.5, height=4)
}


df <- read.csv(paste0("../results/bytime-xsoc.csv"))
df$hour <- sapply(as.POSIXlt(df$datetime, format="%Y-%m-%dT%H:%M"), function(dt) dt$hour) - 1
df$today <- substring(df$datetime, 1, 10) == "2023-07-17"
df$plugged <- df$hour < 8 | df$hour > 16

df2 <- read.csv(paste0("../results/bytime-xsoc-sol2.csv"))
df2$hour <- sapply(as.POSIXlt(df2$datetime, format="%Y-%m-%dT%H:%M"), function(dt) dt$hour) - 1
df2$today <- substring(df2$datetime, 1, 10) == "2023-07-17"
df2$plugged <- df2$hour < 8 | df2$hour > 16

gp <- ggplot(subset(df, today), aes(group=hour)) +
    coord_cartesian(xlim=c(1, 26)) +
    geom_blank()

for (hh in 0:23) {
    gp <- gp +
        geom_density(data=subset(df, today & plugged & hour == hh),
                     aes(y=soc_plugged, x=..density.. / 10),
                     position=position_nudge(x=hh + 1),
                     fill=scales::hue_pal()(24)[hh + 1],
                     alpha=0.6, adjust=1) +
        geom_density(data=subset(df2, today & plugged & hour == hh),
                     aes(y=soc_plugged, x=..density.. / 10),
                     position=position_nudge(x=hh + 1),
                     fill=NA, linetype='dashed', adjust=1)
}

gp + theme_bw() + scale_y_continuous("Plugged-in SOC (%)", labels=scales::percent) + xlab("Hour of the day")
ggsave(paste0("../figures/uniform-ncar-both.pdf"), width=6.5, height=4)

