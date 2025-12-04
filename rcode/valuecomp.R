## setwd("~/research/ev2grid/ev2grid")

library(tidyr)
library(dplyr)
library(ggplot2)

df <- read.csv("results/bytime.csv")

df2 <- df %>% filter(as.Date(datetime) == "2023-07-17") %>%
    pivot_longer(cols=c('valuep', 'valuepns', 'valuee', 'valuer'))
df2$label <- unlist(list('valuep'="Power costs", 'valuee'="Energy level",
                         'valuepns'="State penalty", 'valuer'="Regulation payment")[df2$name])

ggplot(df2, aes(as.POSIXct(datetime, format = "%Y-%m-%dT%H:%M:%S") + 60*30, value, fill=label)) +
    facet_wrap(~ Approach, ncol=1) +
    geom_col() +
    geom_vline(xintercept=c(as.POSIXct("2023-07-17 09:00:00"),
                            as.POSIXct("2023-07-17 17:00:00"))) +
    scale_x_datetime(date_labels="%H:%M", expand=c(0, 0),
                     limits=c(as.POSIXct("2023-07-17 00:00:00"),
                              as.POSIXct("2023-07-18 00:00:00"))) + theme_bw() + xlab(NULL) +
    scale_fill_discrete(name="Value component:")

ggplot(df %>% group_by(datetime, Approach) %>% summarize(value=valuep), aes(as.POSIXct(datetime, format = "%Y-%m-%dT%H:%M:%S") + 60*30, value, colour=Approach)) +
    geom_line() +
    geom_vline(xintercept=c(as.POSIXct("2023-07-17 09:00:00"),
                            as.POSIXct("2023-07-17 17:00:00"))) +
    scale_x_datetime(date_labels="%H:%M", expand=c(0, 0),
                     limits=c(as.POSIXct("2023-07-17 00:00:00"),
                              as.POSIXct("2023-07-18 00:00:00"))) + theme_bw() + xlab(NULL) +
    scale_fill_discrete(name="Value component:")

df <- read.csv("results/bytime-xstart.csv")

ggplot(df %>% filter(as.Date(datetime) == "2023-07-17") %>% group_by(start_hour, Approach) %>%
       summarize(value=sum(valuep)), aes(start_hour, value, colour=Approach)) +
    geom_line() +
    geom_vline(xintercept=c(12, 20)) +
    scale_x_continuous("Driving start hour", expand=c(0, 0),
                       breaks=seq(4, 24, by=4), labels=c(paste(seq(4, 11, by=4), "AM"), paste(seq(12, 24, by=4), "PM"))) +
                                                    theme_bw() + xlab(NULL) +
    ylab("Value of energy arbitrage")

df2 <- df %>% filter(as.Date(datetime) == "2023-07-17" & start_hour == 12) %>%
    pivot_longer(cols=c('valuep', 'valuepns', 'valuee', 'valuer'))
df2$label <- unlist(list('valuep'="Power costs", 'valuee'="Energy level",
                         'valuepns'="State penalty", 'valuer'="Regulation payment")[df2$name])

ggplot(df2, aes(as.POSIXct(datetime, format = "%Y-%m-%dT%H:%M:%S") + 60*30, value, fill=label)) +
    facet_wrap(~ Approach, ncol=1) +
    geom_col() +
    geom_vline(xintercept=c(as.POSIXct("2023-07-17 09:00:00"),
                            as.POSIXct("2023-07-17 17:00:00"))) +
    scale_x_datetime(date_labels="%H:%M", expand=c(0, 0),
                     limits=c(as.POSIXct("2023-07-17 00:00:00"),
                              as.POSIXct("2023-07-18 00:00:00"))) + theme_bw() + xlab(NULL) +
    scale_fill_discrete(name="Value component:")
