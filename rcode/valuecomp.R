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
    scale_x_datetime() + theme_bw() + xlab(NULL) +
    scale_fill_discrete(name="Value component:")
