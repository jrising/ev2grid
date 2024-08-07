setwd("~/research/ev2grid/ev2grid")

library(dplyr)
library(lubridate)
library(ggplot2)
library(lfe)
library(timeDate)
library(FNN)
library(reshape2)
source("~/projects/research-common/R/felm-tools.R")

year.min <- 2018
year.max <- 2024

df <- data.frame()
for (year in year.min:(year.max-1))
    df <- rbind(df, read.csv(paste0("prices-", year, "-", (year+1) %% 100, ".csv")))
df$datetime <- parse_date_time(df$datetime_beginning_ept, "mdY HMS Op", tz="America/New_York")
df$yyyymmdd <- format(df$datetime, "%Y%m%d")

## WILMINGTON NEW CASTLE CO AIRPORT, DE US
## https://www.ncei.noaa.gov/pub/data/ghcn/daily/by_station/USW00013781.csv.gz
weather <- read.csv("USW00013781.csv", header=F, col.names=c('id', 'yyyymmdd', 'element', 'value', 'mf', 'qf', 'sf', 'time'))
weather2 <- dcast(weather, yyyymmdd ~ element)
weather2$yyyymmdd <- as.character(weather2$yyyymmdd)

df$holiday <- F
for (datename in c('USNewYearsDay', 'USInaugurationDay', 'USMLKingsBirthday', 'USLincolnsBirthday', 'USWashingtonsBirthday', 'USMemorialDay', 'USIndependenceDay', 'USLaborDay', 'USColumbusDay', 'USElectionDay', 'USVeteransDay', 'USThanksgivingDay', 'USChristmasDay', 'USCPulaskisBirthday', 'USGoodFriday'))
    df$holiday[as.Date(df$datetime) %in% as.Date(holiday(year.min:year.max, datename)@Data)] <- T
df$weekday <- weekdays(df$datetime)
df$yday <- yday(df$datetime)
df$hour <- hour(df$datetime)
df$yday.cos <- cos(df$yday * 2*pi / 365.25)
df$yday.sin <- sin(df$yday * 2*pi / 365.25)
df$hour.cos <- cos(df$hour * 2*pi / 24)
df$hour.sin <- sin(df$hour * 2*pi / 24)

df2 <- df %>% left_join(weather2) %>%
    arrange(datetime) %>% mutate(lag1=lag(rmccp, 1), lag2=lag(rmccp, 2),
                                 lag3=lag(rmccp, 3), lag4=lag(rmccp, 4),
                                 lag1d=lag(rmccp, 24), lag2d=lag(rmccp, 48),
                                 lag3d=lag(rmccp, 24*3), lag4d=lag(rmccp, 24*4),
                                 lag5d=lag(rmccp, 24*5), lag6d=lag(rmccp, 24*6),
                                 lag7d=lag(rmccp, 24*7),
                                 pcplag2d=lag(rmpcp, 24*2),
                                 pcplag3d=lag(rmpcp, 24*3), pcplag4d=lag(rmpcp, 24*4),
                                 pcplag5d=lag(rmpcp, 24*5), pcplag6d=lag(rmpcp, 24*6),
                                 pcplag7d=lag(rmpcp, 24*7),
                                 rtlag2d=lag(total_pjm_rt_load_mwh, 24*2),
                                 loclag2d=lag(total_pjm_loc_credit, 24*2),
                                 rplag2d=lag(total_pjm_reg_purchases, 24*2),
                                 ssrlag2d=lag(total_pjm_self_sched_reg, 24*2),
                                 arlag2d=lag(total_pjm_assigned_reg, 24*2),
                                 crlag2d=lag(total_pjm_rmccp_cr, 24*2),
                                 pcpcrlag2d=lag(total_pjm_rmpcp_cr, 24*2),
                                 PRCP1d=lag(PRCP, 1), PRCP2d=lag(PRCP, 2),
                                 TMAX1d=lag(TMAX, 1), TMAX2d=lag(TMAX, 2),
                                 TMIN1d=lag(TMIN, 1), TMIN2d=lag(TMIN, 2))
df2$rmccp[df2$rmccp == 0] <- 0.004 # 0.01 is lowest otherwise

ggplot(df2, aes(datetime, rmccp)) +
    geom_line() + scale_x_datetime() + scale_y_log10()

ggplot(df2, aes(datetime, rmccp)) +
    geom_line() + scale_x_datetime(limits=as.POSIXct(c("2022-06-01", "2022-06-02"))) + scale_y_log10()

summary(felm(log(rmccp) ~ datetime + holiday + weekday + lag1 + lag1d + lag2d + lag3d + lag4d + lag5d + lag6d + lag7d + yday.cos + yday.sin + hour.cos + hour.sin, data=df2))

summary(felm(log(rmccp) ~ datetime + holiday + weekday + lag1 + lag1d + lag2d + lag3d + lag4d + lag5d + lag6d + lag7d | factor(yday) + factor(hour), data=df2))

summary(felm(log(rmccp) ~ datetime + holiday + weekday + lag1 + lag2 + lag1d + lag2d + lag3d + lag4d + lag5d + lag6d + lag7d | factor(yday) + factor(hour), data=df2))

summary(felm(log(rmccp) ~ datetime + holiday + weekday + lag1 + lag2 + lag1d + lag2d + lag3d + lag4d + lag5d + lag6d + lag7d + pcplag2d + pcplag3d + pcplag4d + pcplag5d + pcplag6d + pcplag7d + rtlag2d + loclag2d + rplag2d + ssrlag2d + arlag2d + crlag2d + pcpcrlag2d + TMAX + TMAX1d | factor(yday) + factor(hour), data=df2))

df2$set <- rep(rep(1:4, each=24), ceiling(nrow(df2) / 96))[1:nrow(df2)]

est.knn <- function(df2, covars, knum) {
    knnvals <- rep(NA, nrow(df2))

    knnmat <- df2[, covars]
    if ('weekday' %in% covars) {
        knnmat$weekend <- df2$weekday %in% c('Saturday', 'Sunday')
        knnmat$weekday <- as.numeric(factor(knnmat$weekday, levels=c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")))
    }
    knnmat2 <- as.data.frame(lapply(knnmat, function(x) ((x - min(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T)))))
    knnvalid <- rowSums(is.na(knnmat2)) == 0
    for (trainset in 1:4) {
        trainmat <- knnmat2[df2$set == trainset & knnvalid, ]
        testmat <- knnmat2[df2$set == ((trainset - 1 + 2) %% 4) + 1 & knnvalid, ]
        trainyy <- df2$rmccp[df2$set == trainset & knnvalid]
        knn.res <- knn.reg(trainmat, testmat, trainyy, knum)
        knnvals[df2$set == ((trainset - 1 + 2) %% 4) + 1 & knnvalid] <- knn.res$pred
    }
    knnvals
}

df2$holiday.f <- as.numeric(df2$holiday)

df2$knn <- est.knn(df2, c('holiday.f', 'lag2d', 'lag3d', 'lag4d', 'lag5d', 'lag6d', 'lag7d', 'pcplag2d', 'pcplag3d', 'pcplag4d', 'pcplag5d', 'pcplag6d', 'pcplag7d', 'rtlag2d', 'loclag2d', 'rplag2d', 'ssrlag2d', 'arlag2d', 'crlag2d', 'pcpcrlag2d', 'yday.cos', 'yday.sin', 'hour.cos', 'hour.sin', 'weekday', 'TMAX2d', 'TMIN2d', 'PRCP2d'), 5)

summary(felm(log(rmccp) ~ datetime + holiday + weekday + lag1 + lag2 + lag1d + lag2d + lag3d + lag4d + lag5d + lag6d + lag7d + pcplag2d + pcplag3d + pcplag4d + pcplag5d + pcplag6d + pcplag7d + rtlag2d + loclag2d + rplag2d + ssrlag2d + arlag2d + crlag2d + pcpcrlag2d + TMAX1d + knn | factor(yday) + factor(hour), data=df2))

project.48.felm <- function(df2, covars) {
    felm(as.formula(paste("log(rmccp) ~", paste(covars, collapse=" + "), "| factor(yday) + factor(hour) + factor(weekday)")), data=df2)
}

## Project forward from set i to set i + 2
project.48 <- function(df2, covars, get.se=F, use.mod=NULL) {
    if (is.null(use.mod)) {
        mod <- project.48.felm(df2, covars)
    } else {
        mod <- use.mod
    }

    modfe <- getfe(mod)
    modfe.yday <- subset(modfe, fe == 'factor(yday)')
    modfe.yday$idx <- as.numeric(as.character(modfe.yday$idx))
    modfe.hour <- subset(modfe, fe == 'factor(hour)')
    modfe.hour$idx <- as.numeric(as.character(modfe.hour$idx))
    modfe.weekday <- subset(modfe, fe == 'factor(weekday)')
    modfe.weekday$idx <- as.character(modfe.weekday$idx)

    finalpredicted <- rep(NA, nrow(df2))
    if (get.se)
        finalse <- rep(NA, nrow(df2))

    daybegins <- seq(1, nrow(df2), by=24)
    daybeginsets <- ((1:length(daybegins)) - 1) %% 4 + 1
    for (trainset in 1:4) {
        predicted <- rep(NA, nrow(df2))
        predicted[df2$set == trainset] <- df2$rmccp[df2$set == trainset]
        if (get.se)
            thisse <- rep(NA, nrow(df2))

        predstart <- daybegins[daybeginsets == ((trainset - 1 + 1) %% 4) + 1]
        if (trainset == 4)
            predstart <- predstart[-1]

        for (tt in 0:47) {
            preddf <- df2[predstart + tt, covars]
            if ('lag1' %in% covars)
                preddf$lag1 <- predicted[predstart + tt - 1]
            if ('lag2' %in% covars)
                preddf$lag2 <- predicted[predstart + tt - 2]
            if ('lag1d' %in% covars)
                preddf$lag1d <- predicted[predstart + tt - 24]
            predvalid <- rowSums(is.na(preddf)) == 0
            ydayfe <- sapply(df2$yday[predstart + tt], function(yday) modfe.yday$effect[modfe.yday$idx == yday][1])
            hourfe <- sapply(df2$hour[predstart + tt], function(hour) modfe.hour$effect[modfe.hour$idx == hour][1])
            weekdayfe <- sapply(df2$weekday[predstart + tt], function(weekday) modfe.weekday$effect[modfe.weekday$idx == weekday][1])

            if (get.se) {
                fitse <- predict.felm(mod, preddf, se.fit=T)
                predicted[predstart[predvalid] + tt] <- fitse$fit$fit + (ydayfe + hourfe + weekdayfe)[predvalid]

                thisse[predstart[predvalid] + tt] <- fitse$se.fit
            } else {
                predicted[predstart[predvalid] + tt] <- predict.felm(mod, preddf)$fit + (ydayfe + hourfe + weekdayfe)[predvalid]
            }
        }

        finalpredicted[df2$set == ((trainset - 1 + 2) %% 4) + 1] <- predicted[df2$set == ((trainset - 1 + 2) %% 4) + 1]
        if (get.se)
            finalse[df2$set == ((trainset - 1 + 2) %% 4) + 1] <- thisse[df2$set == ((trainset - 1 + 2) %% 4) + 1]
    }

    if (get.se) {
        return(list(predicted=finalpredicted, se=finalse))
    } else {
        return(finalpredicted)
    }
}

df2$predicted <- project.48(df2, c('datetime', 'holiday.f', 'lag1', 'lag2', 'lag1d', 'lag2d', 'lag3d', 'lag4d', 'lag5d', 'lag6d', 'lag7d', 'pcplag2d', 'pcplag3d', 'pcplag4d', 'pcplag5d', 'pcplag6d', 'pcplag7d', 'rtlag2d', 'loclag2d', 'rplag2d', 'ssrlag2d', 'arlag2d', 'crlag2d', 'pcpcrlag2d', 'TMAX1d', 'knn'))

ggplot(df2, aes(log(rmccp), predicted)) +
    geom_point() + geom_abline(yintercept=0, slope=1, colour='#808080') +
    theme_bw()

allcovars.knn <- c('holiday', 'lag2d', 'lag3d', 'lag4d', 'lag5d', 'lag6d', 'lag7d', 'pcplag2d', 'pcplag3d', 'pcplag4d', 'pcplag5d', 'pcplag6d', 'pcplag7d', 'rtlag2d', 'loclag2d', 'rplag2d', 'ssrlag2d', 'arlag2d', 'crlag2d', 'pcpcrlag2d', 'yday.cos', 'yday.sin', 'hour.cos', 'hour.sin', 'weekday', 'TMAX2d', 'TMIN2d', 'PRCP2d')
allcovars.lfe <- c('datetime', 'holiday.f', 'lag1', 'lag2', 'lag1d', 'lag2d', 'lag3d', 'lag4d', 'lag5d', 'lag6d', 'lag7d', 'pcplag2d', 'pcplag3d', 'pcplag4d', 'pcplag5d', 'pcplag6d', 'pcplag7d', 'rtlag2d', 'loclag2d', 'rplag2d', 'ssrlag2d', 'arlag2d', 'crlag2d', 'pcpcrlag2d', 'TMAX2d', 'TMIN2d', 'PRCP2d', 'knn')

df2$knn <- est.knn(df2, allcovars.knn, 5)
df2$predicted <- project.48(df2, allcovars.lfe)
rmse <- sqrt(mean((log(df2$rmccp) - df2$predicted)^2, na.rm=T))

best.covars.knn <- allcovars.knn
best.covars.lfe <- allcovars.lfe
best.knum <- 5
best.rmse <- rmse
for (ii in 1:1000) {
    print(ii)
    ## Try to mutate: drop 1 and add up to 2
    try.covars.knn <- unique(c(sample(best.covars.knn, length(best.covars.knn)-1), sample(allcovars.knn, 2)))
    try.covars.lfe <- unique(c(sample(best.covars.lfe, length(best.covars.lfe)-1), sample(allcovars.lfe, 2)))
    try.knum <- best.knum + sample(-1:1, 1)
    if (identical(sort(try.covars.knn), sort(best.covars.knn)) && identical(sort(try.covars.lfe), sort(best.covars.lfe))) {
        ## Randomly select
        try.covars.knn <- sample(allcovars.knn, 2 + floor(runif(1) * (length(allcovars.knn) - 2)))
        try.covars.lfe <- sample(allcovars.lfe, 2 + floor(runif(1) * (length(allcovars.lfe) - 2)))
        try.knum <- ceiling(runif(1) * 20)
    }

    if ('knn' %in% try.covars.lfe)
        df2$knn <- est.knn(df2, try.covars.knn, try.knum)
    df2$predicted <- project.48(df2, try.covars.lfe)
    try.rmse <- sqrt(mean((log(df2$rmccp) - df2$predicted)^2, na.rm=T))

    if (try.rmse < best.rmse) {
        print("Improved!")
        best.covars.knn <- try.covars.knn
        best.covars.lfe <- try.covars.lfe
        best.knum <- try.knum
        best.rmse <- try.rmse
    }
}

## load("predprice.RData")

allcovars.knn[!(allcovars.knn %in% best.covars.knn)]
allcovars.lfe[!(allcovars.lfe %in% best.covars.lfe)]

df2$knn <- est.knn(df2, best.covars.knn, best.knum)
df2$predicted <- project.48(df2, best.covars.lfe)
sqrt(mean((log(df2$rmccp) - df2$predicted)^2, na.rm=T))

ggplot(df2, aes(log(rmccp), predicted)) +
    geom_point(alpha=.025) + geom_abline(intercept=0, slope=1, colour='#808080') +
    theme_bw() + xlab("Log of realized regulation price") + ylab("Log of predicted regular price") +
    xlim(-3, 7) + ylim(0, 6)
ggsave("predprice.png", width=6.5, height=5)

save(best.covars.knn, best.knum, best.covars.lfe, file="predprice.RData")

## Add levels version of predicted
resvar = var(log(df2$rmccp) - df2$predicted, na.rm=T)
df2$predpe = exp(df2$predicted + resvar / 2)

## Bootstrap full result
for (bs in 1:20) {
    print(bs)
    df2.bs <- df2[sample(1:nrow(df2), nrow(df2), replace=T),]
    df2.bs$knn <- est.knn(df2.bs, best.covars.knn, best.knum)
    mod <- project.48.felm(df2.bs, best.covars.lfe)
    predlog <- project.48(df2, best.covars.lfe, use.mod=mod)
    resvar = var(log(df2$rmccp) - predlog, na.rm=T)
    df2[, paste0("predbs", bs)] <- exp(predlog + resvar / 2)
}

write.csv(df2[!is.na(df2$predpe), c('datetime', 'rmccp', 'predpe', paste0("predbs", 1:20))], "predprice.csv", row.names=F)
