# Power analysis
# author: Coilin/Olga
# date: 15-03-2022
# area: RtD-2013

set.seed(1234)
library(ggplot2)
library(Distance)
library(gratia)
library(mgcv)

###---------------
# change these lines
area.year <- "RtD-2013"
# total number transects/number survey days
ntrans <- round(129 / 6)
###---------------

if (grepl("coilin", getwd())) {
  path <- "./"
} else{
  path <- "../phase1-distance-modelling/rds/"
}

## create the folders
folders <- c("rds", "figs")

for(i in folders){
    if(!file.exists(i)){
        system(paste("mkdir", i))
    }
}


## load dist model
file.name <- paste0("detfc-hn-cos-", area.year, ".rds")
detfc.hn.cos <- readRDS(paste0(path, file.name))

## load gam model
file.name <- paste0("dsm-ours-", area.year, ".rds")
mod <- readRDS(paste0(path, file.name))

## extract the standard deviations of the random effects
re_var_df <- variance_comp(mod)
sd_day <- as.numeric(re_var_df[re_var_df$component == "s(fDay)", "std_dev"])
sd_trans <- as.numeric(re_var_df[re_var_df$component == "s(fTransect)", "std_dev"])
b0 <- coef(mod)["(Intercept)"]
area <- detfc.hn.cos$dht$individuals$summary$Area
conversion.factor <- convert_units("meter", "meter", "square kilometer")

source("rzipch_lambda.R")

## half-normal
hn <- function(y, sigma) {
  exp(-y ^ 2 / (2 * sigma ^ 2))
}

sim_daily_obs_data <- function(nd, nt, p_change) {
  ##-------------------------------
  ## nd is the number of days
  ## nt is the number of transects
  ## p_change change
  ##-------------------------------
  re_day <- rnorm(nd, sd = sd_day)
  dat <- expand.grid(day = 1:nd, transect = 1:nt)
  dat$re_day <- re_day[dat$day]
  dat$re_trans <- rnorm(nd * nt, sd = sd_trans)
  dat$offset <- median(mod$offset)
  dat$b0 <- coef(mod)["(Intercept)"]
  ## linear predictor
  dat$lp <- with(dat, b0 + re_day + re_trans + offset)
  ## FIX 1
  # dat$lp[dat$lp > 5] <- 5
  # dat$lp[dat$lp < -2.5] <- -2.5
  ## parameters or the zero-inflated poisson
  ## USE optim
  dat$count <- rzipch2(gamma = dat$lp, p_change)
  ## add in the observation process
  dat_zero <- subset(dat, count == 0)
  dat_pos <- subset(dat, count > 0)
  ## repeat for individuals
  dat_pos_long <- dat_pos[rep(1:nrow(dat_pos), times = dat_pos$count), ]
  dat_pos_long$count <- 1
  dat_long <- rbind(dat_zero, dat_pos_long)
  ## add in detection
  lnscale <- as.numeric(summary(detfc.hn.cos)$ds$coeff$key.scale["estimate"])
  scale <- exp(lnscale)
  dat_long$distance <- runif(nrow(dat_long), 0, 300)
  ## probability of detection
  dat_long$pdetect <- hn(dat_long$distance, scale)
  ## sample detections
  dat_long$obscount <- rbinom(nrow(dat_long), size = 1, prob = dat_long$pdetect)
  dat_long$obscount[dat_long$count == 0] <- 0
  ## replace zeros with NA
  dat_long$size <- dat_long$obscount
  dat_long$size[dat_long$size == 0] <- NA
  dat_long$distance[is.na(dat_long$size)] <- NA
  # check value
  dat_long$Effort <- 100 # avg transect size
  dat_long$Area <- area
  return(dat_long)
}

# Estimate the significance of the slopes by resampling from the means and standard errors.
# robust to two or more means but will add some computational time for the resample

test_slope <- function(mu, se, nrep = 1e3) {
  ##---------------------
  ## mu is a vector of means
  ## se is a vector of standard errors
  ## nrep is the number of replicates used to
  ## build the distribution of the slope
  ##---------------------
  if (length(mu) != length(se)) {
    stop("length of mu must equal length of se")
  }
  ## simulate the data by drawing randomly from the means
  m <- length(mu)
  y <- sapply(1:m, function(x) {
    rnorm(nrep, mean = mu[x], sd = se[x])
  })
  ## calculate the slopes
  slopes <- apply(y, 1, function(x) {
    coef(lm(x ~ I(1:m)))[2]
  })
  ## one-sided test if the slope is less than zero
  ## slopes less than zero?
  sl0 <- factor(slopes < 0, levels = c("TRUE", "FALSE"))
  res <- ifelse(table(sl0)["TRUE"] / nrep > 0.95, 1, 0)
  return(as.numeric(res))
}


#################################
# tp - total % change between the first year and the last year (e.g., 20% decline is -0.2)
# sf - survey frequency in reporting period
# nd - number of days per survey
# nt - number of transects per survey day
# ny - reporting period - time span in years (e.g., 2020-2025 for 6 years)
# sy - survey years
# nsim - number of simulations
#################################

nsim <- 200 ## change to 500
outer_df <- expand.grid(
    sf = c(2, 3, 6),
    nd = c(6, 8, 10),
    nt = c(ntrans, round(ntrans * 1.2)),
    tp = seq(-0.01,-0.81, by = -0.1),
    ##nd = c(6, 10),
    ## number transects + 0.20%
    iter = 1:nsim,
    ny = 6 # total number
)
n <- nrow(outer_df)

######################
## parallel component
######################
library(doSNOW)
start_time <- Sys.time()

## register nodes 
## change as required
n_node <- 10
##cl <- makeSOCKcluster(n_node)
cl <- makeCluster(n_node, outfile = "")
registerDoSNOW(cl)

## ## progress bar
## pb <- txtProgressBar(max = n, style = 3)
## progress <- function(n)
##   setTxtProgressBar(pb, n)
## opts <- list(progress = progress)

cat("", file = "progress.txt")

## parallel loop
##system.time({
    all_sim_res <-
        ##  foreach(i = 1:n,
        foreach(i = 1:n,
                .combine = rbind
                ##.options.snow = opts) %dopar% {
                ) %dopar% {
                    set.seed(i)
                    library(Distance)
                    library(gratia)
                    library(mgcv)
                    library(R.utils)
                    res <- try({
                        ## equation of line is through points (1,0) and (ny, tp)
                                        # yearly % change
                        slope <- outer_df[i, ]$tp / (outer_df[i, ]$ny - 1)
                        intercept <- -outer_df[i, ]$tp / (outer_df[i, ]$ny - 1)
                                        # change!
                        sy <- seq(1, outer_df[i, ]$ny, length = outer_df[i, ]$sf)
                        seqp <- intercept + slope * sy
                        sim <- sapply(seqp, function(z) {
                            sim_daily_obs_data(nd = outer_df[i, ]$nd,
                                        # all the same
                                               nt = outer_df[i, ]$nt,
                                               p_change = z)
                        }, simplify = FALSE)
                        ## detection
                        nhat_list <- lapply(sim, function(z) {
                            ## fit <- ds(
                            ##   data = z,
                            ##   transect = "line",
                            ##   key = "hn",
                            ##   formula = ~ 1,
                            ##   adjustment = "cos",
                            ##   convert.units = conversion.factor,
                            ##   debug.level = 0,
                            ##   truncation = 300,
                            ##   quiet = TRUE
                            ## )
                            ## to avoid running for too long
                            ## catch a long run > 30 seconds and give it an NA for change
                            fit <- withTimeout({
                                ds(
                                    data = z,
                                    transect = "line",
                                    key = "hn",
                                    formula = ~ 1,
                                    adjustment = "cos",
                                    convert.units = conversion.factor,
                                    debug.level = 0,
                                    truncation = 300,
                                    quiet = TRUE
                                )
                            }, timeout = 30, onTimeout = "silent")
                            if(is.null(fit)){
                                nhat <- NA
                                se <- NA
                            }else{
                                nhat <- summary(fit$ddf)$Nhat
                                se <- summary(fit$ddf)$Nhat.se
                            }
                            res <- data.frame(nhat = nhat, se = se)
                            return(res)
                        })
                        nhat_df <- do.call(rbind, nhat_list)
                        ## significance
                        if(!any(is.na(nhat_df$nhat))){
                            change <- test_slope(mu = nhat_df$nhat,
                                                 se = nhat_df$se,
                                                 nrep = 1e3)
                        }else{
                            change <- NA
                        }
                        ## show progress
                        ## cat(paste0(round(i / n * 100, 3), '% completed'))
                        cat(paste(i, "\n"), file = "progress.txt", append = TRUE)
                        res <- data.frame(i = i, det = change)
                        res
                    })                    
                    if (class(res) == "try-error") {
                        res <- data.frame(i = i, det = NA)
                    }
                    res
                }
##})
end_time <- Sys.time()
end_time - start_time
    
# close(pb)
stopCluster(cl)

## check that the index correct
all(diff(all_sim_res$i) == 1)

outer_df <- cbind(outer_df, all_sim_res)
saveRDS(outer_df, file = paste0("rds/outer-df-", area.year, ".rds"))

outer_agg <- aggregate(det ~ tp + nd + nt + sf, data = outer_df, sum)
outer_agg$power <- outer_agg$det / nsim

saveRDS(outer_agg, file = paste0("rds/outer-agg-", area.year, ".rds"))





