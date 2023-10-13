library(foreach) 
library(doParallel)
library(ProFound)
library(ProSpect)
library(Cairo)
library(data.table)
library(celestial)
library(dftools)
library(matrixStats)
library(dplyr)
library(magicaxis)
library(hyper.fit)

## make CSFH samples for Kevin Croker using same methods in analysis_linfit.R

set.seed(667)
registerDoParallel(cores = 5)
## Analysis script
## Actually do the analysis here 

galaxies = data.frame( fread("~/Documents/PROJ2_JWST_CSFH/data/catalogues/high_z_profound.csv") )

load_with_agn = fread("~/Documents/PROJ2_JWST_CSFH/data/catalogues/ProSpect_highz_withAGN.csv")
load_without_agn = fread("~/Documents/PROJ2_JWST_CSFH/data/catalogues/ProSpect_highz_withoutAGN.csv")


zedges = c(3.5, 6.5, 9.5, 12.5)

## make table of redshift info 
zlist = c(
  median(load_with_agn$z[load_with_agn$z >= zedges[1] & load_with_agn$z < zedges[2]]),
  median(load_with_agn$z[load_with_agn$z >= zedges[2] & load_with_agn$z < zedges[3]]),
  median(load_with_agn$z[load_with_agn$z >= zedges[3] & load_with_agn$z < zedges[4]])
)

zlist_min = c(
  min(load_with_agn$z[load_with_agn$z >= zedges[1] & load_with_agn$z < zedges[2]]),
  min(load_with_agn$z[load_with_agn$z >= zedges[2] & load_with_agn$z < zedges[3]]),
  min(load_with_agn$z[load_with_agn$z >= zedges[3] & load_with_agn$z < zedges[4]])
)

zlist_max = c(
  max(load_with_agn$z[load_with_agn$z >= zedges[1] & load_with_agn$z < zedges[2]]),
  max(load_with_agn$z[load_with_agn$z >= zedges[2] & load_with_agn$z < zedges[3]]),
  max(load_with_agn$z[load_with_agn$z >= zedges[3] & load_with_agn$z < zedges[4]])
)

zlist_lo = c(
  quantile(load_with_agn$z[load_with_agn$z >= zedges[1] & load_with_agn$z < zedges[2]], 0.16),
  quantile(load_with_agn$z[load_with_agn$z >= zedges[2] & load_with_agn$z < zedges[3]], 0.16),
  quantile(load_with_agn$z[load_with_agn$z >= zedges[3] & load_with_agn$z < zedges[4]], 0.16)
)

zlist_hi = c(
  quantile(load_with_agn$z[load_with_agn$z >= zedges[1] & load_with_agn$z < zedges[2]], 0.84),
  quantile(load_with_agn$z[load_with_agn$z >= zedges[2] & load_with_agn$z < zedges[3]], 0.84),
  quantile(load_with_agn$z[load_with_agn$z >= zedges[3] & load_with_agn$z < zedges[4]], 0.84)
)

redshift_data = data.frame(list(
  z = zlist,
  z16 = zlist - zlist_lo,
  z84 = zlist_hi - zlist,
  zmin = zlist_min,
  zmax = zlist_max
))

redshift_data$z

dlogM = 0.01
test_mstar = seq(5.0, 12.0, dlogM) #with agn

smf_z5 = dfmodel(test_mstar, p = c(log10(5.16e-5), 10.97, -1.70)) # https://arxiv.org/pdf/1507.05636.pdf
smf_z6 = dfmodel(test_mstar, p = c(-4.09, 10.24, -1.88)) # https://arxiv.org/pdf/2103.16571.pdf
smf_z7 = dfmodel(test_mstar, p = c(-4.14, 10.04, -1.73))
smf_z7 = dfmodel(test_mstar, p = c(-4.69, 9.98, -1.82))
smf_z9 = dfmodel(test_mstar, p = c(-5.12, 9.50, -2.00))
smf_z10 = dfmodel(test_mstar, p = c(-6.13, 9.50, -2.00))

spline_fit = function(zedge_min, zedge_max, data, key, df, smf){
  ## keys are fixed
  
  idx = data$z >= zedge_min & data$z < zedge_max
  data = data[idx, ]
  xdata = c(log10(data[["Mstar50"]]) )
  xerr_floor =  0
  xdata_err = c( (data[["Mstar16"]] + data[["Mstar84"]])/(log(10)*data[["Mstar50"]]) )
  xdata_err = sqrt((xerr_floor)^2 + xdata_err^2)
  # xdata_err = rep(0.4, length(xdata))
  
  med_key = paste0(key, "50")
  up_key = paste0(key, "84")
  down_key = paste0(key, "16")
  
  ydata = c( log10(data[[med_key]]) )
  # ydata_err = rep(0.4, length(ydata))
  yerr_floor = 0
  ydata_err = c( (data[[down_key]] + data[[up_key]])/(log(10)*data[[med_key]]) )
  ydata_err = sqrt((yerr_floor)^2 + ydata_err^2)
  
  Nsamples = 1001
  spline_data = foreach(i = 1:Nsamples, .combine = "c")%dopar%{
    midx = xdata >= 8
    
    xx = xdata[midx] + rnorm(length(xdata[midx]), mean = 0, sd = xdata_err[midx])
    yy = ydata[midx] + rnorm(length(ydata[midx]), mean = 0, sd = ydata_err[midx])
    
    fit = smooth.spline(
      xx, yy, df = df, w = 1/(ydata_err[midx]^2 + xdata_err[midx]^2), lambda = 0.8
    )
    out = predict(fit, test_mstar)
    sum(10^out$y * smf) * dlogM
    }
  
  print(smf)
  return(
    list(
      "csfh_samples" = spline_data
      )
    )
  
}

df = 2
fit_data = list()
fit_data$withAGN$z5=spline_fit(zedge_min = zedges[1], zedge_max = zedges[2], 
                               data = load_with_agn, key = "SFRBurst",
                               df = df, smf = smf_z5)

fit_data$withAGN$z7=spline_fit(zedge_min = zedges[2], zedge_max = zedges[3], 
                               data = load_with_agn, key = "SFRBurst",
                               df = df, smf = smf_z7)

fit_data$withAGN$z10=spline_fit(zedge_min = zedges[3], zedge_max = zedges[4], 
                                data = load_with_agn, key = "SFRBurst",
                                df = df, smf = smf_z10)

fit_data$withoutAGN$z5=spline_fit(zedge_min = zedges[1], zedge_max = zedges[2], 
                                  data = load_without_agn, key = "SFRBurst",
                                  df = df, smf = smf_z5)

fit_data$withoutAGN$z7=spline_fit(zedge_min = zedges[2], zedge_max = zedges[3], 
                                  data = load_without_agn, key = "SFRBurst",
                                  df = df, smf = smf_z7)

fit_data$withoutAGN$z10=spline_fit(zedge_min = zedges[3], zedge_max = zedges[4], 
                                   data = load_without_agn, key = "SFRBurst",
                                   df = df, smf = smf_z10)

chabrer_conversion = c(1.0/1.53, 1.0/1.53, 1.0/1.53)


csfh_samples = data.frame(
  "withAGN_z5" = fit_data$withAGN$z5$csfh_samples * 1/1.53,
  "withAGN_z7" = fit_data$withAGN$z7$csfh_samples* 1/1.53,
  "withAGN_z10" = fit_data$withAGN$z10$csfh_samples* 1/1.53,
  "withoutAGN_z5" = fit_data$withoutAGN$z5$csfh_samples* 1/1.53,
  "withoutAGN_z7" = fit_data$withoutAGN$z7$csfh_samples* 1/1.53,
  "withoutAGN_z10" = fit_data$withoutAGN$z10$csfh_samples* 1/1.53
)

redone_csfh = data.frame(
  "withAGN" = c(
    log10(quantile(csfh_samples$withAGN_z5, 0.5)), 
    log10(quantile(csfh_samples$withAGN_z7, 0.5)), 
    log10(quantile(csfh_samples$withAGN_z10, 0.5))
    ),
  "withAGN_q16" = c(
    (quantile(csfh_samples$withAGN_z5, 0.5) - quantile(csfh_samples$withAGN_z5, 0.16))/(log(10)*quantile(csfh_samples$withAGN_z5, 0.5)),
    (quantile(csfh_samples$withAGN_z7, 0.5) - quantile(csfh_samples$withAGN_z7, 0.16))/(log(10)*quantile(csfh_samples$withAGN_z7, 0.5)),
    (quantile(csfh_samples$withAGN_z10, 0.5) - quantile(csfh_samples$withAGN_z10, 0.16))/(log(10)*quantile(csfh_samples$withAGN_z10, 0.5))
    ),
  "withAGN_q84" = c(
    (quantile(csfh_samples$withAGN_z5, 0.84) - quantile(csfh_samples$withAGN_z5, 0.50))/(log(10)*quantile(csfh_samples$withAGN_z5, 0.5)),
    (quantile(csfh_samples$withAGN_z7, 0.84) - quantile(csfh_samples$withAGN_z7, 0.50))/(log(10)*quantile(csfh_samples$withAGN_z7, 0.5)),
    (quantile(csfh_samples$withAGN_z10, 0.84) - quantile(csfh_samples$withAGN_z10, 0.50))/(log(10)*quantile(csfh_samples$withAGN_z10, 0.5))
    ),
  "withoutAGN" = c(
    log10(quantile(csfh_samples$withoutAGN_z5, 0.5)),
    log10(quantile(csfh_samples$withoutAGN_z7, 0.5)),
    log10(quantile(csfh_samples$withoutAGN_z10, 0.5))
  ),
  "withoutAGN_q16" = c(
    (quantile(csfh_samples$withoutAGN_z5, 0.5) - quantile(csfh_samples$withoutAGN_z5, 0.16))/(log(10)*quantile(csfh_samples$withoutAGN_z5, 0.5)),
    (quantile(csfh_samples$withoutAGN_z7, 0.5) - quantile(csfh_samples$withoutAGN_z7, 0.16))/(log(10)*quantile(csfh_samples$withoutAGN_z7, 0.5)),
    (quantile(csfh_samples$withoutAGN_z10, 0.5) - quantile(csfh_samples$withoutAGN_z10, 0.16))/(log(10)*quantile(csfh_samples$withoutAGN_z10, 0.5))
  ),
  "withoutAGN_q84" = c(
    (quantile(csfh_samples$withoutAGN_z5, 0.84) - quantile(csfh_samples$withoutAGN_z5, 0.50))/(log(10)*quantile(csfh_samples$withoutAGN_z5, 0.5)),
    (quantile(csfh_samples$withoutAGN_z7, 0.84) - quantile(csfh_samples$withoutAGN_z7, 0.50))/(log(10)*quantile(csfh_samples$withoutAGN_z7, 0.5)),
    (quantile(csfh_samples$withoutAGN_z10, 0.84) - quantile(csfh_samples$withoutAGN_z10, 0.50))/(log(10)*quantile(csfh_samples$withoutAGN_z10, 0.5))
  )
)


magplot(
  read_csfh$z, 
  read_csfh$CSFH_with_AGN,
  col="red",
  xlim = c(0,16),
  ylim = c(-4.5, 6)
)
magerr(
  read_csfh$z, 
  read_csfh$CSFH_with_AGN,
  ylo=read_csfh$CSFH_with_AGN_q16,
  yhi=read_csfh$CSFH_with_AGN_q84,
  col="red"
)
points(
  read_csfh$z, 
  read_csfh$CSFH_without_AGN,
  col="blue"
)
magerr(
  read_csfh$z, 
  read_csfh$CSFH_without_AGN,
  ylo=read_csfh$CSFH_without_AGN_q16,
  yhi=read_csfh$CSFH_without_AGN_q84,
  col="blue"
)
points(
  read_csfh$z+0.05,
  redone_csfh$withAGN,
  pch = 16,
  col="red"
)
magerr(
  read_csfh$z+0.05,
  redone_csfh$withAGN,
  ylo = redone_csfh$withAGN_q16,
  yhi = redone_csfh$withAGN_q84,
  pch = 16,
  col="red"
)

points(
  read_csfh$z+0.05,
  redone_csfh$withoutAGN,
  pch = 16,
  col="blue"
)
magerr(
  read_csfh$z+0.05,
  redone_csfh$withoutAGN,
  ylo = redone_csfh$withoutAGN_q16,
  yhi = redone_csfh$withoutAGN_q84,
  pch = 16,
  col="blue"
)

fwrite(csfh_samples, "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/csfh_samples.csv")
