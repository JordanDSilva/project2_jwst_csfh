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

set.seed(667)
registerDoParallel(cores = 5)
## Analysis script
## Actually do the analysis here 

galaxies = data.frame( fread("~/Documents/PROJ2_JWST_CSFH/data/catalogues/high_z_profound.csv") )

load_with_agn = fread("~/Documents/PROJ2_JWST_CSFH/data/catalogues/ProSpect_highz__withAGN.csv")
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
  z16 = zlist_lo,
  z74 = zlist_hi,
  zmin = zlist_min,
  zmax = zlist_max
))

redshift_data$z

dlogM = 0.01
test_mstar = seq(5.0, 12.0, dlogM) #with agn

spline_fit = function(zedge_min, zedge_max, data, key, df){
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
  spline_data = foreach(i = 1:Nsamples, .combine = "rbind")%dopar%{
    midx = xdata >= 8

    xx = xdata[midx] + rnorm(length(xdata[midx]), mean = 0, sd = xdata_err[midx])
    yy = ydata[midx] + rnorm(length(ydata[midx]), mean = 0, sd = ydata_err[midx])

    fit = smooth.spline(
      xx, yy, df = df, w = 1/(ydata_err[midx]^2 + xdata_err[midx]^2), lambda = 0.8
    )
    out = predict(fit, test_mstar)
    out$y
  }
  
  print(dim(spline_data))
  med_fit = colQuantiles(matrix(spline_data, nrow = Nsamples, ncol = length(test_mstar)), 
                         probs = 0.5)
  q16_fit = colQuantiles(matrix(spline_data, nrow = Nsamples, ncol = length(test_mstar)), 
                         probs = 0.16)
  q84_fit = colQuantiles(matrix(spline_data, nrow = Nsamples, ncol = length(test_mstar)), 
                         probs = 0.84)
  return(
    list(
      "test_mstar" = test_mstar,
      "med_fit" = med_fit,
      "q16_fit" = q16_fit,
      "q84_fit" = q84_fit,
      "input_data" = list(
        mstar = xdata,
        mstar_err = xdata_err,
        ydata = ydata,
        ydata_err = ydata_err
      )
    )
  )
}

df = 2
fit_data = list()
fit_data$withAGN$z5=spline_fit(zedge_min = zedges[1], zedge_max = zedges[2], 
                               data = load_with_agn, key = "SFRBurst",
                               df = df)
fit_data$withAGN$z7=spline_fit(zedge_min = zedges[2], zedge_max = zedges[3], 
                                 data = load_with_agn, key = "SFRBurst",
                                 df = df)
fit_data$withAGN$z10=spline_fit(zedge_min = zedges[3], zedge_max = zedges[4], 
                                 data = load_with_agn, key = "SFRBurst",
                                 df = df)

fit_data$withoutAGN$z5=spline_fit(zedge_min = zedges[1], zedge_max = zedges[2], 
                                 data = load_without_agn, key = "SFRBurst",
                                 df = df)
fit_data$withoutAGN$z7=spline_fit(zedge_min = zedges[2], zedge_max = zedges[3], 
                                 data = load_without_agn, key = "SFRBurst",
                                 df = df)
fit_data$withoutAGN$z10=spline_fit(zedge_min = zedges[3], zedge_max = zedges[4], 
                                 data = load_without_agn, key = "SFRBurst",
                                 df = df)

smf_z5 = dfmodel(fit_data$withAGN$z5$test_mstar, p = c(log10(5.16e-5), 10.97, -1.70)) # https://arxiv.org/pdf/1507.05636.pdf
smf_z6 = dfmodel(fit_data$withAGN$z5$test_mstar, p = c(-4.09, 10.24, -1.88)) # https://arxiv.org/pdf/2103.16571.pdf
smf_z7 = dfmodel(fit_data$withAGN$z5$test_mstar, p = c(-4.14, 10.04, -1.73))
smf_z7 = dfmodel(fit_data$withAGN$z5$test_mstar, p = c(-4.69, 9.98, -1.82))
smf_z9 = dfmodel(fit_data$withAGN$z5$test_mstar, p = c(-5.12, 9.50, -2.00))
smf_z10 = dfmodel(fit_data$withAGN$z5$test_mstar, p = c(-6.13, 9.50, -2.00))

smf_z5_cosmos = dfmodel(fit_data$withAGN$z5$test_mstar, p = c(log10(0.14e-3), 10.30, -1.46)) # https://arxiv.org/pdf/2212.02512.pdf
smf_z6_cosmos = dfmodel(fit_data$withAGN$z5$test_mstar, p = c(log10(0.06e-3), 10.14, -1.46)) # Weaver+23
smf_z7_cosmos = dfmodel(fit_data$withAGN$z5$test_mstar, p = c(log10(0.03e-3), 10.18, -1.46)) # Uses CHABRIER IMF


smfs = list(smf_z5, smf_z7, smf_z10)
# smfs = list(smf_z5_cosmos, smf_z7_cosmos, smf_z10)
# chabrer_conversion = c(1.0, 1.0, 1.0/1.53)
chabrer_conversion = c(1.0/1.53, 1.0/1.53, 1.0/1.53)

png("/Users/22252335/Desktop/smfs.png", width = 12, height = 6, res = 240, units="in")
par(mfrow = c(1,2), oma = c(3.5, 3.5, 1.5, 1.5), mar = rep(2.5, 4))
magplot(
  test_mstar, smf_z5 * 1/1.53, log="y", type="l", xlim = c(6.5, 12.5), ylim = 10^c(-6.5,0.5),
  main = "z=5"
)
lines(
  test_mstar, smf_z5_cosmos, lty=2
)

magplot(
  test_mstar, smf_z7 * 1/1.53, log="y", type="l", xlim = c(6.5, 12.5), ylim = 10^c(-6.5,0.5),
    main = "z=7"
)
lines(
  test_mstar, smf_z7_cosmos, lty=2
)
legend(
  x = "topright",
  c("Song+15 (z=5), Stefanon+21 (z=7)",
    "Weaver=23 COSMOS2020"),
  lty = c(1,2)
)
dev.off()




# devilsd10 = fread("/Volumes/JordanData/PhD/GAMA-DEVILS/CSFRD Compendium/GAMA-DEVILS/devilsd10_super.csv")

calc_csfh = function(data, smf, convert_chabrier=1){
  
  idx = data$test_mstar >= 0.0
  
  med_csfh = sum( smf[idx] * 10^data$med_fit[idx] * dlogM )  
  q16_csfh = sum( smf[idx] * 10^data$q16_fit[idx] * dlogM ) 
  q84_csfh = sum( smf[idx] * 10^data$q84_fit[idx] * dlogM ) 
  
  list(
    # "med_csfh" = log10(med_csfh) + log10(1/1.53),
    "med_csfh" = log10(med_csfh * convert_chabrier),
    "q16_csfh" = (med_csfh - q16_csfh)/(log(10)*med_csfh),
    "q84_csfh" = (q84_csfh - med_csfh)/(log(10)*med_csfh) 
  )
}

redshift_data

csfh = list()

csfh$withAGN$z5 = calc_csfh(fit_data$withAGN$z5, smfs[[1]], convert_chabrier = chabrer_conversion[1])
csfh$withAGN$z7 = calc_csfh(fit_data$withAGN$z7, smfs[[2]], convert_chabrier = chabrer_conversion[2])
csfh$withAGN$z10 = calc_csfh(fit_data$withAGN$z10, smfs[[3]], convert_chabrier = chabrer_conversion[3])

csfh$withoutAGN$z5 = calc_csfh(fit_data$withoutAGN$z5, smfs[[1]], convert_chabrier = chabrer_conversion[1])
csfh$withoutAGN$z7 = calc_csfh(fit_data$withoutAGN$z7, smfs[[2]], convert_chabrier = chabrer_conversion[2])
csfh$withoutAGN$z10 = calc_csfh(fit_data$withoutAGN$z10, smfs[[3]], convert_chabrier = chabrer_conversion[3])

csfh$withAGN$med = c(
  csfh$withAGN$z5$med_csfh,
  csfh$withAGN$z7$med_csfh,
  csfh$withAGN$z10$med_csfh
)
csfh$withAGN$q16 = c(
  csfh$withAGN$z5$q16_csfh,
  csfh$withAGN$z7$q16_csfh,
  csfh$withAGN$z10$q16_csfh
)
csfh$withAGN$q84 = c(
  csfh$withAGN$z5$q84_csfh,
  csfh$withAGN$z7$q84_csfh,
  csfh$withAGN$z10$q84_csfh
)

csfh$withoutAGN$med = c(
  csfh$withoutAGN$z5$med_csfh,
  csfh$withoutAGN$z7$med_csfh,
  csfh$withoutAGN$z10$med_csfh
)
csfh$withoutAGN$q16 = c(
  csfh$withoutAGN$z5$q16_csfh,
  csfh$withoutAGN$z7$q16_csfh,
  csfh$withoutAGN$z10$q16_csfh
)
csfh$withoutAGN$q84 = c(
  csfh$withoutAGN$z5$q84_csfh,
  csfh$withoutAGN$z7$q84_csfh,
  csfh$withoutAGN$z10$q84_csfh
)

harikane_22_constant_sfe = function(z){
  ( 61.7 * (1+z)^-3.13 +
    1.0*10^(0.22*(1+z)) +
    2.4 * 10^(0.5*(1+z)-3.0) )^-1
}

# bouwens_2015 = fread("/Volumes/JordanData/PhD/GAMA-DEVILS/CSFRD Compendium/processed/obs/rest_UV/bouwens2015.csv")

# dsilva23_csfh = fread("/Volumes/JordanData/PhD/GAMA-DEVILS/CSFRD Compendium/processed/obs/GAMA-DEVILS/CSFH_DSILVA+23_Ver_Final.csv")

spline_fit_mstar = function(zedge_min, zedge_max, with_agn, without_agn, df = 2){
  ## keys are fixed
  
  idx = with_agn$z >= zedge_min & with_agn$z < zedge_max
  with_agn = with_agn[idx, ]
  xdata = log10(with_agn[["Mstar50"]])
  xerr_floor =  0
  xdata_err = (with_agn[["Mstar16"]] + with_agn[["Mstar84"]])/(log(10)*with_agn[["Mstar50"]])
  xdata_err = sqrt((xerr_floor)^2 + xdata_err^2)
  # xdata_err = pmax(xdata_err, 1)
  
  idx = without_agn$z >= zedge_min & without_agn$z < zedge_max 
  without_agn = without_agn[idx, ]
  ydata = log10(without_agn[["Mstar50"]])
  yerr_floor = 0
  ydata_err = (without_agn[["Mstar16"]] + without_agn[["Mstar84"]])/(log(10)*without_agn[["Mstar50"]])
  ydata_err = sqrt((yerr_floor)^2 + ydata_err^2)
  # ydata_err = pmax(ydata_err, 1)
  
  Nsamples = 1001
  spline_data = foreach(i = 1:Nsamples, .combine = "rbind")%dopar%{

    midx = xdata>=8 & ydata>=8

    xx = xdata[midx] + rnorm(length(xdata[midx]), mean = 0, sd = xdata_err[midx])
    yy = ydata[midx] + rnorm(length(ydata[midx]), mean = 0, sd = ydata_err[midx])
     # idx = xx>=8 & yy>=8
    
    fit = smooth.spline(
      xx, yy, df = df, w = 1/(ydata_err[midx]^2 + xdata_err[midx]^2)
    )
    out = predict(fit, test_mstar)
    out$y
  }
  

  print(dim(spline_data))
  med_fit = colQuantiles(matrix(spline_data, nrow = Nsamples, ncol = length(test_mstar)), 
                         probs = 0.5)
  q16_fit = colQuantiles(matrix(spline_data, nrow = Nsamples, ncol = length(test_mstar)), 
                         probs = 0.16)
  q84_fit = colQuantiles(matrix(spline_data, nrow = Nsamples, ncol = length(test_mstar)), 
                         probs = 0.84)
  return(
    list(
      "test_mstar" = test_mstar,
      "fit_mstar" = med_fit,
      "q16_fit" = q16_fit,
      "q84_fit" = q84_fit,
      
      "delta_mstar" = test_mstar - med_fit,
      "delta_mstar_q16_fit" = test_mstar - q16_fit,
      "delta_mstar_q84_fit" = test_mstar - q84_fit,
      
      # "hf_fit" = hf$func(cbind(test_mstar)),
      # "hf_fit_q16" = hf$func(cbind(test_mstar)) - hf$parm[3],
      # "hf_fit_q84" = hf$func(cbind(test_mstar)) + hf$parm[3],
      # 
      "input_data" = list(
        with_agn = xdata,
        with_agn_err = xdata_err,
        without_agn = ydata,
        without_agn_err = ydata_err,
        delta_mstar = xdata - ydata,
        delta_mstar_err = sqrt(ydata_err^2 + xdata_err^2)
      ),
      "spline_fun" = splinefun(x = test_mstar, med_fit)
    )
  )
}

df_mstar = 3
mstar_fit_data = list()
mstar_fit_data$z5=spline_fit_mstar(3.5, 6.5, load_with_agn, load_without_agn, df_mstar)
mstar_fit_data$z7=spline_fit_mstar(6.5, 9.5, load_with_agn, load_without_agn, df_mstar)
mstar_fit_data$z10=spline_fit_mstar(9.5, 12.5, load_with_agn, load_without_agn, df_mstar)


fwrite(
  redshift_data,
  "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/redshift_info.csv"
)

csfh_data=data.frame(
  "CSFH_with_AGN" = csfh$withAGN$med,
  "CSFH_with_AGN_q16" = csfh$withAGN$q16,
  "CSFH_with_AGN_q84" = csfh$withAGN$q84,
  
  "CSFH_without_AGN" = csfh$withoutAGN$med,
  "CSFH_without_AGN_q16" = csfh$withoutAGN$q16,
  "CSFH_without_AGN_q84" = csfh$withoutAGN$q84,
  
  redshift_data
  
)

fwrite(csfh_data,
       "/Users/22252335/Documents/PROJ2_JWST_CSFH//data/save/jwst_csfh_data.csv")

mstar_mstar_fits = data.frame(
  "test_mstar" = mstar_fit_data$z5$test_mstar,
  
  "z5_fit" = mstar_fit_data$z5$fit_mstar,
  "z5_q16" = mstar_fit_data$z5$q16_fit,
  "z5_q84" = mstar_fit_data$z5$q84_fit,
  
  "z7_fit" = mstar_fit_data$z7$fit_mstar,
  "z7_q16" = mstar_fit_data$z7$q16_fit,
  "z7_q84" = mstar_fit_data$z7$q84_fit,
  
  "z10_fit" = mstar_fit_data$z10$fit_mstar,
  "z10_q16" = mstar_fit_data$z10$q16_fit,
  "z10_q84" = mstar_fit_data$z10$q84_fit
)
fwrite(
  mstar_mstar_fits,
  "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/mstar_mstar_fits.csv"
)

mstar_z5 = data.frame(
    mstar_fit_data$z5$input_data
)
fwrite(
  mstar_z5,
  "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/mstar_z5.csv"
)

mstar_z7 = data.frame(
  mstar_fit_data$z7$input_data
)
fwrite(
  mstar_z7,
  "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/mstar_z7.csv"
)

mstar_z10 = data.frame(
  mstar_fit_data$z10$input_data
)
fwrite(
  mstar_z10,
  "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/mstar_z10.csv"
)

sfs_fits = data.frame(
  "mstar" = fit_data$withAGN$z5$test_mstar,
  
  "z5_withAGN" = fit_data$withAGN$z5$med_fit,
  "z5_withAGN_q16" = fit_data$withAGN$z5$q16_fit,
  "z5_withAGN_q84" = fit_data$withAGN$z5$q84_fit,
  
  "z7_withAGN" = fit_data$withAGN$z7$med_fit,
  "z7_withAGN_q16" = fit_data$withAGN$z7$q16_fit,
  "z7_withAGN_q84" = fit_data$withAGN$z7$q84_fit,
  
  "z10_withAGN" = fit_data$withAGN$z10$med_fit,
  "z10_withAGN_q16" = fit_data$withAGN$z10$q16_fit,
  "z10_withAGN_q84" = fit_data$withAGN$z10$q84_fit,
  
  "z5_withoutAGN" = fit_data$withoutAGN$z5$med_fit,
  "z5_withoutAGN_q16" = fit_data$withoutAGN$z5$q16_fit,
  "z5_withoutAGN_q84" = fit_data$withoutAGN$z5$q84_fit,
  
  "z7_withoutAGN" = fit_data$withoutAGN$z7$med_fit,
  "z7_withoutAGN_q16" = fit_data$withoutAGN$z7$q16_fit,
  "z7_withoutAGN_q84" = fit_data$withoutAGN$z7$q84_fit,
  
  "z10_withoutAGN" = fit_data$withoutAGN$z10$med_fit,
  "z10_withoutAGN_q16" = fit_data$withoutAGN$z10$q16_fit,
  "z10_withoutAGN_q84" = fit_data$withoutAGN$z10$q84_fit
)
fwrite(
  sfs_fits,
  "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/sfs_fits.csv"
)

z5_sfs_data = data.frame(
  "withAGN_mstar" = fit_data$withAGN$z5$input_data$mstar,
  "withAGN_mstar_err" = fit_data$withAGN$z5$input_data$mstar_err,
  "withAGN_sfr" = fit_data$withAGN$z5$input_data$ydata,
  "withAGN_sfr_err" = fit_data$withAGN$z5$input_data$ydata_err,
  
  "withoutAGN_mstar" = fit_data$withoutAGN$z5$input_data$mstar,
  "withoutAGN_mstar_err" = fit_data$withoutAGN$z5$input_data$mstar_err,
  "withoutAGN_sfr" = fit_data$withoutAGN$z5$input_data$ydata,
  "withoutAGN_sfr_err" = fit_data$withoutAGN$z5$input_data$ydata_err
)
fwrite(
  z5_sfs_data,
  "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/sfs_z5.csv"
)

z7_sfs_data = data.frame(
  "withAGN_mstar" = fit_data$withAGN$z7$input_data$mstar,
  "withAGN_mstar_err" = fit_data$withAGN$z7$input_data$mstar_err,
  "withAGN_sfr" = fit_data$withAGN$z7$input_data$ydata,
  "withAGN_sfr_err" = fit_data$withAGN$z7$input_data$ydata_err,
  
  "withoutAGN_mstar" = fit_data$withoutAGN$z7$input_data$mstar,
  "withoutAGN_mstar_err" = fit_data$withoutAGN$z7$input_data$mstar_err,
  "withoutAGN_sfr" = fit_data$withoutAGN$z7$input_data$ydata,
  "withoutAGN_sfr_err" = fit_data$withoutAGN$z7$input_data$ydata_err
)
fwrite(
  z7_sfs_data,
  "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/sfs_z7.csv"
)

z10_sfs_data = data.frame(
  "withAGN_mstar" = fit_data$withAGN$z10$input_data$mstar,
  "withAGN_mstar_err" = fit_data$withAGN$z10$input_data$mstar_err,
  "withAGN_sfr" = fit_data$withAGN$z10$input_data$ydata,
  "withAGN_sfr_err" = fit_data$withAGN$z10$input_data$ydata_err,
  
  "withoutAGN_mstar" = fit_data$withoutAGN$z10$input_data$mstar,
  "withoutAGN_mstar_err" = fit_data$withoutAGN$z10$input_data$mstar_err,
  "withoutAGN_sfr" = fit_data$withoutAGN$z10$input_data$ydata,
  "withoutAGN_sfr_err" = fit_data$withoutAGN$z10$input_data$ydata_err
)
fwrite(
  z10_sfs_data,
  "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/sfs_z10.csv"
)

smf_data = data.frame(
  "mstar" = fit_data$withAGN$z5$test_mstar,
  "z5_smf" = smf_z5,
  "z6_smf" = smf_z6,
  "z7_smf" = smf_z7,
  "z7_smf" = smf_z7,
  "z9_smf" = smf_z9,
  "z10_smf" = smf_z10
)
#
# smf_data = data.frame(
#   "mstar" = fit_data$withAGN$z5$test_mstar,
#   "z5_smf" = smf_z5_cosmos,
#   "z6_smf" = smf_z6_cosmos,
#   "z7_smf" = smf_z7_cosmos,
#   "z7_smf" = smf_z7,
#   "z9_smf" = smf_z9,
#   "z10_smf" = smf_z10
# )
fwrite(
  smf_data,
  "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/smf_data.csv"
)

# harikane_23_phot = data.frame(
#   z = c(9,12,16),
#   csfh = c(-2.61, -3.23, -3.59) + log10(1/1.53),
#   hi = c(0.18, 0.29, 0.33),
#   lo = c(0.16, 0.27, 2.83)
# )
# fwrite(
#   harikane_23_phot,
#   file = "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/literature/harikane_JWST_2023.csv"
# )
#
# fwrite(
#   bouwens_2015,
#   file = "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/literature/bouwens_2015.csv"
# )
#
# bouwens_2023 = data.frame(
#   redshift = c(8.7, 10.5, 12.6, 14.7),
#   csfh = c(-3.00, -3.82, -3.24, -2.59),
#   hi = c(0.24, 0.30, 0.37, 0.6), #last will be upper limits
#   lo = c(0.24, 0.30, 0.48, 0.0)
# ) #https://arxiv.org/pdf/2211.02607.pdf
# fwrite(
#   bouwens_2023,
#   file = "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/literature/bouwens_2023.csv"
# )

