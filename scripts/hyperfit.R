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

spline_fit = function(zedge_min, zedge_max, data, key, df){
  ## keys are fixed
  
  idx = data$z >= zedge_min & data$z < zedge_max
  data = data[idx, ]
  xdata = c(log10(data[["Mstar50"]]) )
  xerr_floor =  0.0
  xdata_err = c( (data[["Mstar16"]] + data[["Mstar84"]])/(log(10)*data[["Mstar50"]]) )
  xdata_err = sqrt((xerr_floor)^2 + xdata_err^2)
  # xdata_err = rep(0.4, length(xdata))
  
  med_key = paste0(key, "50")
  up_key = paste0(key, "84")
  down_key = paste0(key, "16")
  
  ydata = c( log10(data[[med_key]]) )
  # ydata_err = rep(0.4, length(ydata))
  yerr_floor = 0.0
  ydata_err = c( (data[[down_key]] + data[[up_key]])/(log(10)*data[[med_key]]) )
  ydata_err = sqrt((yerr_floor)^2 + ydata_err^2)

  midx = xdata >= 8

  break_mass = 9.0
  Nsamples = 5001
  spline_data = foreach(i = 1:Nsamples, .combine = "rbind")%dopar%{

    xx = (xdata[midx]-break_mass) + rnorm(length(xdata[midx]), mean = 0, sd = xdata_err[midx])
    yy = ydata[midx] + rnorm(length(ydata[midx]), mean = 0, sd = ydata_err[midx])

    fit = smooth.spline(
      xx, yy, df = df, w = 1/(ydata_err[midx]^2 + xdata_err[midx]^2)
    )
    out = predict(fit, test_mstar)
    lm_fit = lm(out$y~test_mstar)
    c(lm_fit$coefficients[2], lm_fit$coefficients[1])
  }

  med_fit_spline = colQuantiles(matrix(spline_data, nrow = Nsamples, ncol = 2),
                         probs = 0.5)
  q16_fit_spline = colQuantiles(matrix(spline_data, nrow = Nsamples, ncol = 2),
                         probs = 0.16)
  q84_fit_spline = colQuantiles(matrix(spline_data, nrow = Nsamples, ncol = 2),
                         probs = 0.84)

  Niter_max = 500000

  hf = hyper.fit(
    X = cbind(
      xdata[midx]-break_mass, ydata[midx]
    ),
    parm = c(med_fit_spline[1], med_fit_spline[2],0.3),
    vars = cbind(xdata_err[midx],
                 ydata_err[midx])^2,
    algo.func = "LD",
    algo.method = "CHARM",
    itermax = Niter_max,
    Specs = list(alpha.star = 0.44),
    Thinning = 5,
    #prior = function(parm){
      # dnorm(parm[2], mean = 1.0, sd = 1.0, log = T) +
      # dnorm(parm[1], mean = med_fit_spline[1], sd = 1.0, log=T)
    # }
    prior = function(parm){
      dnorm(parm[1], mean = med_fit_spline[1], sd = 0.5, log=T) +
        dnorm(parm[2], mean = med_fit_spline[2], sd = 1.0, log=T)
    }
  )

  # fit_samples = foreach(i = 1:Niter_max, .combine = "rbind")%do%{
  #
  #   xx = rnorm(length(xdata[midx]), xdata[midx], xdata_err[midx])
  #   yy = rnorm(length(ydata[midx]), ydata[midx], ydata_err[midx])
  #   if(i %% 100 == 0){
  #     print(i)
  #   }
  #
  #
  #   hf$parm[1]*test_mstar + hf$parm[2]
  # }

  N_posterior = dim(hf$fit$Posterior1)[1]
  posterior_samples = hf$fit$Posterior1[ceiling((0.68*N_posterior)):(N_posterior), ]
  med_fit_parm = colQuantiles(matrix(posterior_samples, nrow = dim(posterior_samples)[1], ncol = dim(posterior_samples)[2]),
                         probs = 0.5)
  # med_fit_parm = hf$parm
  q16_fit_parm = colQuantiles(matrix(posterior_samples, nrow = dim(posterior_samples)[1], ncol = dim(posterior_samples)[2]),
                         probs = 0.16)
  q84_fit_parm = colQuantiles(matrix(posterior_samples, nrow = dim(posterior_samples)[1], ncol = dim(posterior_samples)[2]),
                         probs = 0.84)

  uncertainty_parm = q84_fit_parm-q16_fit_parm

  fit_samples = foreach(i = 1:Nsamples, .combine = "rbind")%do%{
    if(i %% 1000 == 0){
      print(i)
    }
    samp = rnorm(n = 3,
                 mean = c(
                   med_fit_parm[1],
                   med_fit_parm[2],
                   med_fit_parm[3]
                 ),
                 sd = c(
                   uncertainty_parm[1],
                   uncertainty_parm[2],
                   uncertainty_parm[3]
                 ))
    (samp[1])*(test_mstar-break_mass) + (samp[2])
  }
  med_fit = colQuantiles(matrix(fit_samples, nrow = dim(fit_samples)[1], ncol = dim(fit_samples)[2]),
                         probs = 0.5)
  q16_fit = colQuantiles(matrix(fit_samples, nrow = dim(fit_samples)[1], ncol = dim(fit_samples)[2]),
                         probs = 0.16)
  q84_fit = colQuantiles(matrix(fit_samples, nrow = dim(fit_samples)[1], ncol = dim(fit_samples)[2]),
                         probs = 0.84)

  model = med_fit_parm[1]*(xdata[midx] - break_mass) + med_fit_parm[2]
  resid = (ydata[midx] - model)/ydata_err[midx]
  outliers = list(
    "1sigma" = sum(resid>1)/sum(midx),
    "2sigma" = sum(resid>2)/sum(midx),
    "3sigma" = sum(resid>3)/sum(midx)
  )

  chi_sq = sum(resid^2)
  red_chi_sq = chi_sq/(sum(midx)-2)

  # med_fit = med_fit_parm
  # q16_fit = q16_fit_parm
  # q84_fit = q84_fit_parm
  # med_fit = med_fit_parm[1]*(test_mstar) + med_fit_parm[2]
  # q16_fit = q16_fit_parm[1]*(test_mstar) + q16_fit_parm[2]
  # q84_fit = q84_fit_parm[1]*(test_mstar) + q84_fit_parm[2]

  return(
    list(
      "med_fit" = med_fit,
      "q16_fit" = q16_fit,
      "q84_fit" = q84_fit,
      "med_fit_parm" = med_fit_parm,
      "q16_fit_parm" = q16_fit_parm,
      "q84_fit_parm" = q84_fit_parm,
      "test_mstar" = test_mstar,
      "hf" = hf,
      "resid" = resid,
      "chisq" = chi_sq,
      "red_chi_sq" = red_chi_sq,
      "Posterior.Thinned.BurnedIn" = posterior_samples,
      "init_parm_prior" = med_fit_spline,
      "sd_prior" = q84_fit_spline - q16_fit_spline,
      "outliers" = outliers,
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

saveRDS(fit_data,
        "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/fit_data_master.rds")

pdf("/Users/22252335/Documents/PROJ2_JWST_CSFH/plots/triangle_plots/z5_withAGN.pdf")
magtri(
  fit_data$withAGN$z5$Posterior.Thinned.BurnedIn
)
dev.off()
pdf("/Users/22252335/Documents/PROJ2_JWST_CSFH/plots/triangle_plots/z7_withAGN.pdf")
magtri(
  fit_data$withAGN$z7$Posterior.Thinned.BurnedIn
)
dev.off()
pdf("/Users/22252335/Documents/PROJ2_JWST_CSFH/plots/triangle_plots/z10_withAGN.pdf")
magtri(
  fit_data$withAGN$z10$Posterior.Thinned.BurnedIn
)
dev.off()

pdf("/Users/22252335/Documents/PROJ2_JWST_CSFH/plots/triangle_plots/z5_withoutAGN.pdf")
magtri(
  fit_data$withoutAGN$z5$Posterior.Thinned.BurnedIn
)
dev.off()
pdf("/Users/22252335/Documents/PROJ2_JWST_CSFH/plots/triangle_plots/z7_withoutAGN.pdf")
magtri(
  fit_data$withoutAGN$z7$Posterior.Thinned.BurnedIn
)
dev.off()
pdf("/Users/22252335/Documents/PROJ2_JWST_CSFH/plots/triangle_plots/z10_withoutAGN.pdf")
magtri(
  fit_data$withoutAGN$z10$Posterior.Thinned.BurnedIn
)
dev.off()

plot(density(fit_data$withAGN$z7$resid))
curve(dnorm, add = T)

smf_z5 = dfmodel(fit_data$withAGN$z5$test_mstar, p = c(log10(5.16e-5), 10.97, -1.70)) # https://arxiv.org/pdf/1507.05636.pdf
# smf_z5 = dfmodel(fit_data$withAGN$z5$test_mstar, p = c(log10(0.14e-3), 10.30, -1.70)) # https://arxiv.org/pdf/1507.05636.pdf
# smf_z5 = dfmodel(fit_data$withAGN$z5$test_mstar, p = c(log10(0.14e-3), 10.30, -1.46))
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

fwrite(
  smf_data,
  "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/smf_data.csv"
)
