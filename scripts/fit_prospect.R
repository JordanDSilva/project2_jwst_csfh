library(data.table)
library(ProSpect)
library(dplyr)
library(celestial)
library(Highlander)
library(magicaxis)
library(ProFound)
library(doParallel)
library(foreach)
library(stringr)
library(Rwcs)
library(Rfits)

proj_dir = "/Users/22252335/Documents/PROJ2_JWST_CSFH/"

galaxies = data.frame( fread(paste0(proj_dir, "/data/catalogues/high_z_profound.csv") ))
fincat = galaxies

wave_temp = sapply(c(cenwave$cenwave, EAZY_filters$cenwave), list)
names(wave_temp) = c(cenwave$filter, EAZY_filters$info)

grep("F.+fluxt$", names(galaxies), value=T)

lambda_temp = c(wave_temp$`hst/ACS_update_sep07/wfc_f606w_t81.dat`,
                wave_temp$`hst/ACS_update_sep07/wfc_f814w_t81.dat`,
                wave_temp$`hst/wfc3/IR/f125w.dat`,
                wave_temp$`hst/wfc3/IR/f140w.dat`,
                wave_temp$`hst/wfc3/IR/f160w.dat`,
                wave_temp$F115W_JWST,
                wave_temp$F150W_JWST,
                wave_temp$F200W_JWST,
                wave_temp$F277W_JWST,
                wave_temp$F356W_JWST,
                wave_temp$F410M_JWST,
                wave_temp$F444W_JWST,
                wave_temp$`hst/wfc3/IR/f105w.dat`,
                wave_temp$`hst/ACS_update_sep07/wfc_f775w_t81.dat`,
                wave_temp$F090W_JWST,
                wave_temp$`hst/ACS_update_sep07/wfc_f435w_t81.dat`,
                wave_temp$`hst/ACS_update_sep07/wfc_f475w_t81.dat`
                
)
lambda = sort(lambda_temp)
filter_names_temp = sapply(c(cenwave$filter, EAZY_filters$info), list)
filter_temp = c(filter_names_temp$`hst/ACS_update_sep07/wfc_f606w_t81.dat`,
                filter_names_temp$`hst/ACS_update_sep07/wfc_f814w_t81.dat`,
                filter_names_temp$`hst/wfc3/IR/f125w.dat`,
                filter_names_temp$`hst/wfc3/IR/f140w.dat`,
                filter_names_temp$`hst/wfc3/IR/f160w.dat`,
                filter_names_temp$F115W_JWST,
                filter_names_temp$F150W_JWST,
                filter_names_temp$F200W_JWST,
                filter_names_temp$F277W_JWST,
                filter_names_temp$F356W_JWST,
                filter_names_temp$F410M_JWST,
                filter_names_temp$F444W_JWST,
                filter_names_temp$`hst/wfc3/IR/f105w.dat`,
                filter_names_temp$`hst/ACS_update_sep07/wfc_f775w_t81.dat`,
                filter_names_temp$F090W_JWST,
                filter_names_temp$`hst/ACS_update_sep07/wfc_f435w_t81.dat`,
                filter_names_temp$`hst/ACS_update_sep07/wfc_f475w_t81.dat`
)
filters = filter_temp[order(lambda_temp)]
filtout = sapply(filters, function(x){approxfun(getfilt(x))})

correct_filters = !grepl("F127M|F139M|F153M", names(fincat))

fluxc_names = grep("F.+fluxc$", names(fincat)[correct_filters], value=T)[order(lambda_temp)]
fluxc_err_names = grep("F.+scaled_fluxc_err$", names(fincat)[correct_filters], value=T)[order(lambda_temp)]

fluxt_names = grep("F.+fluxt$", names(fincat)[correct_filters], value=T)[order(lambda_temp)]
fluxt_err_names = grep("F.+scaled_fluxt_err$", names(fincat)[correct_filters], value=T)[order(lambda_temp)]

flux_names_corrt = paste0(fluxt_names, "_corr")
flux_err_names_corrt = paste0(fluxt_err_names, "_corr")

flux_names_corrc = paste0(fluxc_names, "_corr")
flux_err_names_corrc = paste0(fluxc_err_names, "_corr")

registerDoParallel(cores = 7)

do_fit = function(){
 message("Running with AGN")
withAGN=foreach(II = 1:dim(galaxies)[1])%dopar%{
# withAGN=foreach(II = c(1))%dopar%{

  print(paste0("Running: gal_", II))
  Redshift = galaxies$z[II]

  flux_data = list()
  flux_data$filter = filters
  flux_data$cenwave = lambda
  flux_data$flux = unlist(galaxies[II, flux_names_corrt])
  flux_data$fluxerr = unlist(galaxies[II, flux_err_names_corrt])
  flux_data = data.frame(flux_data)

  pfunc=function(parm){
    dnorm(parm['tau_birth'],mean=-0.2,sd=0.5,log=TRUE)+
      (-20*erf(parm['tau_screen']-2))+
      dnorm(parm['alpha_SF_birth'],mean=2,sd=1,log=TRUE)+
      dnorm(parm['alpha_SF_screen'],mean=2,sd=1,log=TRUE)+
      (100*erf(parm['mperiod']+2) - 100)
  }

  #Setting mpeak limitis
  LookbackTime = cosdistTravelTime(Redshift, ref = '737')*1e9
  agemax = cosdistUniAgeNow(ref="737")*1e9 - LookbackTime
  upperlimit = (agemax/1e9)
  lowerlimit = (0 - LookbackTime)/1e9

  Data=list(flux=flux_data,
            arglist=list(
              z=Redshift,
              massfunc=massfunc_snorm_trunc,
              Z=Zfunc_massmap_lin,
              LumDist_Mpc = cosdistLumDist(z = Redshift, ref = '737'),
              agemax=agemax,
              Zagemax = (agemax/1e9),
              Zstart = 1e-6,
              magemax=(agemax/1e9), #, yield = 0.02
              emission = T,
              AGNrm = 60,
              AGNbe = -0.5,
              AGNal = 4.0,
              AGNct = 100,
              IGMabsorb = pnorm(Redshift, mean = 3.8, sd = 1.2),
              ref="737",
              H0 = 70,
              OmegaM = 0.3
            ),
            speclib=BC03hr,
            Dale=Dale_NormTot,
            filtout=filtout,
            SFH=SFHfunc,
            verbose=FALSE,
            AGN = Fritz,
            Dale_M2L_func=Dale_M2L_func,
            parm.names=c('mSFR','mpeak','mperiod','mskew','tau_birth','tau_screen','alpha_SF_birth','alpha_SF_screen','Zfinal','AGNlum', 'AGNan', 'AGNta'),
            mon.names=c("LP","masstot","dustmass.birth", "dustmass.screen", "dustmass.total", "dustlum.birth", "dustlum.screen",
                        "dustlum.total", "SFRburst", paste("flux.",filters,sep='')),
            logged=c(T,F,T,F,T,T,F,F,T,T, F, T),
            intervals=list(lo=c(-3,-2, log10(0.3),-0.5,-2.5,-5,0,0,-4, 35, 0.001, -1),
                           hi=c(4,upperlimit,2,1,1.5,1,4,4,-1.3,49, 89.990, 1)),
            fit = 'LD',
            N=length(filters),
            prior=pfunc
  )

  Data$flux$cenwave = flux_data$cenwave

  startpoint = (Data$intervals$lo+Data$intervals$hi)/2
  startpoint[10] = 44
  testHigh = Highlander(startpoint, Data,
                        ProSpectSEDlike, Niters=c(2000,2000),  NfinalMCMC = 2000, lower=Data$intervals$lo,
                        upper=Data$intervals$hi, seed=666, optim_iters = 2, likefunctype = 'LD'
  )

  Data$fit = 'check'
  bestfit=ProSpectSEDlike(testHigh$par, Data=Data)

  names(testHigh$parm) = names(testHigh$LD_last$Summary1[1:length(startpoint), 1])
  par_list = list(
    mean_par = testHigh$LD_last$Summary1[1:length(startpoint), 1],
    sd_par = testHigh$LD_last$Summary1[1:length(startpoint), 2],
    best_par = testHigh$parm
  )

  sfr_samples = {}
  mstar_samples = {}
  agn_samples = {}
  Nerr_samples = 201

  for(i in 1:Nerr_samples){
    mean_par = testHigh$LD_last$Summary1[1:length(startpoint), 1]
    sd_par = testHigh$LD_last$Summary1[1:length(startpoint), 2]
    sample_par = rnorm(
      n = length(startpoint),
      mean = mean_par,
      sd = sd_par
    )
    sample_fit=ProSpectSEDlike(sample_par, Data=Data)

    star_func = SMstarfunc(
      massfunc = sample_fit$Data$arglist$massfunc,
      speclib = sample_fit$Data$speclib,
      stellpop = "BC03hr",
      ref = "737",
      z = sample_fit$Data$arglist$z,
      Z = sample_fit$Data$arglist$Z,
      Zfinal = (10^sample_fit$parm["Zfinal"]),

      #mass_func_snorm_params
      mSFR = 10^(sample_fit$parm["mSFR"]),
      mpeak = (sample_fit$parm["mpeak"]),
      mperiod = 10^(sample_fit$parm["mperiod"]),
      mskew = (sample_fit$parm["mskew"]),
      magemax = sample_fit$Data$arglist$agemax/1e9
    )

    sfr_samples = c(sample_fit$SEDout$Stars$SFRburst, sfr_samples)
    mstar_samples = c(star_func["TotSMstar"], mstar_samples)
    agn_samples = c(10^sample_fit$parm["AGNlum"], agn_samples)
  }

  sfr_med = quantile(sfr_samples, 0.5, na.rm = T)
  sfr_16 = sfr_med - quantile(sfr_samples, 0.16, na.rm = T)
  sfr_84 = quantile(sfr_samples, 0.84, na.rm = T) - sfr_med

  mstar_med = quantile(mstar_samples, 0.5, na.rm = T)
  mstar_16 = mstar_med - quantile(mstar_samples, 0.16, na.rm = T)
  mstar_84 = quantile(mstar_samples, 0.84, na.rm = T) - mstar_med

  agn_med = quantile(agn_samples, 0.5, na.rm = T)
  agn_16 = agn_med - quantile(agn_samples, 0.16, na.rm = T)
  agn_84 = quantile(agn_samples, 0.84, na.rm = T) - agn_med


  quantity_table = data.frame(
    "SFRBurst50" = sfr_med,
    "SFRBurst16" = sfr_16,
    "SFRBurst84" = sfr_84,
    "AGNlum50" = agn_med,
    "AGNlum16" = agn_16,
    "AGNlum84" = agn_84,
    "Mstar50" = mstar_med,
    "Mstar16" = mstar_16,
    "Mstar84" = mstar_84
  )

  saveRDS(
    list(
      "highlander" = testHigh,
      "bestfit"=bestfit,
      "quantities" = quantity_table,
      "pars" = par_list
    ),
    paste0(proj_dir, "/ProSpectOut/withAGN/gal_", II, "_", galaxies$VID[II], "_", galaxies$MODULE[II], ".rds")
  )

  gc()
} ## Do the fitting here

load_with_agn = foreach(II = 1:dim(galaxies)[1], .combine = bind_rows)%do%{
  print(II)
  foo = readRDS(paste0(proj_dir, "/ProSpectOut/withAGN/gal_", II, "_", galaxies$VID[II], "_", galaxies$MODULE[II], ".rds"))
  testHigh = foo$highlander
  startpoint = foo$bestfit$parm

  par_list = list(
    mean_par = testHigh$LD_last$Summary1[1:length(startpoint), 1],
    sd_par = testHigh$LD_last$Summary1[1:length(startpoint), 3], #Monte Carlo Standard Error
    best_par = testHigh$parm
  )

  sfr_samples = {}
  mstar_samples = {}
  agn_samples = {}
  Nerr_samples = 201

  pdf(paste0(proj_dir, "/plots/ProSpect/withAGN/gal_", II, "_", galaxies$VID[II], "_", galaxies$MODULE[II], ".pdf"),
      width = 10, height = 10)
  plot(foo$bestfit, xlim = 10^c(3.0, 5.0), ylim = 10^c(-10, -3))
  plot(foo$bestfit$SEDout, xlim = 10^c(3.0, 5.0), ylim = 10^c(3, 9))
  par(mfrow = c(1,1), mar = rep(1.0,4), oma = rep(1.0,4))
  magplot(NA, xlim = 10^c(3.0, 5.0), ylim = 10^c(-10, -3), log = "xy")
  legend(x='topright', legend = paste0("z=",galaxies$z[II]))
  for(ii in 1:Nerr_samples){
    sample_par = rnorm(
      n = length(startpoint),
      mean = par_list$best_par,
      sd = par_list$sd_par
    )
    sample_fit=ProSpectSEDlike(sample_par, Data=foo$bestfit$Data)

    lines(sample_fit$SEDout$FinalFlux$wave, sample_fit$SEDout$FinalFlux$flux, col = "grey")

    star_func = SMstarfunc(
      massfunc = sample_fit$Data$arglist$massfunc,
      speclib = sample_fit$Data$speclib,
      stellpop = "BC03hr",
      ref = "737",
      z = sample_fit$Data$arglist$z,
      Z = sample_fit$Data$arglist$Z,
      Zfinal = (10^sample_fit$parm["Zfinal"]),

      #mass_func_snorm_params
      mSFR = 10^(sample_fit$parm["mSFR"]),
      mpeak = (sample_fit$parm["mpeak"]),
      mperiod = 10^(sample_fit$parm["mperiod"]),
      mskew = (sample_fit$parm["mskew"]),
      magemax = sample_fit$Data$arglist$agemax/1e9
    )

    sfr_samples = c(sample_fit$SEDout$Stars$SFRburst, sfr_samples)
    mstar_samples = c(star_func["TotSMstar"], mstar_samples)
    agn_samples = c(10^sample_fit$parm["AGNlum"], agn_samples)
  }

  lines(foo$bestfit$SEDout$FinalFlux, log = "xy", type = "l", xlim = 10^c(3, 5.5), ylim = 10^c(-8, -6), lwd=5)
  points(foo$bestfit$Data$flux$cenwave, foo$bestfit$Data$flux$flux, pch = 16, cex = 2, col="red")
  magerr(foo$bestfit$Data$flux$cenwave, foo$bestfit$Data$flux$flux, ylo = foo$bestfit$Data$flux$fluxerr, col="red")
  dev.off()

  sfr_med = quantile(sfr_samples, 0.5, na.rm = T)
  sfr_16 = sfr_med - quantile(sfr_samples, 0.16, na.rm = T)
  sfr_84 = quantile(sfr_samples, 0.84, na.rm = T) - sfr_med

  mstar_med = quantile(mstar_samples, 0.5, na.rm = T)
  mstar_16 = mstar_med - quantile(mstar_samples, 0.16, na.rm = T)
  mstar_84 = quantile(mstar_samples, 0.84, na.rm = T) - mstar_med

  agn_med = quantile(agn_samples, 0.5, na.rm = T)
  agn_16 = agn_med - quantile(agn_samples, 0.16, na.rm = T)
  agn_84 = quantile(agn_samples, 0.84, na.rm = T) - agn_med

  quan = data.frame(
    "SFRBurst50" = sfr_med,
    "SFRBurst16" = sfr_16,
    "SFRBurst84" = sfr_84,
    "AGNlum50" = agn_med,
    "AGNlum16" = agn_16,
    "AGNlum84" = agn_84,
    "Mstar50" = mstar_med,
    "Mstar16" = mstar_16,
    "Mstar84" = mstar_84,
    "RA" = galaxies$RAmax[II],
    "DEC" = galaxies$Decmax[II],
    "VISITID" = paste0(galaxies$VID[II]),
    "MODULE" = galaxies$MODULE[II],
    "COLID" = II
  )

  quan$z = foo$bestfit$Data$arglist$z
  quan$LP = foo$bestfit$LP

  return(quan)
}
fwrite(load_with_agn, paste0(proj_dir, "/data/catalogues/ProSpect_highz_withAGN.csv") )

message("Running without AGN")

withoutAGN=foreach(II = 1:dim(galaxies)[1])%dopar%{
  message(paste0("Running: gal_", II))
  Redshift = galaxies$z[II]

  flux_data = list()
  flux_data$filter = filters
  flux_data$cenwave = lambda
  flux_data$flux = unlist(galaxies[II, flux_names_corrt])
  flux_data$fluxerr = unlist(galaxies[II, flux_err_names_corrt])
  flux_data = data.frame(flux_data)

  pfunc=function(parm){
    dnorm(parm['tau_birth'],mean=-0.2,sd=0.5,log=TRUE)+
      (-20*erf(parm['tau_screen']-2))+
      dnorm(parm['alpha_SF_birth'],mean=2,sd=1,log=TRUE)+
      dnorm(parm['alpha_SF_screen'],mean=2,sd=1,log=TRUE)+
      (100*erf(parm['mperiod']+2) - 100)
  }

  #Setting mpeak limitis
  LookbackTime = cosdistTravelTime(Redshift, ref = '737')*1e9
  agemax = cosdistUniAgeNow(ref="737")*1e9 - LookbackTime
  upperlimit = (agemax/1e9)
  lowerlimit = (0 - LookbackTime)/1e9

  Data=list(flux=flux_data,
            arglist=list(
              z=Redshift,
              massfunc=massfunc_snorm_trunc,
              Z=Zfunc_massmap_lin,
              LumDist_Mpc = cosdistLumDist(z = Redshift, ref = '737'),
              agemax=agemax,
              Zagemax = (agemax/1e9),
              Zstart = 1e-6,
              magemax=(agemax/1e9), #, yield = 0.02
              emission = T,
              IGMabsorb = pnorm(Redshift, mean = 3.8, sd = 1.2),
              ref="737",
              H0 = 70,
              OmegaM = 0.3
            ),
            speclib=BC03hr,
            Dale=Dale_NormTot,
            filtout=filtout,
            SFH=SFHfunc,
            verbose=FALSE,
            # AGN = Fritz,
            Dale_M2L_func=Dale_M2L_func,
            parm.names=c('mSFR','mpeak','mperiod','mskew','tau_birth','tau_screen','alpha_SF_birth','alpha_SF_screen','Zfinal'),
            mon.names=c("LP","masstot","dustmass.birth", "dustmass.screen", "dustmass.total", "dustlum.birth", "dustlum.screen",
                        "dustlum.total", "SFRburst", paste("flux.",filters,sep='')),
            logged=c(T,F,T,F,T,T,F,F,T),
            intervals=list(lo=c(-3,-2, log10(0.3),-0.5,-2.5,-5,0,0,-4),
                           hi=c(4,upperlimit,2,1,1.5,1,4,4,-1.3)),
            fit = 'LD',
            N=length(filters),
            prior=pfunc
  )

  Data$flux$cenwave = flux_data$cenwave

  startpoint = (Data$intervals$lo+Data$intervals$hi)/2
  testHigh = Highlander(startpoint, Data,
                        ProSpectSEDlike, Niters=c(2000,2000),  NfinalMCMC = 2000, lower=Data$intervals$lo,
                        upper=Data$intervals$hi, seed=666, optim_iters = 2, likefunctype = 'LD')

  Data$fit = 'check'
  bestfit=ProSpectSEDlike(testHigh$par, Data=Data)

  names(testHigh$parm) = names(testHigh$LD_last$Summary1[1:length(startpoint), 1])
  par_list = list(
    mean_par = testHigh$LD_last$Summary1[1:length(startpoint), 1],
    sd_par = testHigh$LD_last$Summary1[1:length(startpoint), 2],
    best_par = testHigh$parm
  )

  sfr_samples = {}
  mstar_samples = {}
  Nerr_samples = 201

  for(i in 1:Nerr_samples){
    mean_par = testHigh$LD_last$Summary1[1:length(startpoint), 1]
    sd_par = testHigh$LD_last$Summary1[1:length(startpoint), 2]
    sample_par = rnorm(
      n = length(startpoint),
      mean = mean_par,
      sd = sd_par
    )
    sample_fit=ProSpectSEDlike(sample_par, Data=Data)

    star_func = SMstarfunc(
      massfunc = sample_fit$Data$arglist$massfunc,
      speclib = sample_fit$Data$speclib,
      stellpop = "BC03hr",
      ref = "737",
      z = sample_fit$Data$arglist$z,
      Z = sample_fit$Data$arglist$Z,
      Zfinal = (10^sample_fit$parm["Zfinal"]),

      #mass_func_snorm_params
      mSFR = 10^(sample_fit$parm["mSFR"]),
      mpeak = (sample_fit$parm["mpeak"]),
      mperiod = 10^(sample_fit$parm["mperiod"]),
      mskew = (sample_fit$parm["mskew"]),
      magemax = sample_fit$Data$arglist$agemax/1e9
    )

    sfr_samples = c(sample_fit$SEDout$Stars$SFRburst, sfr_samples)
    mstar_samples = c(star_func["TotSMstar"], mstar_samples)
  }

  sfr_med = quantile(sfr_samples, 0.5, na.rm = T)
  sfr_16 = sfr_med - quantile(sfr_samples, 0.16, na.rm = T)
  sfr_84 = quantile(sfr_samples, 0.84, na.rm = T) - sfr_med

  mstar_med = quantile(mstar_samples, 0.5, na.rm = T)
  mstar_16 = mstar_med - quantile(mstar_samples, 0.16, na.rm = T)
  mstar_84 = quantile(mstar_samples, 0.84, na.rm = T) - mstar_med

  quantity_table = data.frame(
    "SFRBurst50" = sfr_med,
    "SFRBurst16" = sfr_16,
    "SFRBurst84" = sfr_84,
    "Mstar50" = mstar_med,
    "Mstar16" = mstar_16,
    "Mstar84" = mstar_84
  )

  saveRDS(
    list(
      "highlander" = testHigh,
      "bestfit"=bestfit,
      "quantities" = quantity_table,
      "pars" = par_list
    ),
    paste0(proj_dir, "/ProSpectOut/withoutAGN/gal_", II, "_", galaxies$VID[II], "_", galaxies$MODULE[II], ".rds")
  )

  gc()
} ## Do the fitting here

load_without_agn = foreach(II = 1:dim(galaxies)[1], .combine = bind_rows)%do%{
  print(II)
  foo = readRDS(paste0(proj_dir, "/ProSpectOut/withoutAGN/gal_", II, "_", galaxies$VID[II], "_", galaxies$MODULE[II], ".rds"))
  testHigh = foo$highlander
  startpoint = foo$bestfit$parm

  par_list = list(
    mean_par = testHigh$LD_last$Summary1[1:length(startpoint), 1],
    sd_par = testHigh$LD_last$Summary1[1:length(startpoint), 3], #Monte Carlo Standard Error
    best_par = testHigh$parm
  )

  sfr_samples = {}
  mstar_samples = {}
  agn_samples = {}
  Nerr_samples = 201

  pdf(paste0(proj_dir, "/plots/ProSpect/withoutAGN/gal_", II, "_", galaxies$VID[II], "_", galaxies$MODULE[II], ".pdf"),
      width = 10, height = 10)
  plot(foo$bestfit, xlim = 10^c(3.0, 5.0), ylim = 10^c(-10, -3))
  plot(foo$bestfit$SEDout, xlim = 10^c(3.0, 5.0), ylim = 10^c(3, 9))
  par(mfrow = c(1,1), mar = rep(1.0,4), oma = rep(1.0,4))
  magplot(NA, xlim = 10^c(3.0, 5.0), ylim = 10^c(-10, -3), log = "xy")
  legend(x='topright', legend = paste0("z=",galaxies$z[II]))
  for(ii in 1:Nerr_samples){
    sample_par = rnorm(
      n = length(startpoint),
      mean = par_list$best_par,
      sd = par_list$sd_par
    )
    sample_fit=ProSpectSEDlike(sample_par, Data=foo$bestfit$Data)

    lines(sample_fit$SEDout$FinalFlux$wave, sample_fit$SEDout$FinalFlux$flux, col = "grey")

    star_func = SMstarfunc(
      massfunc = sample_fit$Data$arglist$massfunc,
      speclib = sample_fit$Data$speclib,
      stellpop = "BC03hr",
      ref = "737",
      z = sample_fit$Data$arglist$z,
      Z = sample_fit$Data$arglist$Z,
      Zfinal = (10^sample_fit$parm["Zfinal"]),

      #mass_func_snorm_params
      mSFR = 10^(sample_fit$parm["mSFR"]),
      mpeak = (sample_fit$parm["mpeak"]),
      mperiod = 10^(sample_fit$parm["mperiod"]),
      mskew = (sample_fit$parm["mskew"]),
      magemax = sample_fit$Data$arglist$agemax/1e9
    )

    sfr_samples = c(sample_fit$SEDout$Stars$SFRburst, sfr_samples)
    mstar_samples = c(star_func["TotSMstar"], mstar_samples)
  }

  lines(foo$bestfit$SEDout$FinalFlux, log = "xy", type = "l", xlim = 10^c(3, 5.5), ylim = 10^c(-10, -4), lwd=5)
  points(foo$bestfit$Data$flux$cenwave, foo$bestfit$Data$flux$flux, pch = 16, cex = 2, col="red")
  magerr(foo$bestfit$Data$flux$cenwave, foo$bestfit$Data$flux$flux, ylo = foo$bestfit$Data$flux$fluxerr, col="red")
  dev.off()

  sfr_med = quantile(sfr_samples, 0.5, na.rm = T)
  sfr_16 = sfr_med - quantile(sfr_samples, 0.16, na.rm = T)
  sfr_84 = quantile(sfr_samples, 0.84, na.rm = T) - sfr_med

  mstar_med = quantile(mstar_samples, 0.5, na.rm = T)
  mstar_16 = mstar_med - quantile(mstar_samples, 0.16, na.rm = T)
  mstar_84 = quantile(mstar_samples, 0.84, na.rm = T) - mstar_med


  quan = data.frame(
    "SFRBurst50" = sfr_med,
    "SFRBurst16" = sfr_16,
    "SFRBurst84" = sfr_84,
    "Mstar50" = mstar_med,
    "Mstar16" = mstar_16,
    "Mstar84" = mstar_84,
    "RA" = galaxies$RAmax[II],
    "DEC" = galaxies$Decmax[II],
    "VISITID" = paste0(galaxies$VID[II]),
    "MODULE" = galaxies$MODULE[II],
    "COLID" = II
  )

  quan$z = foo$bestfit$Data$arglist$z
  quan$LP = foo$bestfit$LP

  return(quan)
}
fwrite(load_without_agn, paste0(proj_dir, "/data/ProSpect_highz_withoutAGN.csv") )

message("Fin")
stopImplicitCluster()

}

