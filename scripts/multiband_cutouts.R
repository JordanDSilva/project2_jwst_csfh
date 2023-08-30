library(data.table)
library(foreach)
library(doParallel)
library(Cairo)
library(ProSpect)
library(ProFound)
library(Rwcs)
library(pracma)
library(scales)
library(dplyr)
library(matrixStats)

registerDoParallel(cores = 8)

galaxies = data.frame( fread("/Users/22252335/Documents/PROJ2_JWST_CSFH/data/catalogues/high_z_profound.csv") )
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
                wave_temp$`hst/ACS_update_sep07/wfc_f475_t81.dat`
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


prospect_SED_plot_widget = function(prospect_data){

  wavelength = ProSpectSEDlike(prospect_data$pars$best_par,
                                Data=prospect_data$bestfit$Data)$SEDout$FinalFlux$wave
  best_sample_fit=ProSpectSEDlike(prospect_data$pars$best_par,
                                Data=prospect_data$bestfit$Data)$SEDout$FinalFlux$flux


  # sed_samples = foreach(ii = 1:200, .combine = "rbind")%do%{
  #     sample_par = rnorm(
  #       n = length(prospect_data$pars$best_par),
  #       mean = prospect_data$pars$best_par,
  #       sd = prospect_data$pars$sd_par
  #     )
  #     sample_fit=ProSpectSEDlike(sample_par, Data=prospect_data$bestfit$Data)$SEDout$FinalFlux$flux
  # }

  # hi_sample_fit = apply(sed_samples, 2, function(x)quantile(x, 0.84))
  # lo_sample_fit = apply(sed_samples, 2, function(x)quantile(x, 0.16))
  return(
    list(
      "wave" = wavelength,
      "best" = best_sample_fit
      # "hi" = hi_sample_fit,
      # "lo" = lo_sample_fit
    )
  )


}


# foreach(ii = 1:dim(galaxies)[1])%dopar%{
#   files = list.files(
#     paste0("/Volumes/RAIDY/JWST/ProFound/Measurements/",
#            galaxies$VID[ii], "/", galaxies$MODULE[ii], "/"),
#     full.names = T,
#     recursive = T,
#     pattern = ".rds"
#   )
#
#   idx = unlist(sapply(substr(fluxt_names, 1,5), function(x){which(grepl(x, files))}))
#   multiband_segim_images = lapply(files[idx],
#                                   function(x){
#                                     foo = readRDS(x)
#                                     list(image = foo$image, segim = foo$segim)}
#   )
#
#   xcen = galaxies$xmax[ii]
#   ycen = galaxies$ymax[ii]
#   RA = galaxies$RAcen[ii]
#   Dec = galaxies$Deccen[ii]
#   CairoPNG(paste0("/Users/22252335/Documents/PROJ2_JWST_CSFH/plots/segimmaps/", ii, "_", galaxies$VID[ii], "_",
#                   galaxies$MODULE[ii], ".png"),
#            width = 10, height = 10, units = "in", res = 240)
#   par(mfrow = c(5,5), mar = rep(0,4), oma = rep(0,4))
#
#   box = 50
#
#   for(i in 1:length(idx)){
#     imcut = magcutout(multiband_segim_images[[i]]$image, loc = c(xcen,ycen), box = box)$image
#     if(sum(is.na(imcut)) >= box^2){
#       imcut = matrix(0, box, box)
#     }
#     profoundSegimPlot(image = imcut,
#                       segim = magcutout(multiband_segim_images[[i]]$segim,
#                                         loc = c(xcen,ycen),
#                                         box = box)$image, flip = T)
#     legend(x = "topright",names(idx)[i], col="red", cex = 1.0)
#     # text(x = box/2, y = box/4, names(idx)[i], col="red", cex = 1.3
#   }
#
#   dev.off()
#   print("end")
#
# }

for(k in 1:dim(galaxies)[1]){

  print(k)

  R = Rfits_read_image(paste0("/Volumes/RAIDY/JWST/Patch_Stacks/stack_patch_",galaxies$VID[k], "_F444W_", galaxies$MODULE[k], "_long.fits"))
  G = Rfits_read_image(paste0("/Volumes/RAIDY/JWST/Patch_Stacks/stack_patch_",galaxies$VID[k], "_F277W_", galaxies$MODULE[k], "_long.fits"))
  B = Rfits_read_image(paste0("/Volumes/RAIDY/JWST/Patch_Stacks/stack_patch_",galaxies$VID[k], "_F150W_", galaxies$MODULE[k], "_long.fits"))
  
  detect_image = Rfits_read_image(paste0("/Volumes/RAIDY/JWST/ProFound/Detects/",galaxies$VID[k], "/", galaxies$MODULE[k], "/", 
                                         galaxies$VID[k], "_", galaxies$MODULE[k], "_profound_stack.fits"))

  with_agn_sed = readRDS(
    paste0("~/Documents/PROJ2_JWST_CSFH/ProSpectOut/withAGN/gal_", k, "_", galaxies$VID[k], "_", galaxies$MODULE[k], ".rds")
  )
  
  without_agn_sed = readRDS(
    paste0("~/Documents/PROJ2_JWST_CSFH/ProSpectOut/withoutAGN/gal_", k, "_", galaxies$VID[k], "_", galaxies$MODULE[k], ".rds")
  )

  withAGN_data = prospect_SED_plot_widget(with_agn_sed)
  withoutAGN_data = prospect_SED_plot_widget(without_agn_sed)
  
  det_segim = readRDS(
    paste0("/Volumes/RAIDY/JWST/ProFound/Detects/", galaxies$VID[k], "/",
           galaxies$MODULE[k], "/", galaxies$VID[k], "_", galaxies$MODULE[k], "_profound.rds")
  )$segim

  segim_wcs = Rfits_create_image(
    data = det_segim,
    keyvalues = detect_image$keyvalues
  )

  segim_cutout = segim_wcs[galaxies$RAmax[k], galaxies$Decmax[k], box = 100, type="coord"]$imDat

  png(paste0("~/Documents/PROJ2_JWST_CSFH/plots/RGB_cutouts/gal_", k, "_", galaxies$VID[k], "_", galaxies$MODULE[k], ".png"), height = 500, width = 1000)
  
  par(mfrow = c(1,3), mar = rep(2.0,4), oma = rep(0,4))
  Rwcs_imageRGB(
    R = R[galaxies$RAmax[k], galaxies$Decmax[k], box = 100, type="coord"],
    G = G[galaxies$RAmax[k], galaxies$Decmax[k], box = 100, type="coord"],
    B = B[galaxies$RAmax[k], galaxies$Decmax[k], box = 100, type="coord"],
    
    sparse = 1,
    locut = 5e-4, 
    hicut = 5e-3,
    type="num",
    decorate = F
  )
  profoundSegimPlot(
     NULL,
    segim = segim_cutout,
    add = T
  )
  legend(x = "topleft", "F150W+F277W+F444W")
  plot(detect_image[galaxies$RAmax[k], galaxies$Decmax[k], box = 100, type="coord"],
       flip=T)
  profoundSegimPlot(
    NULL,
    segim = segim_cutout,
    add = T
  )
  magplot(
    NA,
    # log="y",
    ylim = c(-8.3, -6),
    xlim = c(0, 5)
  )
  legend(x = "topleft", legend = paste0("Redshift=", galaxies$z[k]))
  lines(
    withAGN_data$wave/10000, log10(withAGN_data$best), lwd = 2, col="red"
  )
  lines(
    withoutAGN_data$wave/10000,
    log10(withoutAGN_data$best),
    lwd = 2, col="cornflowerblue"
  )

  for(filtout_func in with_agn_sed$bestfit$Data$filtout){

    wave = seq(0,6,0.01)*1e4

    transmission_func = (
      filtout_func(
        wave
      )
    )
    transmission_func = (transmission_func/max(transmission_func, na.rm = T)) - 8.5

    # magplot(wave/1e4, transmission_func)

    polygon(
      c( wave/1e4, rev(wave/1e4) ),
      c( (transmission_func), rev(transmission_func) ),
      col = alpha("grey", 0.2),
      border = NA
    )
  }

  points(
    unlist(lambda)/1e4,
    log10(unlist(galaxies[k, flux_names_corrt])),
    pch = 16
  )
  legend(x = "topleft", legend = paste0("Redshift=", galaxies$z[k]))
  magerr(
    unlist(lambda)/1e4,
    log10(unlist(galaxies[k, flux_names_corrt])),
    ylo = unlist(galaxies[k, flux_err_names_corrt]) / (log(10)*unlist(galaxies[k, flux_names_corrt]))
  )
  dev.off()
}



# k = 106
# detect_image = Rfits_read_image(paste0("/Volumes/RAIDY/JWST/ProFound/Detects/",galaxies$VID[k], "/", galaxies$MODULE[k], "/",
#                                          galaxies$VID[k], "_", galaxies$MODULE[k], "_profound_stack.fits"))
# det_segim = readRDS(
#     paste0("/Volumes/RAIDY/JWST/ProFound/Detects/", galaxies$VID[k], "/", galaxies$MODULE[k], "/", galaxies$VID[k], "_", galaxies$MODULE[k], "_profound.rds")
#   )$segim
# segim_wcs = Rfits_create_image(
#   data = det_segim,
#   keyvalues = detect_image$keyvalues
# )
# cutout = detect_image[galaxies$RAmax[k], galaxies$Decmax[k], box = 100, type="coord"]
# segim_cutout = segim_wcs[galaxies$RAmax[k], galaxies$Decmax[k], box = 100, type="coord"]$imDat
# plot(cutout,
#        flip=T)
# profoundSegimPlot(
#     NULL,
#     segim = segim_cutout,
#     add = T
#   )

