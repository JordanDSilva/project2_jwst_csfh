library(data.table)
library(foreach)
library(doParallel)
library(Cairo)
library(ProSpect)
library(ProFound)
library(Rwcs)

registerDoParallel(cores = 8)

galaxies = data.frame( fread("/Users/22252335/Documents/PROJ2_JWST_CSFH/data/catalogues/high_z_profound.csv") )

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

foreach(ii = 1:dim(galaxies)[1])%dopar%{
  files = list.files(
    paste0("/Volumes/RAIDY/JWST/ProFound/Measurements/",
           galaxies$VID[ii], "/", galaxies$MODULE[ii], "/"),
    full.names = T,
    recursive = T,
    pattern = ".rds"
  )
  
  idx = unlist(sapply(substr(fluxt_names, 1,5), function(x){which(grepl(x, files))}))
  multiband_segim_images = lapply(files[idx],
                                  function(x){
                                    foo = readRDS(x)
                                    list(image = foo$image, segim = foo$segim)}
  )
  
  xcen = galaxies$xcen[ii]
  ycen = galaxies$ycen[ii]
  RA = galaxies$RAcen[ii]
  Dec = galaxies$Deccen[ii]
  CairoPNG(paste0("/Users/22252335/Documents/PROJ2_JWST_CSFH/plots/segimmaps/", ii, "_", galaxies$VID[ii], "_",
                  galaxies$MODULE[ii], ".png"),
           width = 10, height = 10, units = "in", res = 240)
  par(mfrow = c(5,5), mar = rep(0,4), oma = rep(0,4))
  
  box = 50
  
  for(i in 1:length(idx)){
    imcut = magcutout(multiband_segim_images[[i]]$image, loc = c(xcen,ycen), box = box)$image
    if(sum(is.na(imcut)) >= box^2){
      imcut = matrix(0, box, box)
    }
    profoundSegimPlot(image = imcut,
                      segim = magcutout(multiband_segim_images[[i]]$segim,
                                        loc = c(xcen,ycen),
                                        box = box)$image, flip = T)
    legend(x = "topright",names(idx)[i], col="red", cex = 1.0)
    # text(x = box/2, y = box/4, names(idx)[i], col="red", cex = 1.3
  }
  
  dev.off()
  print("end")
  
}


k = 45
detect_image = Rfits_read_image(paste0("/Volumes/RAIDY/JWST/ProFound/Detects/",galaxies$VID[k], "/", galaxies$MODULE[k], "/",
                                         galaxies$VID[k], "_", galaxies$MODULE[k], "_profound_stack.fits"))
det_segim = readRDS(
    paste0("/Volumes/RAIDY/JWST/ProFound/Detects/", galaxies$VID[k], "/", galaxies$MODULE[k], "/", galaxies$VID[k], "_", galaxies$MODULE[k], "_profound.rds")
  )$segim
segim_wcs = Rfits_create_image(
  data = det_segim,
  keyvalues = detect_image$keyvalues
)
cutout = detect_image[galaxies$RAmax[k], galaxies$Decmax[k], box = 100, type="coord"]
segim_cutout = segim_wcs[galaxies$RAmax[k], galaxies$Decmax[k], box = 100, type="coord"]$imDat
plot(cutout,
       flip=T)
profoundSegimPlot(
    NULL,
    segim = segim_cutout,
    add = T
  )