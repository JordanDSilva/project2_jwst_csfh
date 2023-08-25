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
library(dftools)
library(scales)

detect_names = c(
  "1345001001_NRCA",
  "1345001001_NRCB",
  "1345002001_NRCA",
  "1345002001_NRCB",
  "1345003001_NRCA",
  "1345003001_NRCB",
  "1345004001_NRCA",
  "1345004001_NRCB",
  "1345061001_NRCA",
  "1345061001_NRCB",
  "1345064001_NRCA",
  "1345064001_NRCB",
  "1345066001_NRCA",
  "1345066001_NRCB",
  "1345068001_NRCA",
  "1345068001_NRCB",
  "1345069001_NRCA",
  "1345069001_NRCB",
  "1345071001_NRCA",
  "1345071001_NRCB",
  "1324001001_NRCA"
)

## Prepare data
load_catalogue = data.frame(
  cbind(
    bind_rows(
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345001001/NRCA/1345001001_NRCA_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345001001/NRCB/1345001001_NRCB_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345002001/NRCA/1345002001_NRCA_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345002001/NRCB/1345002001_NRCB_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345003001/NRCA/1345003001_NRCA_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345003001/NRCB/1345003001_NRCB_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345004001/NRCA/1345004001_NRCA_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345004001/NRCB/1345004001_NRCB_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345061001/NRCA/1345061001_NRCA_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345061001/NRCB/1345061001_NRCB_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345064001/NRCA/1345064001_NRCA_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345064001/NRCB/1345064001_NRCB_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345066001/NRCA/1345066001_NRCA_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345066001/NRCB/1345066001_NRCB_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345068001/NRCA/1345068001_NRCA_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345068001/NRCB/1345068001_NRCB_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345069001/NRCA/1345069001_NRCA_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345069001/NRCB/1345069001_NRCB_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345071001/NRCA/1345071001_NRCA_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1345071001/NRCB/1345071001_NRCB_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1324001001/NRCA/1324001001_NRCA_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1324001001/NRCB/1324001001_NRCB_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1324006001/NRCA/1324006001_NRCA_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1324006001/NRCB/1324006001_NRCB_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1324007001/NRCA/1324007001_NRCA_segstats.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Detects/1324007001/NRCB/1324007001_NRCB_segstats.csv")
    ),
    bind_rows(
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345001001/NRCA/1345001001_NRCA_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345001001/NRCB/1345001001_NRCB_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345002001/NRCA/1345002001_NRCA_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345002001/NRCB/1345002001_NRCB_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345003001/NRCA/1345003001_NRCA_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345003001/NRCB/1345003001_NRCB_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345004001/NRCA/1345004001_NRCA_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345004001/NRCB/1345004001_NRCB_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345061001/NRCA/1345061001_NRCA_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345061001/NRCB/1345061001_NRCB_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345064001/NRCA/1345064001_NRCA_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345064001/NRCB/1345064001_NRCB_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345066001/NRCA/1345066001_NRCA_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345066001/NRCB/1345066001_NRCB_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345068001/NRCA/1345068001_NRCA_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345068001/NRCB/1345068001_NRCB_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345069001/NRCA/1345069001_NRCA_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345069001/NRCB/1345069001_NRCB_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345071001/NRCA/1345071001_NRCA_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1345071001/NRCB/1345071001_NRCB_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1324001001/NRCA/1324001001_NRCA_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1324001001/NRCB/1324001001_NRCB_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1324006001/NRCA/1324006001_NRCA_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1324006001/NRCB/1324006001_NRCB_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1324007001/NRCA/1324007001_NRCA_photometry.csv"),
      fread("/Volumes/RAIDY/JWST/ProFound/Measurements/1324007001/NRCB/1324007001_NRCB_photometry.csv")
    )
  )
)

fincat = load_catalogue
wave_temp = sapply(c(cenwave$cenwave, EAZY_filters$cenwave), list)
names(wave_temp) = c(cenwave$filter, EAZY_filters$info)

grep("F.+fluxt$" , names(fincat), value=T)

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

arrabol_haro = data.frame(fread("/Users/22252335/Documents/PROJ2_JWST_CSFH/data/literature/ArrabalHaro+23b_Table6.csv"))
nakajima = data.frame(fread("/Users/22252335/Documents/PROJ2_JWST_CSFH/data/literature/nakajima.csv", stringsAsFactors = T, skip = 1))
harikane23 = data.frame(list(
  RA =  c(
    hms2deg(14, 20, 2.81),
    hms2deg(14, 20, 52.50),
    hms2deg(14, 19, 37.59),
    3.499
  ),
  Dec = c(
    dms2deg(52, 59, 12.9),
    dms2deg(53, 4, 11.5),
    dms2deg(52, 56, 56.43),
    -30.40072
  ),
  z = c(
    8.876, #ceers
    8.612,
    11.04,
    9.76 #glass
  )
))

glass_z12 = data.frame(list(
  RA = 3.499,
  Dec = -30.32475,
  z = 12.117
))
highz = data.frame(
  RA = c(arrabol_haro$ra, nakajima$deg, glass_z12$RA, harikane23$RA),
  DEC = c(arrabol_haro$dec, nakajima$deg.1,glass_z12$Dec, harikane23$Dec),
  z = c(arrabol_haro$z_spec, nakajima$V4, glass_z12$z, harikane23$z)
)
match = coordmatch(coordref = highz[,c("RA","DEC")],
                   coordcompare = fincat[,c("RAmax", "Decmax")])

galaxies = fincat[match$bestmatch$compareID, ]
galaxies$z = highz$z[match$bestmatch$refID]

i=10
magplot(
  lambda, unlist(galaxies[i, fluxt_names])
)
magerr(
  lambda, unlist(galaxies[i, fluxt_names]),
  ylo = unlist(galaxies[i, fluxt_err_names])
)


#Apply galactic extinction correction and 10% error floor
cat_for_dust_map = galaxies[, c("RAmax", "Decmax")]
names(cat_for_dust_map) = c("ra", "dec")
fwrite(cat_for_dust_map, "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/catalogues/cat_for_dustmap.csv")

extinction_table_names = data.frame(read.table("/Users/22252335/Documents/PROJ2_JWST_CSFH/data/catalogues/extinction.tbl", skip = 13, nrows =1, sep="|"))
extinction_table = data.frame(read.table("/Users/22252335/Documents/PROJ2_JWST_CSFH/data/catalogues/extinction.tbl", skip = 16))
names(extinction_table) =  extinction_table_names[3:17]

klambda = function(lambda){
  c1 = -0.175
  c2 = 0.807
  c3 = 2.991 
  c4 = 0.319 
  c5 = 6.097
  x0 = 4.592 
  gamma = 0.922
  O1 = 2.055
  O2 = 1.322
  O3 = 0.000
  RV = 3.001
  KIR = 1.057
  
  if (lambda < 2700){
    x = 1/(lambda*1e-4)
    drude = x^2 / ((x^2-x0^2)^2 + (x*gamma)^2)
    if (x<=c5){
      return(c1+c2*x+c3*drude)
    }else{
      return(c1 + c2*x + c3*drude + c4*(x-c5)^2)
    }
  }else{
    x = 1/(lambda*1e-4)
    lambda_UV = c(2600, 2699.99)*1e-4
    lambda_optical = c(3300, 4000, 5530)*1e-4
    lambda_IR = 1/ c(1e-10, 0.25, 0.50, 0.75, 1)
    xtest = c(1/lambda_UV, 1/(lambda_optical), 1e-323, 0.25, 0.50, 0.75, 1)
    ytest = c(3.77672, 3.435908, O1, O2, O3, KIR*(lambda_IR^-1.84)-RV)
    foo = cubicspline(sort(xtest), ytest[order(xtest)], x)
    return(foo)
  }
  
}

for(k in 1:dim(galaxies)[1]){
  Rv = 3.1
  E_bv = extinction_table$`E_B_V_SandF`[k]
  extinction_corr = (foreach(x = lambda, .combine = "c")%do%{klambda(x)} + Rv)*E_bv
  
  ext_corrt_flux = galaxies[k, fluxt_names] *( 1e-6 / 10^(-0.4*extinction_corr))
  ext_corrt_flux_err = sqrt( (galaxies[k, fluxt_err_names] *( 1e-6 / 10^(-0.4*extinction_corr)))^2 + (0.1 * galaxies[k, fluxt_names] *( 1e-6 / 10^(-0.4*extinction_corr)))^2 )
  
  #only for EAZY
  ext_corrc_flux = galaxies[k, fluxc_names] *( 1e-6 / 10^(-0.4*extinction_corr))
  ext_corrc_flux_err = sqrt( (galaxies[k, fluxc_err_names] *( 1e-6 / 10^(-0.4*extinction_corr)))^2)
  ext_corrc_flux[is.na(ext_corrc_flux)] = -999
  ext_corrc_flux_err[is.na(ext_corrc_flux_err)] = 0
  
  flux_names_corrt = paste0(fluxt_names, "_corr")
  flux_err_names_corrt = paste0(fluxt_err_names, "_corr")
  
  flux_names_corrc = paste0(fluxc_names, "_corr")
  flux_err_names_corrc = paste0(fluxc_err_names, "_corr")
  
  galaxies[k, flux_names_corrt] = ext_corrt_flux
  galaxies[k, flux_err_names_corrt] = ext_corrt_flux_err
  
  galaxies[k, flux_names_corrc] = ext_corrc_flux
  galaxies[k, flux_err_names_corrc] = ext_corrc_flux_err
}

flux_names_corrt = paste0(fluxt_names, "_corr")
flux_err_names_corrt = paste0(fluxt_err_names, "_corr")

flux_names_corrc = paste0(fluxc_names, "_corr")
flux_err_names_corrc = paste0(fluxc_err_names, "_corr")

fwrite(galaxies, file = "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/catalogues/high_z_profound.csv")
galaxies = data.frame( fread("/Users/22252335/Documents/PROJ2_JWST_CSFH/data/catalogues/high_z_profound.csv") )

## prepare EAZY photoz
filter_eazy_loc = c(233,
                    236, 
                    238, 
                    239, 
                    363,
                    202,
                    364,
                    203,
                    204,
                    365,
                    205,
                    366,
                    375,
                    376,
                    383,
                    377
)

eazy_translate = data.frame(cbind(
  c(flux_names_corrc, flux_err_names_corrc),
  c(paste0("F", filter_eazy_loc), paste0("E", filter_eazy_loc))
))

fwrite(eazy_translate, "/Users/22252335/Documents/eazy/eazy-photoz/translates/csfh_cagnh_highz_study.translate", sep = " ")
eazy_cat = galaxies[, c("RAmax", "Decmax", flux_names_corrc, flux_err_names_corrc)]
eazy_cat$ID = paste0(galaxies$visitid, ifelse(galaxies$module=="NRCA", "000", "001"))
fwrite(eazy_cat, "/Users/22252335/Documents/eazy/eazy-photoz/cats/csfh_cagnh_highz_study.csv")
zout = Rfits_read_table("/Users/22252335/Documents/eazy/eazy-photoz/csfh_cagnh_highz_study.zout.fits")

galaxies$zphot = zout$z_ml

magplot(galaxies$z, galaxies$zphot)
fwrite(galaxies, file = "~/Documents/CSFH_CAGNH_JWST/Data/highz_photometry_cat.csv")
galaxies = data.frame( fread("~/Documents/CSFH_CAGNH_JWST/Data//highz_photometry_cat.csv") )
