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
library(RColorBrewer)


## Do a comparison with the literature results
literature_galaxies = data.frame(
  list(
    RA  = c(
      hms2deg(14,20, 57.76),
      hms2deg(0, 14, 19.27),
      hms2deg(0, 14, 18.52),
      hms2deg(14,19, 14.19),
      hms2deg(14, 20, 42.77),
      hms2deg(14, 19, 33.52),
      hms2deg(14, 19, 17.63),
      hms2deg(14, 19, 20.69),
      hms2deg(14, 29, 19.54),
      hms2deg(14, 20, 34.87),
      215.0353914
    ),
    DEC = c(
      dms2deg(53, 2, 9.8),
      dms2deg(-30, 25, 27.8),
      dms2deg(-30, 25, 21.3),
      dms2deg(52, 52, 6.5),
      dms2deg(53, 3, 33.7),
      dms2deg(52, 49, 58.7),
      dms2deg(52, 49, 49),
      dms2deg(52, 52, 57.7),
      dms2deg(52, 58, 19.9),
      dms2deg(52, 58, 2.2),
      52.8906618
    ),
    AGNLum = c(
      5.8e44,
      3.4e44,
      2.8e44,
      3.3e45,
      1.4e45,
      1.3e45,
      8.3e44,
      8.3e44,
      2.1e44,
      5.4e43,
      10^(45.1)
    ),
    AGNLum_hi = c(
      9.7e44,
      14.8e44,
      3.0e44,
      0.3e45,
      1.5e45,
      0.1e45,
      4.3e44,
      34.3e44,
      12.7e44,
      5.6e43,
      0.2 * log(10) * 10^45.1
    ),
    AGNLum_lo = c(
      1.7e44,
      2.6e44,
      2.2e44,
      0.6e45,
      1.1e45,
      0.4e45,
      4.2e44,
      6.7e44,
      1.3e44,
      1.7e43,
      0.2 * log(10) * 10^45.1
    ),
    Mstar = c(
      8.63,
      8.82,
      9.1,
      9.11,
      9.92,
      9.01,
      9.35,
      9.36,
      9.61,
      8.94,
      9.5
    )
  )
)

withAGN = data.frame(
  fread(
    "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/catalogues/ProSpect_highz_withAGN.csv"
  )
)

withoutAGN = data.frame(
  fread(
    "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/catalogues/ProSpect_highz_withoutAGN.csv"
  )
)

match_literature = coordmatch(
  coordref = literature_galaxies[, c("RA", "DEC")],
  coordcompare = withAGN[, c("RA", "DEC")]
)

match_withAGN = withAGN[match_literature$bestmatch$compareID, ]
match_withoutAGN = withoutAGN[match_literature$bestmatch$compareID, ]
match_harikane = literature_galaxies[match_literature$bestmatch$refID, ]

magplot(
  log10(match_withAGN$AGNlum50),
  log10(match_harikane$AGNLum),
  xlim = c(35, 48),
  ylim = c(35, 48),
  pch = 16,
  cex = 1.2,
  col = brewer.pal(6, "Dark2")
)
magerr(
  log10(match_withAGN$AGNlum50),
  log10(match_harikane$AGNLum),
  xlo = match_withAGN$AGNlum16/(log(10)*match_withAGN$AGNlum50),
  xhi = match_withAGN$AGNlum84/(log(10)*match_withAGN$AGNlum50),
  ylo = match_harikane$AGNLum_lo/(log(10)*match_harikane$AGNLum),
  yhi = match_harikane$AGNLum_hi/(log(10) * match_harikane$AGNLum),
  col = brewer.pal(6, "Dark2")
  )
abline(0,1)
legend(
  x = "bottomright",
  paste(
    match_withAGN$COLID, "-", match_withAGN$VISITID, "-", match_withAGN$MODULE
  ),
  pch = 16,
  col = brewer.pal(6, "Dark2")

)

filelist = paste0(
  "/Volumes/RAIDY/JWST/ProFound/Detects/",
  match_withAGN$VISITID, "/", match_withAGN$MODULE, "/", match_withAGN$VISITID, "_", match_withAGN$MODULE,
  "_profound_stack.fits"
)

detect_list = Rfits_make_list(
  filelist = filelist
)

png("/Users/22252335/Documents/PROJ2_JWST_CSFH/plots/AGNcutouts.png",
    width = 10, height = 8, units = "in", res = 240)
par(mfrow = c(2,3), mar = rep(2,4), oma = rep(2,4))
for(i in 1:6){
  magimage(detect_list[[i]][match_withAGN$RA[i], match_withAGN$DEC[i], box = 50, type="coord"]$imDat,
       flip = T)
  legend(
    x = "topright",
    legend = paste0(match_withAGN$COLID[i], "_", match_withAGN$VISITID[i], "_", match_withAGN$MODULE[i])
  )
}
dev.off()





nakajima = data.frame(fread("/Users/22252335/Documents/PROJ2_JWST_CSFH/data/literature/nakajima.csv", stringsAsFactors = F, skip = 1))
match_literature = coordmatch(
  coordref = nakajima[, 2:3],
  coordcompare = withAGN[, c("RA", "DEC")]
)


nakajima_mass = nakajima[[8]][match_literature$bestmatch$refID]
match_mass_withoutAGN = log10(withoutAGN$Mstar50[match_literature$bestmatch$compareID])
match_mass_err_withoutAGN = ( withoutAGN$Mstar16[match_literature$bestmatch$compareID]
+ withoutAGN$Mstar84[match_literature$bestmatch$compareID]) / (log(10) * withoutAGN$Mstar50[match_literature$bestmatch$compareID])

match_mass_withAGN = log10(withAGN$Mstar50[match_literature$bestmatch$compareID])
match_mass_err_withAGN = ( withAGN$Mstar16[match_literature$bestmatch$compareID]
+ withAGN$Mstar84[match_literature$bestmatch$compareID]) / (log(10) * withAGN$Mstar50[match_literature$bestmatch$compareID])


idx_fix = grepl("<", nakajima_mass)
fix_upperlim = sapply(nakajima_mass[idx_fix], function(x)str_split(x,"<")[[1]][2])
nakajima_mass[idx_fix] = fix_upperlim
nakajima_mass = as.double(nakajima_mass)

pdf("/Users/22252335/Documents/PROJ2_JWST_CSFH/plots/mass_vs_mass_nakajima.pdf",
    height = 7, width = 10)
par(mfrow = c(1,2), mar = c(0,0,4,0), oma = c(4.0,4.0,1.0,1.0))
magplot(
  nakajima_mass,
  match_mass_withAGN,
  xlim = c(4.5,10.5),
  ylim = c(4.5,10.5),
  col="red",
  pch=16,
  main = "ProSpect Stellar+AGN"
)
magerr(
  nakajima_mass,
  match_mass_withAGN,
  xlo = nakajima[[9]][match_literature$bestmatch$refID],
  xhi = nakajima[[10]][match_literature$bestmatch$refID],
  ylo = match_mass_err_withAGN,
    col = "red"
)
abline(0,1)

magplot(
  nakajima_mass,
  match_mass_withoutAGN,
  xlim = c(4.5,10.5),
  ylim = c(4.5,10.5),
  col="cornflowerblue",
  pch=16,
  main = "ProSpect Stellar"
)
magerr(
  nakajima_mass,
  match_mass_withoutAGN,
  xlo = nakajima[[9]][match_literature$bestmatch$refID],
  xhi = nakajima[[10]][match_literature$bestmatch$refID],
  ylo = match_mass_err_withoutAGN,
      col = "cornflowerblue"

)
abline(0,1)

mtext("Nakajima JWST mass", side = 1, outer = T, line = 2.0)
mtext("ProSpect stellar mass", side = 2, outer = T, line = 2.0)

dev.off()

## mass vs. mass from EGS HST catalogue
ceers_egs_cat = data.frame(
  fread(
    "/Users/22252335/Downloads/hlsp_candels_hst_wfc3_egs_multi_v1_mass-cat.txt",
    skip = 39
  )
)

match_egs = coordmatch(
  coordref = ceers_egs_cat[, c(2,3)],
  coordcompare = withAGN[, c("RA", "DEC")],
  rad = 1.0
)

match_egs_cat = ceers_egs_cat[match_egs$bestmatch$refID, ]
match_galaxies_withAGN = withAGN[match_egs$bestmatch$compareID, ]
match_galaxies_withoutAGN = withoutAGN[match_egs$bestmatch$compareID, ]

redshift_match_idx = (abs(match_egs_cat$V12 - match_galaxies_withAGN$z)/match_galaxies_withAGN$z) <= 0.3

match_egs_cat = match_egs_cat[redshift_match_idx, ]
match_galaxies_withAGN = match_galaxies_withAGN[redshift_match_idx, ]
match_galaxies_withoutAGN = match_galaxies_withoutAGN[redshift_match_idx, ]

magplot(
  match_egs_cat$V12,
  match_galaxies_withAGN$z
)
abline(0,1)
mtext("HST EGS Cat stellar mass", side = 1, outer = T, line = 2.0)
mtext("ProSpect stellar mass", side = 2, outer = T, line = 2.0)

pdf(
  "/Users/22252335/Documents/PROJ2_JWST_CSFH/plots/mass_vs_mass_hst_egs.pdf",
  width = 10,
  height = 7
)
par(mfrow = c(1,2), mar = c(0,0,4,0), oma = c(4.0,4.0,1.0,1.0))
magplot(
  match_egs_cat$V18,
  log10(match_galaxies_withAGN$Mstar50),
  xlim = c(4.5,10.5),
  ylim = c(4.5,10.5),
  col="red",
  pch=16,
  main = "ProSpect Stellar+AGN"
)
magerr(
  match_egs_cat$V18,
  log10(match_galaxies_withAGN$Mstar50),
  xlo = match_egs_cat$V19,
  ylo = (match_galaxies_withAGN$Mstar16+match_galaxies_withAGN$Mstar84)/(log(10)*match_galaxies_withAGN$Mstar50),
  xlim = c(4.5,10.5),
  ylim = c(4.5,10.5),
  col="red",
)
abline(0,1, lty=2)

magplot(
  match_egs_cat$V18,
  log10(match_galaxies_withoutAGN$Mstar50),
  xlim = c(4.5,10.5),
  ylim = c(4.5,10.5),
  col="cornflowerblue",
  pch=16,
  main = "ProSpect Stellar"
)
magerr(
  match_egs_cat$V18,
  log10(match_galaxies_withoutAGN$Mstar50),
  xlo = match_egs_cat$V19,
  ylo = (match_galaxies_withoutAGN$Mstar16+match_galaxies_withoutAGN$Mstar84)/(log(10)*match_galaxies_withoutAGN$Mstar50),
  xlim = c(4.5,10.5),
  ylim = c(4.5,10.5),
  col="cornflowerblue"
)
abline(0,1, lty=2)
mtext("HST EGS Cat stellar mass", side = 1, outer = T, line = 2.0)
mtext("ProSpect stellar mass", side = 2, outer = T, line = 2.0)

dev.off()