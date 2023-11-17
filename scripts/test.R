library(data.table)
library(foreach)
library(magicaxis)
library(cooltools)

sfs_fits = fread("/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/sfs_fits.csv")
sfs_z5 = fread("/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/sfs_z5.csv")
sfs_z7 = fread("/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/sfs_z7.csv")
sfs_z10 = fread("/Users/22252335/Documents/PROJ2_JWST_CSFH/data/save/sfs_z10.csv")

sfs_data = list(
  sfs_z5, sfs_z7, sfs_z10
)

keys = c("z5", "z7", "z10")

test_mstar = sfs_fits$mstar

foo = function(i){

  y_fit = sfs_fits[[paste0(keys[i], "_withAGN")]]
  y_fit_q16 = sfs_fits[[paste0(keys[i], "_withAGN_q16")]]
  y_fit_q84 = sfs_fits[[paste0(keys[i], "_withAGN_q84")]]

  lin_fit = lm(y_fit~test_mstar)
  lin16_fit = lm(y_fit_q16~test_mstar)
  lin84_fit = lm(y_fit_q84~test_mstar)

  x_data = sfs_data[[i]]$withAGN_mstar
  y_data = sfs_data[[i]]$withAGN_sfr

  x_err = sfs_data[[i]]$withAGN_mstar_err
  y_err = sfs_data[[i]]$withAGN_sfr_err

  y_test = lin_fit$coefficients[2]*x_data + lin_fit$coefficients[1]
  y_test16 = lin16_fit$coefficients[2]*x_data + lin16_fit$coefficients[1]
  y_test84 = lin84_fit$coefficients[2]*x_data + lin84_fit$coefficients[1]

  residuals = ((y_data - y_test)) / sqrt(x_err^2 + y_err^2)

  print(sum(residuals > 1))
  print(length(residuals))

  residuals_unnorm = y_data-y_test
  med_res = quantile(residuals_unnorm, 0.5)
  q5_res = quantile(residuals_unnorm, 0.05)
  q95_res = quantile(residuals_unnorm, 0.95)

  outlier=sum(
    residuals_unnorm < q5_res | residuals_unnorm > q95_res
  )
  print(outlier)
}
bar = function(i){

  y_fit = sfs_fits[[paste0(keys[i], "_withoutAGN")]]
  y_fit_q16 = sfs_fits[[paste0(keys[i], "_withoutAGN_q16")]]
  y_fit_q84 = sfs_fits[[paste0(keys[i], "_withoutAGN_q84")]]

  lin_fit = lm(y_fit~test_mstar)
  lin16_fit = lm(y_fit_q16~test_mstar)
  lin84_fit = lm(y_fit_q84~test_mstar)

  x_data = sfs_data[[i]]$withAGN_mstar
  y_data = sfs_data[[i]]$withAGN_sfr

  x_err = sfs_data[[i]]$withAGN_mstar_err
  y_err = sfs_data[[i]]$withAGN_sfr_err

  y_test = lin_fit$coefficients[2]*x_data + lin_fit$coefficients[1]
  y_test16 = lin16_fit$coefficients[2]*x_data + lin16_fit$coefficients[1]
  y_test84 = lin84_fit$coefficients[2]*x_data + lin84_fit$coefficients[1]

  residuals = ((y_data - y_test)) / sqrt(x_err^2 + y_err^2)

  print(sum(residuals > 1))
  print(length(residuals))

    residuals_unnorm = y_data-y_test
  med_res = quantile(residuals_unnorm, 0.5)
  q5_res = quantile(residuals_unnorm, 0.05)
  q95_res = quantile(residuals_unnorm, 0.95)

  outlier=sum(
    residuals_unnorm < q5_res | residuals_unnorm > q95_res
  )
  print(outlier)

}

foo(3)
bar(2)

foo = function(i){

  y_fit = sfs_fits[[paste0(keys[i], "_withAGN")]]
  y_fit_q16 = sfs_fits[[paste0(keys[i], "_withAGN_q16")]]
  y_fit_q84 = sfs_fits[[paste0(keys[i], "_withAGN_q84")]]

  lin_fit = lm(y_fit~test_mstar)
  lin16_fit = lm(y_fit_q16~test_mstar)
  lin84_fit = lm(y_fit_q84~test_mstar)

  x_data = sfs_data[[i]]$withAGN_mstar
  y_data = sfs_data[[i]]$withAGN_sfr

  x_err = sfs_data[[i]]$withAGN_mstar_err
  y_err = sfs_data[[i]]$withAGN_sfr_err

  y_test = lin_fit$coefficients[2]*x_data + lin_fit$coefficients[1]
  y_test16 = lin16_fit$coefficients[2]*x_data + lin16_fit$coefficients[1]
  y_test84 = lin84_fit$coefficients[2]*x_data + lin84_fit$coefficients[1]

  residuals = ((y_data - y_test))
  magplot(x_data, residuals)

  print(sum(residuals > 1))
  print(length(residuals))

}



magplot(
  sfs_z5$withoutAGN_mstar, sfs_z5$withoutAGN_sfr
)
magerr(
  sfs_z5$withoutAGN_mstar, sfs_z5$withoutAGN_sfr, xlo = sfs_z5$withoutAGN_mstar_err,
  ylo = sfs_z5$withoutAGN_sfr_err
)
lines(
  sfs_fits$mstar, sfs_fits$z5_withoutAGN, lwd=3
)
lines(
  sfs_fits$mstar, sfs_fits$z5_withoutAGN_q84, lwd=3
)
lines(
  sfs_fits$mstar, sfs_fits$z5_withoutAGN_q16, lwd=3
)

magplot(
  sfs_z7$withoutAGN_mstar, sfs_z7$withoutAGN_sfr
)
magerr(
  sfs_z7$withoutAGN_mstar, sfs_z7$withoutAGN_sfr, xlo = sfs_z7$withoutAGN_mstar_err,
  ylo = sfs_z7$withoutAGN_sfr_err
)
lines(
  sfs_fits$mstar, sfs_fits$z7_withoutAGN, lwd=3
)
lines(
  sfs_fits$mstar, sfs_fits$z7_withoutAGN_q84, lwd=3
)
lines(
  sfs_fits$mstar, sfs_fits$z7_withoutAGN_q16, lwd=3
)

magplot(
  sfs_z7$withoutAGN_mstar, sfs_z7$withoutAGN_sfr
)
magerr(
  sfs_z7$withoutAGN_mstar, sfs_z7$withoutAGN_sfr, xlo = sfs_z7$withoutAGN_mstar_err,
  ylo = sfs_z7$withoutAGN_sfr_err
)
lines(
  sfs_fits$mstar, sfs_fits$z7_withoutAGN, lwd=3
)
lines(
  sfs_fits$mstar, sfs_fits$z7_withoutAGN_q84, lwd=3
)
lines(
  sfs_fits$mstar, sfs_fits$z7_withoutAGN_q16, lwd=3
)