library(ProFound)
library(Rfits)
library(Rwcs)
library(foreach)
library(doParallel)

old_stub = "/Volumes/BRAIDY/JWST/1345/CAL/NIRCAM/"
new_stub = "/Volumes/Expansion/JWST/6969/CAL/NIRCAM/"

file_names = list.files(
  new_stub,
  pattern = glob2rx("*1345001001*")
)

pdf("/Volumes/Expansion/JWST/6969/diff_im.pdf", width = 10, height = 10)

hist_pix = foreach(i = 1:length(file_names))%do%{
  
  old = Rfits_read_image(
    filename = paste0(old_stub, "/", file_names[i]), ext = 2
  )$imDat
  
  new = Rfits_read_image(
    filename = paste0(new_stub, "/", file_names[i]), ext = 2
  )$imDat
  
  par(mfrow = c(1,2), oma = rep(1,4), mar = rep(1,4))
  magimage(old-new, qdiff = T)
  text(1000, 100, "old-new")
  text(1000, 2000, file_names[i])
  hh=maghist(old - new)
  legend(x = "topright", 
         legend=paste(names(hh$summary), 
                      round(hh$summary,3)))
  
}
dev.off()
