
# install.packages(c("hexSticker", "magick"))

# library(magick)
library(hexSticker)

img <- file.path(getwd(),"data-raw", "HotgenesLogo.png")
hexSticker_out <- file.path(getwd(), "inst", "logo", "HotgenesLogo_hexSticker.png")

# ?sticker
p<- sticker(img,
  
#h_fill = "#A9A9A9",
# h_color = "#919191",
h_fill = "white",
h_color = "white",
  
h_size = 4,
  p_y = 1, p_x = 1,
  asp = 1,
  s_x = 1, s_y = 1.24,
  s_width = 0.96, s_height = 0.96,
  filename = hexSticker_out, 
  package = ""
  )


print(p)
#?usethis::use_logo()
usethis::use_logo(hexSticker_out, geometry = "240x278", retina = TRUE)

