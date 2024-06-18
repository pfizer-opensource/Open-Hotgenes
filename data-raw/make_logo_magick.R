
# install.packages(c("hexSticker", "magick"))

# library(magick)
library(hexSticker)

img <- file.path(getwd(),"data-raw", "HotgenesLogo.png")
hexSticker_out <- file.path(getwd(), "inst", "logo", "HotgenesLogo_hexSticker.png")

# ?sticker
p <- sticker(img,
  
#h_fill = "#A9A9A9",
# h_color = "#919191",
h_fill = "#FFFB00",
h_color = "black",
  
h_size = 1,
  p_y = 1, p_x = 1,
  asp = 1,
  s_x = 1, s_y = 1.15,
  s_width = 1.4, s_height = 1.4,
  filename = hexSticker_out, 
  package = ""
  )


print(p)

# ?usethis::use_logo()
usethis::use_logo(hexSticker_out, geometry = "240x278", retina = TRUE)

