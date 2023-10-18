img <- magick::image_read(paste0(this.path::this.dir(), "/flocker.png"))

plot(img)
st <- hexSticker::sticker(img, package = "flocker", p_color = "gray95",
              s_width = 1, s_height = 1, s_x = 1, p_y = 1.42, 
              h_color = "dodgerblue3", h_fill = "dodgerblue2",
              filename=paste0(this.path::this.dir(), "flocker_sticker.png"))
plot(st)
