source(here::here("src/global.R"))
library(rmarkdown)
render("master.Rmd", output_file="slides.html")
