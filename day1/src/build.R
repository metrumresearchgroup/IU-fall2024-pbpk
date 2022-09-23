source(here::here("day1/src/global.R"))
library(rmarkdown)
render("master.Rmd", output_file="slides.html")
