
# Install packages ---------------------------------

#' Be sure to set this!

#.libPaths("")

pkg <- readLines("src/packages.txt")

install.packages(pkg)


# x <- list.files(".", pattern = "\\.R[md]*$", full.names=TRUE,recursive=TRUE)
# x <- x[x!="./src/packages.R"]
# code <- lapply(x,readLines) %>% unlist
# code <- code[grepl("library(", code,fixed = TRUE)]
# pkg <- gsub(".*library\\((\\w+)\\)", "\\1", code) %>% unique()
# pkg <- pkg[!(pkg %in% c("parallel"))]
# writeLines(pkg, con = "src/packages.txt")

