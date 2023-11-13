library(purrr)
library(stringr)

if (!dir.exists("data")) {
  dir.create("data")
}

options(timeout = 600)

# download/extract cellranger outputs
download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE237170&format=file",
  "data/GSE237170.tar",
  quiet = TRUE
)
untar("data/GSE237170.tar", exdir = "data/GSE237170")
file.remove("data/GSE237170.tar")
walk(
  list.files("data/GSE237170"),
  function(f) {
    file.rename(
      paste("data/GSE237170", f, sep = "/"),
      paste("data/GSE237170", str_replace(f, "GSM7596425_", ""), sep = "/")
    )
  }
)
