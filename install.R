
# automatic installation of required packages

needed_packages <- c(
  "MASS",
  "fda",
  "fdakma",
  "latex2exp",
  "roahd",
  "plyr",
  "coda",
  "devtools",
  "mcclust.ext",
  "invgamma",
  "pbmcapply",
  "LaplacesDemon"
)

new_packages  <- needed_packages[!(needed_packages %in%installed.packages()[,"Package"])] 

if(length(new_packages))
{
  install.packages(new_packages)
}

