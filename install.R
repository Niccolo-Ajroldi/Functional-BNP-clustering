
# required packages
needed_packages <- c(
  "MASS",
  "fda",
  "latex2exp",
  "roahd",
  "plyr",
  "coda",
  "devtools",
  "invgamma",
  "pbmcapply",
  "LaplacesDemon",
  "devtools"
)

# new packages
new_packages  <- needed_packages[!(needed_packages %in%installed.packages()[,"Package"])] 

# install required packages
if(length(new_packages))
{
  install.packages(new_packages)
}

# install mcclust using devtools
if(!c("mcclust.ext") %in% installed.packages())
{
  devtools::install_github("sarawade/mcclust.ext")
}

