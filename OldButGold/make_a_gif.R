library(fields) # for tim.colors
library(caTools) # for write.gif
m = 400 # grid size
C = complex( real=rep(seq(-1.8,0.6, length.out=m), each=m ), imag=rep(seq(-1.2,1.2, length.out=m), m ) )
C = matrix(C,m,m)

Z = 0
X = array(0, c(m,m,20))
for (k in 1:20) {
  Z = Z^2+C
  X[,,k] = exp(-abs(Z))
}

AAA <- image(X[,,k], col=tim.colors(256))
is(AAA)
dim(I)
# show final image in R
write.gif(X, "Mandelbrot.gif", col=tim.colors(256), delay=100)

#######################################################################################
#######################################################################################

#library(installr)
#install.ImageMagick()
library(magick)

#crearing working directory
dir.create("examples")
setwd("examples") 

#Creating countdown .png files from 10 to "GO!"
png(file="example%02d.png", width=200, height=200)
for (i in c(10:1, "G0!")){
  plot.new()
  text(.5, .5, i, cex = 6)
}
dev.off()

# Converting .png files in one .gif image using ImageMagick
system("-delay 80 *.png example_1.gif")
# Remove .png files from working directory
file.remove(list.files(pattern=".png"))


# C:\Users\nicco\Documents\R\win-library\3.6\BH\include\boost\python\converter

#######################################################################################
#######################################################################################

# https://stackoverflow.com/questions/28142300/r-gif-animation-with-imagemagick

#library(devtools)
#dev_mode(on=T)
#install.packages('animation', repos = 'http://yihui.name/xran')
#library(animation)
#saveGIF({
#  for (i in 1:10) plot(runif(10), ylim = 0:1)
#})
#dev_mode(on=F)

library(animation)

saveGIF(
  expr = {
           for(iter in 1:nn){
             matplot(t(mu_j[[iter]]), type='l', main=paste0("Traceplot of mu(t) for cluster 1, iter=",iter) )
           }
    },
  movie.name = "cluster_1_mu.gif",
  interval = 0.1
  )

# TODO: fissare gli assi
# magari far vedere solo i clusters finali significativi?
# aggiungere iterazionis

is(X)
dim(X)
View(X)


length(mu_j)
is(mu_j)
is(mu_j[[1]])
dim(mu_j[[1]])

YO <- array(0, c())
