
# byrow
aa <- c(1,2)
bb <- c(3,4)
cc <- cbind(aa,bb)
cc
matrix(cc, byrow=TRUE, nrow=2, ncol=2) # DIOCANE


# DIMOSTRAZIONE DELL'ERRORE
matplot(time.grid, t(X), type='l')
X.mat <- matrix(X, byrow=TRUE, nrow=nrow(X), ncol=ncol(X))
matplot(time.grid, t(X.mat), type='l')

#' praticamente riempivamo bene mu_coef
#' calcolavamo giusto la ricostruzione di mu valutata nei tempi
#' ma poi riempivamo la matrice al contrario


# ok: phi sembra giusta invece, perchè rinvgamma è un vettore
yeah <- rinvgamma(n=M*n_time, shape=c, rate=d)
phi.sbagliata <- matrix(yeah, byrow=TRUE, nrow=M, ncol=n_time)
View(yeah)
View(phi.sbagliata)
max(yeah-phi.sbagliata)


# sottrarre il massimo è completamente inutile, perchè p è spesso negativa, peggioriamo solo le cose
>       p_i
[1]  -38.799852  -83.295081 -167.080435 -125.043673
[5]   -9.037072 -183.628832 -234.159006 -173.816861
[9] -131.490358 -111.296659
> p_i - max(p_i)
[1]  -29.76278  -74.25801 -158.04336 -116.00660    0.00000
[6] -174.59176 -225.12193 -164.77979 -122.45329 -102.25959

# ora vorrei provare a sommare il minimo!
# comunque non penso migliori, perchè il problema è poi nel passare all'esponenziale mannaggia ladra
> p_i - min(p_i)
[1] 195.35915 150.86392  67.07857 109.11533 225.12193
[6]  50.53017   0.00000  60.34215 102.66865 122.86235

> exp(p_i)
[1]  1.410713e-17  6.689687e-37  2.740875e-73  4.945644e-55
[5]  1.189185e-04  1.782427e-80 2.023186e-102  3.253092e-76
[9]  7.842654e-58  4.618225e-49
> exp(p_i-max(p_i))
[1] 1.186285e-13 5.625437e-33 2.304834e-69 4.158850e-51
[5] 1.000000e+00 1.498864e-76 1.701321e-98 2.735564e-72
[9] 6.594980e-54 3.883520e-45
> exp(p_i-min(p_i))
[1] 6.972729e+84 3.306511e+65 1.354732e+29 2.444483e+47
[5] 5.877786e+97 8.810000e+21 1.000000e+00 1.607906e+26
[9] 3.876388e+44 2.282650e+53

# vorrei trovare una bella trasformazione monotona da applicare qui. Magari un log in una base potente??

