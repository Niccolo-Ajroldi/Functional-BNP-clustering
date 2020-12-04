setwd('C:/Users/Teresa Bortolotti/Documents/R/bayes_project/Functional-BNP-clustering')
load('smooth_60b_nopenalization.RData')

## ALIGNMENT

library(fda)
library(roahd)
library(fdakma)

t_ax <- 1:1600
y0 <- X_smoothed
y1 <- X_smoothed1

# Non siamo interessati a particolari comportamenti simili di crescita o decrescita delle
# nostre curve, quindi potrebbe avere senso porre la misura di similarità sulle funzioni
# originali e non sulle derivate. 
# A noi importa che le curve siano allineate in fase. Per questo motivo, 
# può avere senso pensare di fare un alignment con warping method "shift", che quindi
# consideri come warping functions solo delle traslazioni.

# Con metodo shift, per motivi di consistency, le distanze che posso provare sono L2 e 
# L2 centered (definizioni in pacchetto fdakma)
# Non capisco bene cosa in cosa differiscano in termini applicativi, comunque le provo
# entrambe e vedo che succede.

# Alla fine i tentativi che faccio sono:
# 1. d0.Pearson con metodo:
#     - affine -> trasformazioni affini dell'asse dei tempi h(t) = m*t + q
#     - shift --> h(t) = t + q
# 2. L2 sulle funzioni originali (metodo shift)
# 3. L2 centered sulle funcioni originali

# d0 Pearson, affine method

d0.pearson.affine <- kma(
  x=t_ax, y0=y0, n.clust = 1, 
  warping.method = 'affine', 
  similarity.method = 'd0.pearson',
  center.method = 'k-means'
)

#save(d0.pearson.affine, file = 'd0_pearson_affine.RData')
load('d0_pearson_affine.RData')

kma.show.results(d0.pearson.affine)
# questo non è male, però non mi dà il picco negativo iniziale, che è quello che vorrei
# vedere nel short latency evoked potential.
# Di buono è che non ha altri picchi oltre i primi due, che potrebbero confondere la
# mia analisi quando faccio clustering

# d0 Pearson, shift method

d0.Pearson.shift <- kma(
  x=t_ax, y0=y0, n.clust = 1, 
  warping.method = 'shift', 
  similarity.method = 'd0.pearson',
  center.method = 'k-means'
)

#save(d0.Pearson.shift, file = 'd0_pearson_shift.RData')
load('d0_pearson_shift.RData')

kma.show.results(d0.Pearson.shift)
# più brutto di quello di prima


# L2, shift method

d0.L2.shift <- kma(
  x=t_ax, y0=y0, n.clust = 1, 
  warping.method = 'shift', 
  similarity.method = 'd0.L2',
  center.method = 'k-means'
)

#save(d0.L2.shift, file = 'd0_L2_shift.RData')
load('d0_L2_shift.RData')

kma.show.results(d0.L2.shift)
# questo il migliore per ora, perché è il primo che mi fa
# vedere un picco negativo prima di quello positivo. Però non mi piace che abbia una terza
# deflessione dopo le prime due.

# L2 centered, shift method

d0.L2.centered.shift <- kma(
  x=t_ax, y0=y0, n.clust = 1, 
  warping.method = 'shift', 
  similarity.method = 'd0.L2.centered',
  center.method = 'k-means'
)

#save(d0.L2.centered.shift, file = 'd0_L2_centered_shift.RData')
load('d0_L2_centered_shift.RData')

kma.show.results(d0.L2.centered.shift)
# questo mi sembra identico a quello di prima

# COMMENTO:

# Per ora, quello che preferisco è un L2 (anche L2 centered ma non cambia nulla) perché
# mi trova il primo picco negativo.

# Questo alignment comunque non è ottimale, perché c'è una terza deflessione
# significativa dopo le prime due, e questa cosa non mi piace.

# Provo a vedere cosa succede se tiro in mezzo la derivata prima

# d1 Pearson
d1.pearson.affine <- kma(
  x=t_ax, y0=y0,y1=y1, n.clust = 1, 
  warping.method = 'affine', 
  similarity.method = 'd1.pearson',
  center.method = 'k-means'
)

#save(d1.pearson.affine, file = 'd1_pearson_affine.RData')
load('d1_pearson_affine.RData')

kma.show.results(d1.pearson.affine)

# L2, shift method

d1.L2.shift <- kma(
  x=t_ax, y0=y0,y1=y1, n.clust = 1, 
  warping.method = 'shift', 
  similarity.method = 'd1.L2',
  center.method = 'k-means'
)

#save(d1.L2.shift, file = 'd1_L2_shift.RData')
load('d1_L2_shift.RData')

kma.show.results(d1.L2.shift)

# No i risultati che ottengo sono peggio

# Applico d0.L2 ai dati di nico (fino a qui erano quelli di bieggie) per vedere se cambia qualcosa

# oggi non c'ho cazzi faccio domani
