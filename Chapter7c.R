# Efeitos individuais com priori hierárquica para dados em painel
# Programa que faz uma ilustração empírica para a terceira parte do capítulo 7
# Amostragem de Gibbs para priori Normal-Gama independente

# Carregar os dados em y e z (todos empilhados por indivíduo)
# Este é na verdade um programa mais geral que permite que alguns
# coeficientes variem entre os indivíduos e outros sejam constantes
# x contém variáveis com coeficientes variáveis --- aqui definido como interceptação

# Observe também que a notação é ligeiramente diferente do livro
# theta é usado para coeficientes que variam entre os indivíduos
# gamma é usado para coeficientes constantes entre os indivíduos

# Carregar dados em y e x (todos empilhados por indivíduo)
set.seed(123)  # Define a semente para garantir a reprodutibilidade dos resultados
library(readr)
setwd("C:/Users/02678854159/Documents/Empirico")
ch7data2 <- read_csv("dadosmatlab.csv")
#ch7data2 <- read_excel("dadosmatlab1.xlsx", 
#                           col_names = FALSE)

tn <- nrow(ch7data2)
t <- 5
n <- 100
y <- ch7data2[, 1]
x <- ch7data2[, 2]
z <- ch7data2[, 3]
xz <- cbind(x, z)
kz <- ncol(z)
kx <- ncol(x)
k <- kz + kx

#--------
# Definir os hiperparâmetros da priori
#--------
# Para coeficientes constantes, Normal com média mu_gamma e variância V_gamma
mu_gamma <- rep(0, kz)
V_gamma <- diag(1, kz)
V_ginv <- solve(V_gamma)

# Priori hierárquica para coeficientes variáveis
# Normal com média mu_theta e variância V_theta
mu_theta <- rep(0, kx)
V_theta <- diag(1, kx)
V_thinv <- solve(V_theta)

# Para a precisão da heterogeneidade, use Wishart
# Grau de liberdade = rho, matriz de escala R --- implicando média = rho*R
# Observe que quando kx = 1, isso se reduz a uma priori Gamma
rho <- 2
R <- 0.5

# Para a precisão do erro, use priori Gamma com média h02 e v0 = grau de liberdade
v0 <- 1
h02 <- 25
s02 <- 1/h02

xz <- as.matrix(xz)
y <- as.numeric(y[[1]])

# Faça OLS e resultados relacionados (assumindo nenhuma heterogeneidade) para obter valores iniciais
bols <- solve(t(xz) %*% xz) %*% t(xz) %*% y
s2 <- t(y - xz %*% bols) %*% (y - xz %*% bols) / (tn - k)


# Escolha um valor inicial para h
hdraw <- 1/s2

# Calcule algumas quantidades fora do loop para uso posterior
xsquare <- t(xz) %*% xz
v1 <- v0 + tn
v0s02 <- v0 * s02
vrho <- rho + n

# Valor inicial para theta
thetdraw <- rep(bols[1], n)
thet0draw <- bols[1]

#Se imlike for igual a 1, então calcular a verossimilhança marginal; caso contrário, não calcular a verossimilhança marginal.

imlike <- 0

# Especificar o número de replicações
# Número de replicações para queima
s0 <- 100
# Número de replicações retidas
s1 <- 500
s <- s0 + s1

# Armazenar todos os resultados nas seguintes matrizes
# Inicializá-las aqui
g_ <- numeric(0)
h_ <- numeric(0)
th0_ <- numeric(0)
sig_ <- numeric(0)
postg <- 0
thmean <- matrix(0, nrow = n, ncol = kx)
thsd <- matrix(0, nrow = n, ncol = kx)

thetdraw <- as.matrix(thetdraw)
thetdraw <- t(thetdraw)
y <- as.matrix(y)
x <- as.matrix(x)
z <- as.matrix(z)
hdraw <- as.vector(hdraw)

#Now start Gibbs loop
for (irep in 1:s) {
  cat("irep: ", irep, "\n")
  
  # Desenhar de sigmainverse usando Wishart
  sigterm <- matrix(0, nrow = kx, ncol = kx)
  for (i in 1:n) {
    sigterm <- sigterm + (thetdraw[, i] - thet0draw) %*% t(thetdraw[, i] - thet0draw)
  }
  sigterm1 <- solve(sigterm + rho * R)
  siginv <- wish_rnd(sigterm1, vrho)
  sigdraw <- solve(siginv)
  
  # Desenhar de Gamma (coeficientes constantes) condicional a outros parâmetros
  capdginv <- matrix(0, nrow = kz, ncol = kz)
  dgamma <- rep(0, kz)
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    zuse <- z[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
    hdraw <- as.vector(hdraw)
    capri <- (xuse %*% sigdraw %*% t(xuse))+((1/hdraw) * diag(t))
    capriinv <- solve(capri)
    capdginv <- capdginv + t(zuse) %*% capriinv %*% zuse
    dgamma <- dgamma + t(zuse) %*% capriinv %*% (yuse - xuse %*% thet0draw)
  }
  capdg <- solve(capdginv + V_ginv)
  dgamma <- dgamma + V_ginv %*% mu_gamma
  gmean <- capdg %*% dgamma
  gdraw <- gmean + norm_rnd(capdg)
  
  # Desenhar de h condicional a outros parâmetros
  s12 <- 0
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    zuse <- z[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
    thetuse <- thetdraw[, i]
    s12 <- s12 + t(yuse - xuse %*% thetuse - zuse %*% gdraw) %*% (yuse - xuse %*% thetuse - zuse %*% gdraw)
  }
  s12 <- (s12 + v0s02) / v1
  hdraw <- gamm_rnd(1,1,.5 * v1, .5 * v1 * s12)
  
  # Agora desenhar theta0 (média na priori hierárquica) condicional a outros parâmetros
  capd0 <- solve(n * siginv + V_thinv)
  dtheta0 <- n * siginv %*% mean(thetdraw) + V_thinv %*% mu_theta
  thet0draw <- capd0 %*% dtheta0 + norm_rnd(capd0)
  
  # Agora desenhar theta-i s condicional a outros parâmetros
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    zuse <- z[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
    capdt <- solve(hdraw * t(xuse) %*% xuse + siginv)
    dt <- hdraw * t(xuse) %*% (yuse - zuse %*% gdraw) + siginv %*% thet0draw
    thetdraw[, i] <- capdt %*% dt + norm_rnd(capdt)
  }
  
  if (irep > s0) {
    # após descartar a queima, armazenar todos os desenhos
    g_ <- cbind(g_, gdraw)
    h_ <- cbind(h_, hdraw)
    th0_ <- cbind(th0_, thet0draw)
    temp <- matrix(sigdraw, nrow = kx^2)
    sig_ <- cbind(sig_, temp)
    thmean <- thmean + t(thetdraw)
    thsd <- thsd + t(thetdraw^2)
    
    if (imlike == 1) {
      # log posterior para gamma avaliado no ponto - usar para marg like
      # veja Chib (1995, JASA) pp. 1314-1316 para justificação
      logpostg <- -.5 * kz * log(2*pi) - .5 * kz * log(det(capdg)) -
        .5 * t(gchib - gmean) %*% solve(capdg) %*% (gchib - gmean)
      postg <- postg + logpostg
    }
  }
}

# Calculando a média
thmean <- thmean / s1

# Calculando o desvio padrão
thsd <- thsd / s1
thsd <- sqrt(thsd - thmean^2)


#Ajustando os paramentros 

g_ <- t(g_)
h_ <- t(h_)
th0_ <- t(th0_)
sig_ <- t(sig_)

colnames(g_) <- c("g_")
colnames(h_) <- c("h_")
colnames(th0_) <- c("th0_")
colnames(sig_) <- c("sig_")


# Concatena as amostras em uma matriz
alldraws <- cbind(g_, h_, th0_, sig_)

# Chama a função momentg para calcular estatísticas dos resultados
result <- momentg(alldraws)
means <- sapply(result, function(x) x$pmean)
stdevs <- sapply(result, function(x) x$pstd)
nse <- sapply(result, function(x) x$nse)
nse1 <- sapply(result, function(x) x$nse1)
nse2 <- sapply(result, function(x) x$nse2)
nse3 <- sapply(result, function(x) x$nse3)



# Resultados posteriores baseados na priori informativa
resultados_posteriores <- list(
  s0 = s0,
  s1 = s1
)

if (imlike == 1) {
  # Log da verossimilhança marginal
  resultados_posteriores$logmlike <- logmlike
}

# Médias posteriores, desvios padrão e NSE para os parâmetros
resumo_posterior <- list(
  medias = means,
  desvios_padrao = stdevs,
  nse = nse
)

# Combina tudo em uma única lista
resultados_analise <- list(
  resultados_posteriores = resultados_posteriores,
  resumo_posterior = resumo_posterior
)

resultados_analise



# Note que o próximo grande trecho de código realiza as execuções auxiliares
# necessárias para fazer o método de Chib para o cálculo da verossimilhança marginal
################################################################################

# Carregar dados em y e x (todos empilhados por indivíduo)
set.seed(123)  # Define a semente para garantir a reprodutibilidade dos resultados
setwd("C:/Users/02678854159/Documents/Empirico")
ch7data2 <- read_csv("dadosmatlab.csv")

tn <- nrow(ch7data2)
t <- 5
n <- 100
y <- ch7data2[, 1]
x <- ch7data2[, 2]
z <- ch7data2[, 3]
xz <- cbind(x, z)
kz <- ncol(z)
kx <- ncol(x)
k <- kz + kx

#--------
# Definir os hiperparâmetros da priori
#--------
# Para coeficientes constantes, Normal com média mu_gamma e variância V_gamma
mu_gamma <- rep(0, kz)
V_gamma <- diag(1, kz)
V_ginv <- solve(V_gamma)

# Priori hierárquica para coeficientes variáveis
# Normal com média mu_theta e variância V_theta
mu_theta <- rep(0, kx)
V_theta <- diag(1, kx)
V_thinv <- solve(V_theta)

# Para a precisão da heterogeneidade, use Wishart
# Grau de liberdade = rho, matriz de escala R --- implicando média = rho*R
# Observe que quando kx = 1, isso se reduz a uma priori Gamma
rho <- 2
R <- 0.5

# Para a precisão do erro, use priori Gamma com média h02 e v0 = grau de liberdade
v0 <- 1
h02 <- 25
s02 <- 1/h02

xz <- as.matrix(xz)
y <- as.numeric(y[[1]])

# Faça OLS e resultados relacionados (assumindo nenhuma heterogeneidade) para obter valores iniciais
bols <- solve(t(xz) %*% xz) %*% t(xz) %*% y
s2 <- t(y - xz %*% bols) %*% (y - xz %*% bols) / (tn - k)


# Escolha um valor inicial para h
hdraw <- 1/s2

# Calcule algumas quantidades fora do loop para uso posterior
xsquare <- t(xz) %*% xz
v1 <- v0 + tn
v0s02 <- v0 * s02
vrho <- rho + n

# Valor inicial para theta
thetdraw <- rep(bols[1], n)
thet0draw <- bols[1]

# Armazenar todos os resultados nas seguintes matrizes
# Inicializá-las aqui
g_ <- numeric(0)
h_ <- numeric(0)
th0_ <- numeric(0)
sig_ <- numeric(0)
postg <- 0
thmean <- matrix(0, nrow = n, ncol = kx)
thsd <- matrix(0, nrow = n, ncol = kx)

thetdraw <- as.matrix(thetdraw)
thetdraw <- t(thetdraw)
y <- as.matrix(y)
x <- as.matrix(x)
z <- as.matrix(z)
hdraw <- as.vector(hdraw)


#Se imlike for igual a 1, então calcular a verossimilhança marginal; caso contrário, não calcular a verossimilhança marginal.
################################################################################

imlike <- 1

# Usar o método de Chib (1995) para o cálculo da verossimilhança marginal
# Isso requer pontos para avaliação de tudo - tente ols ou média posterior
# chibval7c.out é um arquivo onde eu salvei a média posterior de uma execução anterior
if (imlike == 1) {
  chibval7c <- resumo_posterior$medias
  gchib <- chibval7c[1:kz]
  hchib <- chibval7c[kz + 1]
  th0chib <- chibval7c[(kz + 2):(kz + kx + 1)]
  
  
  # sigma inverso maiúsculo
  sigchib <- diag(kx)
  sigchib <- chibval7c[(k + 2):(k + kx + 1)]
  sigichib <- solve(sigchib)
  
  
  # log da priori avaliado neste ponto - use Poirier página 100 para a parte Gamma
  logprior <- -0.5 * v0 * log(2 * h02 / v0) - lgamma(0.5 * v0) +
    0.5 * (v0 - 2) * log(hchib) - 0.5 * v0 * hchib / h02 -
    0.5 * kz * log(2 * pi) - 0.5 * kz * log(det(V_gamma)) -
    0.5 * t(gchib - mu_gamma) %*% V_ginv %*% (gchib - mu_gamma) -
    0.5 * kx * log(2 * pi) - 0.5 * kx * log(det(V_theta)) -
    0.5 * t(th0chib - mu_theta) %*% V_thinv %*% (th0chib - mu_theta)
  
  R <- as.matrix(R)
  sigterm <- logwish_pdf(sigichib, R, rho)
  logprior <- logprior + sigterm
  
  # log da verossimilhança avaliado no ponto
  loglike <- 0
  
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    zuse <- z[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
    capri <- (xuse %*% sigchib %*% t(xuse))+((1/hchib) * diag(t))
    capriinv <- solve(capri)
    loglike <- loglike - .5 * t * log(2*pi) - .5 * log(det(capri)) -
      .5 * t(t(yuse - xuse %*% th0chib - zuse %*% gchib) %*% capriinv %*% (yuse - xuse %*% th0chib - zuse %*% gchib))
  }
}



#Now start Gibbs loop
for (irep in 1:s) {
  cat("irep: ", irep, "\n")
  
  # Desenhar de sigmainverse usando Wishart
  sigterm <- matrix(0, nrow = kx, ncol = kx)
  for (i in 1:n) {
    sigterm <- sigterm + (thetdraw[, i] - thet0draw) %*% t(thetdraw[, i] - thet0draw)
  }
  sigterm1 <- solve(sigterm + rho * R)
  siginv <- wish_rnd(sigterm1, vrho)
  sigdraw <- solve(siginv)
  
  # Desenhar de Gamma (coeficientes constantes) condicional a outros parâmetros
  capdginv <- matrix(0, nrow = kz, ncol = kz)
  dgamma <- rep(0, kz)
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    zuse <- z[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
    hdraw <- as.vector(hdraw)
    capri <- (xuse %*% sigdraw %*% t(xuse))+((1/hdraw) * diag(t))
    capriinv <- solve(capri)
    capdginv <- capdginv + t(zuse) %*% capriinv %*% zuse
    dgamma <- dgamma + t(zuse) %*% capriinv %*% (yuse - xuse %*% thet0draw)
  }
  capdg <- solve(capdginv + V_ginv)
  dgamma <- dgamma + V_ginv %*% mu_gamma
  gmean <- capdg %*% dgamma
  gdraw <- gmean + norm_rnd(capdg)
  
  # Desenhar de h condicional a outros parâmetros
  s12 <- 0
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    zuse <- z[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
    thetuse <- thetdraw[, i]
    s12 <- s12 + t(yuse - xuse %*% thetuse - zuse %*% gdraw) %*% (yuse - xuse %*% thetuse - zuse %*% gdraw)
  }
  s12 <- (s12 + v0s02) / v1
  hdraw <- gamm_rnd(1,1,.5 * v1, .5 * v1 * s12)
  
  # Agora desenhar theta0 (média na priori hierárquica) condicional a outros parâmetros
  capd0 <- solve(n * siginv + V_thinv)
  dtheta0 <- n * siginv %*% mean(thetdraw) + V_thinv %*% mu_theta
  thet0draw <- capd0 %*% dtheta0 + norm_rnd(capd0)
  
  # Agora desenhar theta-i s condicional a outros parâmetros
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    zuse <- z[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
    capdt <- solve(hdraw * t(xuse) %*% xuse + siginv)
    dt <- hdraw * t(xuse) %*% (yuse - zuse %*% gdraw) + siginv %*% thet0draw
    thetdraw[, i] <- capdt %*% dt + norm_rnd(capdt)
  }
  
  if (irep > s0) {
    # após descartar a queima, armazenar todos os desenhos
    g_ <- cbind(g_, gdraw)
    h_ <- cbind(h_, hdraw)
    th0_ <- cbind(th0_, thet0draw)
    temp <- matrix(sigdraw, nrow = kx^2)
    sig_ <- cbind(sig_, temp)
    thmean <- thmean + t(thetdraw)
    thsd <- thsd + t(thetdraw^2)
    
    if (imlike == 1) {
      # log posterior para gamma avaliado no ponto - usar para marg like
      # veja Chib (1995, JASA) pp. 1314-1316 para justificação
      logpostg <- -.5 * kz * log(2*pi) - .5 * kz * log(det(capdg)) -
        .5 * t(gchib - gmean) %*% solve(capdg) %*% (gchib - gmean)
      postg <- postg + logpostg
    }
  }
}

# Calculando a média
thmean <- thmean / s1

# Calculando o desvio padrão
thsd <- thsd / s1
thsd <- sqrt(thsd - thmean^2)

if (imlike == 1) {
  postg <- postg / s1

# Agora faça simulações auxiliares conforme descrito em Chib (1995)
# O primeiro auxiliar é avaliar o componente theta(0) 

postth0 <- 0

for (irep in 1:s) {
  cat("irep: ", irep, "\n")
  
  # Desenhar de sigmainverse usando Wishart
  sigterm <- matrix(0, nrow = kx, ncol = kx)
  for (i in 1:n) {
    sigterm <- sigterm + (thetdraw[, i] - thet0draw) %*% t(thetdraw[, i] - thet0draw)
  }
  sigterm1 <- solve(sigterm + rho * R)
  siginv <- wish_rnd(sigterm1, vrho)
  sigdraw <- solve(siginv)
  
  # Desenhar de h condicional em beta
  s12 <- 0
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    zuse <- z[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
    thetuse <- thetdraw[, i]
    s12 <- s12 + t(yuse - xuse %*% thetuse - zuse %*% gchib) %*% (yuse - xuse %*% thetuse - zuse %*% gchib)
  }
  s12 <- (s12 + v0s02) / v1
  hdraw <- gamm_rnd(1, 1, 0.5 * v1, 0.5 * v1 * s12)
  
  # Agora desenhe theta0 (média na priori hierárquica)
  capd0 <- solve(n * siginv + V_thinv)
  dtheta0 <- n * siginv %*% mean(thetdraw) + V_thinv %*% mu_theta
  th0mean <- capd0 %*% dtheta0
  thet0draw <- th0mean + norm_rnd(capd0)
  
  
  # Agora desenhe theta-i s
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    zuse <- z[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
    capdt <- solve(hdraw * t(xuse) %*% xuse + siginv)
    dt <- hdraw * t(xuse) %*% (yuse - zuse %*% gchib) + siginv %*% thet0draw
    thetdraw[, i] <- capdt %*% dt + norm_rnd(capdt)
  }
  
  if (irep > s0) {
    logpost <- -0.5 * kx * log(2 * pi) - 0.5 * kx * log(det(capd0)) -
      0.5 * t(th0chib - th0mean) %*% solve(capd0) %*% (th0chib - th0mean)
    postth0 <- postth0 + logpost
  }
}

postth0 <- postth0 / s1
# O segundo auxiliar Gibbs run é para avaliar o componente sigma-inverso
# necessário se o método de Chib para cálculo da verossimilhança marginal for usado
postsig <- 0

for (irep in 1:s) {
  cat("irep: ", irep, "\n")
  
  # Desenhar a partir da inversa da matriz sigma usando Wishart
  sigterm <- matrix(0, nrow = kx, ncol = kx)
  for (i in 1:n) {
    sigterm <- sigterm + (thetdraw[, i] - th0chib) %*% t(thetdraw[, i] - th0chib)
  }
  sigterm1 <- solve(sigterm + rho * R)
  siginv <- wish_rnd(sigterm1, vrho)
  sigdraw <- solve(siginv)
  
  # Desenhar de h condicional em beta
  s12 <- 0
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    zuse <- z[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
    thetuse <- thetdraw[, i]
    s12 <- s12 + t(yuse - xuse %*% thetuse - zuse %*% gchib) %*% (yuse - xuse %*% thetuse - zuse %*% gchib)
  }
  s12 <- (s12 + v0s02) / v1
  hdraw <- gamm_rnd(1, 1, 0.5 * v1, 0.5 * v1 * s12)
  
  # Agora desenhar os theta-i's
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    zuse <- z[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
    capdt <- solve(hdraw * t(xuse) %*% xuse + siginv)
    dt <- hdraw * t(xuse) %*% (yuse - zuse %*% gchib) + siginv %*% th0chib
    thetdraw[, i] <- capdt %*% dt + norm_rnd(capdt)
  }
  
  if (irep > s0) {
    post <- logwish_pdf(sigichib, sigterm1, vrho)
    postsig <- postsig + post
  }
}

postsig <- postsig / s1
# O terceiro amostrador auxiliar de Gibbs é para avaliar o componente de precisão do erro
# isso é necessário para avaliar a verossimilhança marginal usando o método Chib
posth <- 0

for (irep in 1:s) {
  irep
  
  # Desenhar de h condicional em beta
  s12 <- 0
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    zuse <- z[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
    thetuse <- thetdraw[, i]
    s12 <- s12 + t(yuse - xuse %*% thetuse - zuse * gchib) %*% (yuse - xuse %*% thetuse - zuse * gchib)
  }
  s12 <- (s12 + v0s02) / v1
  hdraw <- gamm_rnd(1, 1, .5 * v1, .5 * v1 * s12)
  
  # Agora desenhar theta-i s
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    zuse <- z[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
    capdt <- solve(hdraw * t(xuse) %*% xuse + sigichib)
    dt <- hdraw * t(xuse) %*% (yuse - zuse %*% gchib) + sigichib %*% th0chib
    thetdraw[, i] <- capdt %*% dt + norm_rnd(capdt)
  }
  if (irep > s0) {
    logpost <- -.5 * v1 * log(2 / (v1 * s12)) - lgamma(.5 * v1) + .5 * (v1 - 2) * log(hchib) - .5 * v1 * hchib * s12
    posth <- posth + logpost
  }
}
posth <- posth / s1
logmlike <- loglike + logprior - postg - postth0 - postsig - posth

}


#Ajustando os paramentros 

g_ <- t(g_)
h_ <- t(h_)
th0_ <- t(th0_)
sig_ <- t(sig_)

colnames(g_) <- c("g_")
colnames(h_) <- c("h_")
colnames(th0_) <- c("th0_")
colnames(sig_) <- c("sig_")


# Concatena as amostras em uma matriz
alldraws <- cbind(g_, h_, th0_, sig_)

# Chama a função momentg para calcular estatísticas dos resultados
result <- momentg(alldraws)
means <- sapply(result, function(x) x$pmean)
stdevs <- sapply(result, function(x) x$pstd)
nse <- sapply(result, function(x) x$nse)
nse1 <- sapply(result, function(x) x$nse1)
nse2 <- sapply(result, function(x) x$nse2)
nse3 <- sapply(result, function(x) x$nse3)



# Resultados posteriores baseados na priori informativa
resultados_posteriores <- list(
  s0 = s0,
  s1 = s1
)

if (imlike == 1) {
  # Log da verossimilhança marginal
  resultados_posteriores$logmlike <- logmlike
}

# Médias posteriores, desvios padrão e NSE para os parâmetros
resumo_posterior <- list(
  medias = means,
  desvios_padrao = stdevs,
  nse = nse
)

# Combina tudo em uma única lista
resultados_analise <- list(
  resultados_posteriores = resultados_posteriores,
  resumo_posterior = resumo_posterior
)

resultados_analise


# Cria um histograma dos valores de thmean
hist(thmean[,1], breaks = 20, main = 'Histogram of Posterior Means of Alpha(i)s, Hierarchical Prior', xlab = 'Alpha(i)')



