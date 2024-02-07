# Efeitos individuais com priori não hierárquica para dados em painel
# Programa que faz uma ilustração empírica para a segunda parte do capítulo 7
# Amostragem de Gibbs para priori Normal-Gamma independente
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
x <- ch7data2[, 3]

# Agora crie variáveis dummy específicas para cada indivíduo
xdummy <- matrix(0, nrow = tn, ncol = n)
for (i in 1:n) {
  xdummy[((i - 1) * t + 1):(i * t), i] <- rep(1, t)
}
x <- cbind(xdummy, x)
k <- ncol(x)

#--------
#Define os hiperparâmetros priori
#--------

# Hiperparâmetros para priori Normal-Gamma independente
v0 <- 1
b0 <- rep(0, k)
s02 <- 0.04
h02 <- 1 / s02
capv0 <- (1^2) * diag(k)
capv0inv <- solve(capv0)

x <- as.matrix(x)
y <- as.numeric(y[[1]])

# Realiza OLS e resultados relacionados para obter valores iniciais
bols <- solve(t(x) %*% x) %*% t(x) %*% y
s2 <- sum((y - x %*% bols)^2) / (tn - k)

# Escolhe um valor inicial para h
hdraw <- 1 / s2

# Calcula algumas quantidades fora do loop para uso posterior
xsquare <- t(x) %*% x
v1 <- v0 + tn
v0s02 <- v0 * s02

# Se imlike == 1, então calcula a verossimilhança marginal, se não, não calcula
imlike <- 0

# Usa o método de Chib (1995) para o cálculo da verossimilhança marginal
# Isso requer um ponto para avaliar tudo - tente os resultados do OLS ou use médias posteriores
if (imlike == 1) {
  # bchib <- bols
  # hchib <- 1 / s2
  chibval <- resumo_posterior$medias
  bchib <- as.vector(chibval[1:k])  # Convertendo para vetor coluna
  hchib <- chibval[k + 1]
  
  # Log da prior avaliado neste ponto
  logprior <- -0.5 * v0 * log(2 * h02 / v0) - lgamma(0.5 * v0) + 0.5 * (v0 - 2) * log(hchib) -
    0.5 * v0 * hchib / h02 - 0.5 * k * log(2 * pi) - 0.5 * k * log(det(capv0)) -
    0.5 * t(bchib - b0) %*% capv0inv %*% (bchib - b0)
  
  # Log da verossimilhança avaliado neste ponto
  loglike <- -0.5 * tn * log(2 * pi) + 0.5 * tn * log(hchib) -
    0.5 * hchib * t(y - x %*% bchib) %*% (y - x %*% bchib)
}

# Especificar o número de replicatas
# Número de replicatas de burnin
s0 <- 100
# Número de replicatas retidas
s1 <- 500
s <- s0 + s1

# Armazenar todos os sorteios nas seguintes matrizes
# Inicializá-las aqui
b_ <- matrix(numeric(0), nrow = s1, ncol = k)  # Inicializar matriz vazia para b
h_ <- numeric(0)  # Inicializar vetor vazio para h
logpost2 <- 0



# Agora inicie o loop de Gibbs
# beta condicional em h é Normal
# h condicional em beta é Normal

for (irep in 1:s) {
  cat("Iteration:", irep, "\n")
  
  # draw from beta conditional on h
  capv1inv <- capv0inv + hdraw * xsquare
  capv1 <- solve(capv1inv)
  b1 <- capv1 %*% (capv0inv %*% b0 + hdraw * xsquare %*% bols)
  bdraw <- b1 + norm_rnd(capv1)
  
  # draw from h conditional on beta
  s12 <- (sum((y - x %*% bdraw)^2) + v0s02) / v1
  hdraw <- gamm_rnd(1,1,.5*v1,.5*v1*s12)
  hdraw <- as.vector(hdraw)
  
  
  if (irep > s0) {
    # after discarding burnin, store all draws
    b_[irep - s0, ] <- bdraw
    h_ <- c(h_, hdraw)
    
    if (imlike == 1) {
      # log posterior for beta evaluated at point -- use for marg like
      # see Chib (1995, JASA) pp. 1315 for justification
      logpost <- -0.5 * k * log(2 * pi) - 0.5 * k * log(det(capv1)) -
        0.5 * t(bchib - b1) %*% capv1inv %*% (bchib - b1)
      logpost2 <- logpost2 + logpost
    }
  }
}

if (imlike == 1) {
  logpost2 <- logpost2 / s1
  
  # Precisamos de p(beta,h|y) avaliado como ponto
  # No loop, calculamos p(beta|y) agora precisamos de p(h|y,beta) para completar
  s12 <- (sum((y - x %*% bchib)^2) + v0s02) / v1
  logpost1 <- -0.5 * v1 * log(2 / (v1 * s12)) - lgamma(0.5 * v1) + 0.5 * (v1 - 2) * log(hchib) -
    0.5 * v1 * hchib * s12
  logmlike <- loglike + logprior - logpost2 - logpost1
}

# Unir todas as amostras em uma matriz
alldraws <- cbind(b_, h_)

medias <- colMeans(alldraws)
medias <- as.vector(medias)
print(medias)

desvios_padrao <- apply(alldraws, 2, sd)
desvios_padrao <- as.vector(desvios_padrao)
print(desvios_padrao)


summary(b_)
summary(h_)
sd(b_)
sd(h_)

result <- momentg(alldraws)
means <- sapply(result, function(x) x$pmean)
stdevs <- sapply(result, function(x) x$pstd)
nse <- sapply(result, function(x) x$nse)
nse1 <- sapply(result, function(x) x$nse1)
nse2 <- sapply(result, function(x) x$nse2)
nse3 <- sapply(result, function(x) x$nse3)

# Hiperparâmetros para a priori Normal-Gama independente
hiperparametros_priori <- list(
  b0 = b0,
  capv0 = capv0,
  v0 = v0,
  s02 = s02
)

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
  hiperparametros_priori = hiperparametros_priori,
  resultados_posteriores = resultados_posteriores,
  resumo_posterior = resumo_posterior
)

resultados_analise


# Extrai as médias dos parâmetros theta
thmean <- means[1:n]

# Cria um gráfico de histograma das médias posteriores de Alpha(i)s
hist(thmean, breaks = 20, main = "Figura 7.3: Histograma das Médias Posteriores de Alpha(i)s, Prior Não-hierárquica",
     xlab = "Alpha(i)", ylab = "Frequência")


################################################################################
################################################################################
set.seed(123)  # Define a semente para garantir a reprodutibilidade dos resultados
imlike <- 1
################################################################################
################################################################################
# Se imlike==1, então calcula a verossimilhança marginal, se não, não calcula
# Usa o método de Chib (1995) para o cálculo da verossimilhança marginal
# Isso requer pontos para avaliar tudo, tente resultados OLS ou use médias posteriores
if (imlike == 1) {
  # bchib <- bols
  # hchib <- 1 / s2
  chibval <- resumo_posterior$medias
  bchib <- as.vector(chibval[1:k])  # Convertendo para vetor coluna
  hchib <- chibval[k + 1]
  
  # Log da prior avaliado neste ponto
  logprior <- -0.5 * v0 * log(2 * h02 / v0) - lgamma(0.5 * v0) + 0.5 * (v0 - 2) * log(hchib) -
    0.5 * v0 * hchib / h02 - 0.5 * k * log(2 * pi) - 0.5 * k * log(det(capv0)) -
    0.5 * t(bchib - b0) %*% capv0inv %*% (bchib - b0)
  
  # Log da verossimilhança avaliado neste ponto
  loglike <- -0.5 * tn * log(2 * pi) + 0.5 * tn * log(hchib) -
    0.5 * hchib * t(y - x %*% bchib) %*% (y - x %*% bchib)
}

# Armazenar todos os desenhos nas seguintes matrizes

# Inicializar matrizes e variáveis
b_ <- matrix(numeric(0), nrow = s1, ncol = k)  # Inicializar matriz vazia para b
h_ <- numeric(0)  # Inicializar vetor vazio para h
logpost2 <- 0

# Especificar o número de repetições
# Número de repetições de burnin
s0 <- 100
# Número de repetições retidas
s1 <- 500
s <- s0 + s1

# Iniciar o loop Gibbs
for (irep in 1:s) {
  cat("Iteration:", irep, "\n")
  
  # draw from beta conditional on h
  capv1inv <- capv0inv + hdraw * xsquare
  capv1 <- solve(capv1inv)
  b1 <- capv1 %*% (capv0inv %*% b0 + hdraw * xsquare %*% bols)
  bdraw <- b1 + norm_rnd(capv1)
  
  # draw from h conditional on beta
  s12 <- (sum((y - x %*% bdraw)^2) + v0s02) / v1
  hdraw <- gamm_rnd(1,1,.5*v1,.5*v1*s12)
  hdraw <- as.vector(hdraw)
  
  
  if (irep > s0) {
    # after discarding burnin, store all draws
    b_[irep - s0, ] <- bdraw
    h_ <- c(h_, hdraw)
    
    if (imlike == 1) {
      # log posterior for beta evaluated at point -- use for marg like
      # see Chib (1995, JASA) pp. 1315 for justification
      logpost <- -0.5 * k * log(2 * pi) - 0.5 * k * log(det(capv1)) -
        0.5 * t(bchib - b1) %*% capv1inv %*% (bchib - b1)
      logpost2 <- logpost2 + logpost
    }
  }
}

if (imlike == 1) {
  logpost2 <- logpost2 / s1
  
  # Precisamos de p(beta,h|y) avaliado como ponto
  # No loop, calculamos p(beta|y) agora precisamos de p(h|y,beta) para completar
    s12 <- (sum((y - x %*% bchib)^2) + v0s02) / v1
  logpost1 <- -0.5 * v1 * log(2 / (v1 * s12)) - lgamma(0.5 * v1) + 0.5 * (v1 - 2) * log(hchib) -
    0.5 * v1 * hchib * s12
  logmlike <- loglike + logprior - logpost2 - logpost1
}


# Unir todas as amostras em uma matriz
alldraws <- cbind(b_, h_)

summary(b_)
summary(h_)
sd(b_)
sd(h_)


result <- momentg(alldraws)
means <- sapply(result, function(x) x$pmean)
stdevs <- sapply(result, function(x) x$pstd)
nse <- sapply(result, function(x) x$nse)
nse1 <- sapply(result, function(x) x$nse1)
nse2 <- sapply(result, function(x) x$nse2)
nse3 <- sapply(result, function(x) x$nse3)


# Hiperparâmetros para a priori Normal-Gama independente
hiperparametros_priori <- list(
  b0 = b0,
  capv0 = capv0,
  v0 = v0,
  s02 = s02
)

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
  hiperparametros_priori = hiperparametros_priori,
  resultados_posteriores = resultados_posteriores,
  resumo_posterior = resumo_posterior
)

resultados_analise

# Extrai as médias dos parâmetros theta
thmean <- means[1:n]

thmean

# Cria um gráfico de histograma das médias posteriores de Alpha(i)s
hist(thmean, breaks = 20, main = "Figura 7.3: Histograma das Médias Posteriores de Alpha(i)s, Prior Não-hierárquica",
     xlab = "Alpha(i)", ylab = "Frequência")
