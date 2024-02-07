# Efeitos individuais com modelo de coeficientes aleatórios para dados de painel
# Programa que faz uma ilustração empírica para a quarta parte do capítulo 7
# Amostragem de Gibbs para priori hierárquica normal independente

# Carregar os dados em y e x (todos empilhados por indivíduo)

# Observe também que a notação é ligeiramente diferente do livro
# theta é usado para coeficientes que variam entre os indivíduos

# Carregar dados em y e x (todos empilhados por indivíduo)
set.seed(123)  # Define a semente para garantir a reprodutibilidade dos resultados
setwd("C:/Users/02678854159/Documents/Empirico")
ch7data2 <- read.csv("dadosmatlab.csv")
#ch7data2 <- read_excel("dadosmatlab1.xlsx", 
#                           col_names = FALSE)

tn <- nrow(ch7data2)
t <- 5
n <- 100
y <- ch7data2[,1]
x <- ch7data2[,2:3]
kx <- ncol(x)

# --------
# Define os hiperparâmetros da priori
# --------

# Priori hierárquica para coeficientes variáveis
# A notação usada aqui tem theta0 sendo a média de theta na priori hierárquica
# theta0 tem priori normal com média mu_theta e variância V_theta
mu_theta <- rep(0, kx)
V_theta <- diag(1, kx)
V_thinv <- solve(V_theta)

# A notação aqui tem Sigma sendo a variância de theta na priori hierárquica
# Sigma-inversa tem priori que é Wishart
# grau de liberdade = rho, matriz de escala R --- implicando média = rho*R
rho <- 2
R <- 0.5 * diag(kx)
# para a precisão do erro, use priori Gamma com média h02 e v0 = graus de liberdade
v0 <- 1
h02 <- 25
s02 <- 1/h02

x <- as.matrix(x)
y <- as.matrix(y)

# Faça OLS e resultados relacionados (assumindo nenhuma heterogeneidade) para obter valores iniciais
bols <- solve(t(x) %*% x) %*% t(x) %*% y
s2 <- t(y - x %*% bols) %*% (y - x %*% bols) / (tn - kx)

# Escolha um valor inicial para h
hdraw <- 1/s2
# Calcule algumas quantidades fora do loop para uso posterior
xsquare <- t(x) %*% x
v1 <- v0 + tn
v0s02 <- v0 * s02
vrho <- rho + n
# Valor inicial para theta
thetdraw <- matrix(rep(c(bols[1,1], bols[2,1]), each = n), nrow = 2, byrow = TRUE)
thet0draw <- bols

# Se imlike for igual a 1, então calcule a verossimilhança marginal, caso contrário, não calcule
imlike <- 0

# Especifique o número de repetições
# número de repetições de burn-in
s0 <- 100
# número de repetições retidas
s1 <- 500
s <- s0 + s1

# Armazene todos os resultados nas seguintes matrizes
# inicialize-as aqui
h_ <- numeric() 
th0_ <- matrix(numeric(), nrow = kx, ncol = 0)
sig_ <- matrix(numeric(), nrow = kx^2, ncol = 0)
postth0 <- 0
# Apenas calcule a média e o desvio padrão para todos os coeficientes variáveis
# requer muito armazenamento para armazenar todos os resultados
thmean <- matrix(0, nrow = n, ncol = kx)
thsd <- matrix(0, nrow = n, ncol = kx)


#thetdraw <- as.matrix(thetdraw)
#thetdraw <- t(thetdraw)
#y <- as.matrix(y)
#x <- as.matrix(x)

#dtheta0 <- n * siginv %*% rowMeans(thetdraw) + V_thinv %*% mu_theta

# Agora comece o loop de Gibbs
set.seed(123)
for (irep in 1:s) {
  cat("irep: ", irep, "\n")
  
  # Amostragem de Sigma-inversa usando Wishart
  sigterm <- matrix(0, nrow = kx, ncol = kx)
  for (i in 1:n) {
    sigterm <- sigterm + (thetdraw[, i] - thet0draw) %*% t(thetdraw[, i] - thet0draw)
  }
  sigterm1 <- solve(sigterm + rho*R)
  siginv <- wish_rnd(sigterm1, vrho)
  sigdraw <- solve(siginv)
  
  # Amostragem de h condicional a outros parâmetros
  s12 <- 0
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
    thetuse <- thetdraw[, i]
    s12 <- s12 + t(yuse - xuse %*% thetuse) %*% (yuse - xuse %*% thetuse)
  }
  s12 <- (s12 + v0s02) / v1
  hdraw <- gamm_rnd(1,1,.5 * v1, .5 * v1 * s12)
  
  # Agora, desenhe theta0 (média na priori hierárquica) condicional a outros parâmetros
  capd0 <- solve(n * siginv + V_thinv)
  thetdrawt <- t(thetdraw)
  thetdrawtm <- colMeans(thetdrawt)
  thetdrawtm <- as.matrix(thetdrawtm)
  dtheta0 <- n * siginv %*% thetdrawtm + V_thinv %*% mu_theta
  th0mean <- capd0 %*% dtheta0
  thet0draw <- capd0 %*% dtheta0 + norm_rnd(capd0)
  
  # Agora, desenhe os thetas-i condicionais a outros parâmetros
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
    capdt <- solve(as.vector(hdraw) * (t(xuse) %*% xuse) + siginv)
    dt <- as.vector(hdraw) * (t(xuse) %*% yuse) + siginv %*% thet0draw
    thetdraw[, i] <- capdt %*% dt + norm_rnd(capdt)
  }
  
  if (irep > s0) {
    # após descartar o período de aquecimento, armazene todos os resultados
    h_ <- cbind(h_, hdraw)
    th0_ <- cbind(th0_, thet0draw)
    temp <- matrix(sigdraw, nrow = kx^2)
    sig_ <- cbind(sig_, temp)
    thmean <- thmean + t(thetdraw)
    thsd <- thsd + t(thetdraw^2)
    
    if (imlike == 1) {
      
      # Para o método Chib, isso calcula a parcela posterior relacionada a theta0
      logpost <- -0.5 * kx * log(2 * pi) - 0.5 * log(det(capd0)) - 
        0.5 * t(th0chib - th0mean) %*% solve(capd0) %*% (th0chib - th0mean)
      postth0 <- postth0 + logpost
    } 
  }
} 

thmean <- thmean / s1
thsd <- thsd / s1
thsd <- sqrt(thsd - thmean^2)


h_ <- t(h_)
th0_ <- t(th0_)
sig_ <- t(sig_)

alldraws <- cbind(h_, th0_, sig_)
# A função momentg é retirada da caixa de ferramentas de LeSage
# ela recebe todos os desenhos Gibbs e produz a média posterior,
# desvio padrão, nse e rne
# ela calcula o que Geweke chama de S(0) de várias maneiras
# consulte momentg.m para mais detalhes

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


#Se imlike==1 então calcule a probabilidade marginal, se não, então não marglike

imlike <- 1

#Use o método Chib (1995) para cálculo de probabilidade marginal
#isso requer ponto para avaliar tudo em --- try ols ou média posterior
#chibval7d.out é um arquivo onde salvei a média posterior da execução anterior


if (imlike == 1) {
  # Carregar os valores da média posterior de um arquivo
  chibval7d <- resumo_posterior$medias
  chibval7d <- as.vector(chibval7d)
  # Extrair os valores de hchib e th0chib
  hchib <- chibval7d[1]
  th0chib <- chibval7d[2:(kx + 1)]

  # Capital Sigma inverse
  sigchib <- diag(kx)
  sigchib[1:kx, 1] <- chibval7d[(kx + 2):(2 * kx + 1)]
  sigchib[1:kx, 2] <- chibval7d[((2 * kx) + 2):((3 * kx) + 1)]
  sigichib <- solve(sigchib)
  
  # Log da priori avaliada neste ponto
  logprior <- -.5 * v0 * log(2 * h02 / v0) - lgamma(.5 * v0) + .5 * (v0 - 2) * log(hchib) -
    .5 * v0 * hchib / h02 - .5 * kx * log(2 * pi) - .5 * kx * log(det(V_theta)) -
    .5 * t(th0chib - mu_theta) %*% V_thinv %*% (th0chib - mu_theta)
  
  sigtemp <- logwish_pdf(sigichib, R, rho)
  logprior <- logprior + sigtemp
  
  # Log da verossimilhança avaliada neste ponto
  loglike <- 0
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
    
    capri <- xuse %*% sigchib %*% t(xuse) + (1 / hchib) * diag(t)
    capriinv <- solve(capri)
    loglike <- loglike - .5 * t * log(2 * pi) - .5 * log(det(capri)) -
      .5 * t(yuse - xuse %*% th0chib) %*% capriinv %*% (yuse - xuse %*% th0chib)
  }
}
  
  if (imlike == 1) {
    postth0 <- postth0 / s1
  
    #a primeira execução auxiliar de Gibbs é avaliar o componente sigma-inverso
    # necessário se o método Chib para cálculo de verossimilhança marginal for usado

    postsig <- 0
    for (irep in 1:s) {
      # Inicialização da matriz sigterm
      sigterm <- matrix(0, nrow = kx, ncol = kx)
  
      # Loop para calcular sigterm
      for (i in 1:n) {
        sigterm <- sigterm + (thetdraw[,i] - th0chib) %*% t(thetdraw[,i] - th0chib)
      }
  
      # Cálculo de sigterm1
      sigterm1 <- solve(sigterm + rho * R)
  
      # Amostragem da inversa da Wishart
      siginv <- wish_rnd(sigterm1, vrho)
  
      # Cálculo da inversa da matriz amostrada
      sigdraw <- solve(siginv)
  

      #draw from h conditional on beta
      s12 <- 0
  
        # Loop para calcular s12
        for (i in 1:n) {
          # Seleção dos subconjuntos de dados x e y
          xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
          yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
          thetuse <- thetdraw[, i]
    
        # Cálculo da contribuição de cada subconjunto para s12
        s12 <- s12 + t(yuse - xuse %*% thetuse) %*% (yuse - xuse %*% thetuse)
        }
  
    s12 <- (s12 + v0s02) / v1
    hdraw <- gamm_rnd(1,1,.5 * v1, .5 * v1 * s12)
  
    # Loop para desenhar theta_i s
    for (i in 1:n) {
      # Seleção dos subconjuntos de dados x e y
      xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
      yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
      # Cálculo de capdt
      capdt <- solve(as.vector(hdraw) * (t(xuse) %*% xuse) + siginv)
      # Cálculo de dt
      dt <- as.vector(hdraw) * (t(xuse) %*% yuse) + siginv %*% th0chib
      # Desenho de thetadraw
      thetdraw[, i] <- capdt %*% dt + norm_rnd(capdt)
    }
  
    if (irep > s0) {
      post <- logwish_pdf(sigichib, sigterm1, vrho)
      postsig <- postsig + post
    }
  }

  postsig <- postsig / s1

  # O segundo algoritmo de amostragem de Gibbs auxiliar para avaliar o componente de precisão do erro
  # Isso é necessário para calcular a verossimilhança marginal usando o método de Chib

  posth <- 0
  for (irep in 1:s) {
    irep
    # Inicialização da variável s12
    s12 <- 0

      # Loop para calcular s12
      for (i in 1:n) {
      # Seleção dos subconjuntos de dados x e y
      xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
      yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
      thetuse <- thetdraw[, i]
      # Cálculo da contribuição de cada subconjunto para s12
      s12 <- s12 + t(yuse - xuse %*% thetuse) %*% (yuse - xuse %*% thetuse)
      }

    s12 <- (s12 + v0s02) / v1
    hdraw <- gamm_rnd(1,1,.5 * v1, .5 * v1 * s12)

  # Agora, gere as amostras de theta-i
  for (i in 1:n) {
    xuse <- x[((i - 1) * t + 1):(i * t), , drop = FALSE]
    yuse <- y[((i - 1) * t + 1):(i * t), , drop = FALSE]
  
    capdt <- solve(as.vector(hdraw) * (t(xuse) %*% xuse) + sigichib)
    dt <- as.vector(hdraw) * t(xuse) %*% yuse + sigichib %*% th0chib
  
    thetdraw[,i] <- capdt %*% dt + norm_rnd(capdt)
    }

  if (irep > s0) {
    logpost <- -.5 * v1 * log(2 / (v1 * s12)) - lgamma(.5 * v1) + .5 * (v1 - 2) * log(hchib) -
    .5 * v1 * hchib * s12
    posth <- posth + logpost
    }
  }
  posth <- posth / s1
  logmlike <- loglike + logprior - postth0 - postsig - posth
  }
  
  alldraws <- cbind(h_, th0_, sig_)
  # A função momentg é retirada da caixa de ferramentas de LeSage
  # ela recebe todos os desenhos Gibbs e produz a média posterior,
  # desvio padrão, nse e rne
  # ela calcula o que Geweke chama de S(0) de várias maneiras
  # consulte momentg.m para mais detalhes
  
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
  
  thetadraw
  
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
  
  # Iniciar uma nova janela de gráficos (se necessário)
  # par(mfrow = c(1, 1))  # Descomente se quiser dividir a tela em várias janelas de gráficos
  
  # Criar o histograma
  hist(thmean[, 1], breaks = 20, main = 'Figura 1: Histograma das Médias Posteriores de Alpha(i)s',
       xlab = 'Alpha(i)', col = 'lightblue', border = 'black')
  
  # Adicionar título e rótulos aos eixos
  title(main = 'Figura 1: Histograma das Médias Posteriores de Alpha(i)s', xlab = 'Alpha(i)')
  
  # Adicionar uma grade se desejar
  # grid()
  
  # Adicionar outras personalizações conforme necessário
  
  # Se preferir, você pode usar ggplot2 para mais flexibilidade e opções
  # library(ggplot2)
  # ggplot(data.frame(alpha = thmean[, 1]), aes(x = alpha)) +
  #   geom_histogram(binwidth = (max(thmean[, 1]) - min(thmean[, 1]))/20,
  #                  fill = 'lightblue', color = 'black') +
  #   ggtitle('Figura 1: Histograma das Médias Posteriores de Alpha(i)s') +
  #   xlab('Alpha(i)')
  
  
  