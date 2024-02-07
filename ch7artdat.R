set.seed(123)  # Define a semente para garantir a reprodutibilidade dos resultados

n <- 100
t <- 5
tn <- t * n

# Define o tipo de conjunto de dados a ser gerado
dataset <- 2

if (dataset == 1 | dataset == 3) {
  theta <- c(0, 2)
  capsig <- 0.25 * diag(2)
  sigma <- 0.2
}

if (dataset == 2) {
  pidraw <- c(0.75, 0.25)
  theta1 <- c(1, 2)
  theta2 <- c(-1, 2)
  sigma <- 0.2
}

# Gera artificialmente os dados na variável explicativa
x <- cbind(rep(1, tn), matrix(runif(tn), ncol = 1))

# Gera os dados para cada indivíduo
y <- numeric(tn)

for (i in 1:n) {
  if (dataset == 1) {
    thetai <- c(rnorm(1, theta[1], sqrt(capsig[1, 1])), theta[2])
  }
  
  if (dataset == 3) {
    thetai <- rnorm(2, theta, sqrt(diag(capsig)))
  }
  
  if (dataset == 2) {
    if (runif(1) < pidraw[1]) {
      thetai <- theta1
    } else {
      thetai <- theta2
    }
  }
  
    y[(1 + t * (i - 1)):(i * t)] <- x[(1 + t * (i - 1)):(i * t), ] %*% thetai + sigma * rnorm(t)
}

# Empilha os dados em uma matriz
data <- cbind(y, x)

#Definir diretorio
#setwd("C:/Users/jhona/OneDrive/Documentos/Mestrado_UFSC/Econometria Bayesiana/Artigo/Empirico")

# Salva os dados em um arquivo CSV
write.csv(data, file = "ch7data2.csv", row.names = FALSE)


