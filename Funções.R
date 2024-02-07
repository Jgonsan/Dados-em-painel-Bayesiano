rm(list = ls())
set.seed(123)
gamm_rnd <- function(nrow, ncol, m, k) {
  if (missing(nrow) || missing(ncol) || missing(m) || missing(k)) {
    stop("Wrong number of arguments to gamm_rnd")
  }
  
  gb <- matrix(0, nrow = nrow, ncol = ncol)
  
  if (m <= 1) {
    # Use RGS algorithm by Best, p. 426
    c <- 1 / m
    t <- 0.07 + 0.75 * sqrt(1 - m)
    b <- 1 + exp(-t) * m / t
    
    for (i1 in 1:nrow) {
      for (i2 in 1:ncol) {
        accept <- FALSE
        while (!accept) {
          u <- runif(1)
          w <- runif(1)
          v <- b * u
          if (v <= 1) {
            x <- t * (v^c)
            accept <- (w <= ((2 - x) / (2 + x))) | (w <= exp(-x))
          } else {
            x <- -log(c * t * (b - v))
            y <- x / t
            accept <- ((w * (m + y - m * y)) <= 1) | (w <= (y^(m - 1)))
          }
        }
        gb[i1, i2] <- x
      }
    }
  } else {
    # Use Best's rejection algorithm XG, p. 410
    b <- m - 1
    c <- 3 * m - 0.75
    
    for (i1 in 1:nrow) {
      for (i2 in 1:ncol) {
        accept <- FALSE
        while (!accept) {
          u <- runif(1)
          v <- runif(1)
          w <- u * (1 - u)
          y <- sqrt(c / w) * (u - 0.5)
          x <- b + y
          if (x >= 0) {
            z <- 64 * (w^3) * v^2
            accept <- (z <= (1 - 2 * y^2 / x)) | (log(z) <= (2 * (b * log(x / b) - y)))
          }
        }
        gb[i1, i2] <- x
      }
    }
  }
  gb <- gb / k
  return(gb)
}

norm_rnd <- function(sig) {
  if (!is.matrix(sig) || nrow(sig) != ncol(sig)) {
    stop("Input 'sig' must be a square symmetric covariance matrix")
  }
  
  h <- chol(sig)
  x <- rnorm(nrow(sig))
  y <- t(h) %*% x
  
  return(y)
}


wish_rnd <- function(sigma, v) {
  # Função para gerar uma matriz aleatória Wishart
  # ------------------------------------------------
  # USO:     w = wish_rnd(sigma, v)
  # Onde: sigma = matriz simétrica definida positiva de entrada
  #        v = parâmetro de graus de liberdade
  # ------------------------------------------------
  # RETORNA:
  #        w = matriz distribuída aleatoriamente
  #            wishart_n(sigma)
  # ------------------------------------------------
  # REFERÊNCIAS: Gelman, Carlin, Stern, Rubin, Bayesian Data 
  #              Analysis, (1995,96) páginas 474, 480-481.
  
  n <- nrow(sigma)
  k <- ncol(sigma)
  
  if (n != k) {
    stop("wish_rnd: requer uma matriz quadrada")
  } else if (n < k) {
    warning("wish_rnd: n deve ser >= k+1 para uma distribuição finita")
  }
  
  t <- chol(sigma)
  if (any(is.na(t))) {
    stop("wish_rnd: a matriz deve ser definida positiva")
  }
  
  y <- t(t %*% matrix(rnorm(n * v), nrow = n))
  w <- t(y) %*% y
  return(w)
}

logwish_pdf <- function(Z, A, omega) {
  M <- nrow(Z)
  lintcon <- 0.5 * omega * M * log(2) + 0.25 * M * (M - 1) * log(pi)
  for (i in 1:M) {
    lintcon <- lintcon + lgamma(0.5 * (omega + 1 - i))
  }
  pdf <- -lintcon + 0.5 * (omega - M - 1) * log(det(Z)) - 
    0.5 * omega * log(det(A)) - 0.5 * sum(diag(solve(A) %*% Z))
  return(pdf)
}


momentg <- function(draws) {
  # Número de desenhos e variáveis
  ndraw <- nrow(draws)
  nvar <- ncol(draws)
  
  results <- list()
  results$ndraw <- ndraw
  results$meth <- 'momentg'
  results$nvar <- nvar
  
  NG <- 100
  
  if (ndraw < NG) {
    stop('momentg: precisa de um número maior de desenhos')
  }
  
  ntaper <- c(4, 8, 15)
  ns <- floor(ndraw / NG)
  nuse <- ns * NG
  
  for (jf in 1:nvar) {
    if (jf > length(results)) {
      results[[jf]] <- list()
    }
    cnt <- 0
    cn <- numeric(NG)
    cd <- numeric(NG)
    cdn <- numeric(NG)
    cdd <- numeric(NG)
    cnn <- numeric(NG)
    cvar <- numeric(NG)
    
    td <- 0
    tn <- 0
    tdn <- 0
    tdd <- 0
    tnn <- 0
    tvar <- 0
    
    for (ig in 1:NG) {
      gd <- 0
      gn <- 0
      gdn <- 0
      gdd <- 0
      gnn <- 0
      gvar <- 0
      
      for (is in 1:ns) {
        cnt <- cnt + 1
        g <- draws[cnt, jf]
        ad <- 1
        an <- ad * g
        gd <- gd + ad
        gn <- gn + an
        gdn <- gdn + ad * an
        gdd <- gdd + ad * ad
        gnn <- gnn + an * an
        gvar <- gvar + an * g
      }
      
      td <- td + gd
      tn <- tn + gn
      tdn <- tdn + gdn
      tdd <- tdd + gdd
      tnn <- tnn + gnn
      tvar <- tvar + gvar
      
      cn[ig] <- gn / ns
      cd[ig] <- gd / ns
      cdn[ig] <- gdn / ns
      cdd[ig] <- gdd / ns
      cnn[ig] <- gnn / ns
      cvar[ig] <- gvar / ns
    }
    
    eg <- tn / td
    varg <- tvar / td - eg^2
    sdg <- ifelse(varg > 0, sqrt(varg), -1)
    
    results[[jf]]$pmean <- eg
    results[[jf]]$pstd <- sdg
    
    varnum <- (tnn - 2 * eg * tdn + tdd * eg^2) / (td^2)
    sdnum <- ifelse(varnum > 0, sqrt(varnum), -1)
    
    results[[jf]]$nse <- sdnum
    results[[jf]]$rne <- varg / (nuse * varnum)
    
    # Calcular autocovariância das médias agrupadas
    barn <- tn / nuse
    bard <- td / nuse
    
    rnn <- numeric(NG)
    rdd <- numeric(NG)
    rnd <- numeric(NG)
    rdn <- numeric(NG)
    
    for (lag in 0:(NG - 1)) {
      ann <- 0
      add <- 0
      and <- 0
      adn <- 0
      
      for (ig in (lag + 1):NG) {
        ann <- ann + cn[ig] * cn[ig - lag]
        add <- add + cd[ig] * cd[ig - lag]
        and <- and + cn[ig] * cd[ig - lag]
        adn <- adn + cd[ig] * cd[ig - lag]
      }
      
      rnn[lag + 1] <- ann / NG
      rdd[lag + 1] <- add / NG
      rnd[lag + 1] <- and / NG
      rdn[lag + 1] <- adn / NG
    }
    
    for (mm in 1:3) {
      m <- ntaper[mm]
      am <- m
      snn <- rnn[1]
      sdd <- rdd[1]
      snd <- rnd[1]
      
      for (lag in 1:(m - 1)) {
        att <- 1 - lag / am
        snn <- snn + 2 * att * rnn[lag + 1]
        sdd <- sdd + 2 * att * rdd[lag + 1]
        snd <- snd + att * (rnd[lag + 1] + rnd[lag + 1])
      }
      
      varnum <- ns * nuse * (snn - 2 * eg * snd + sdd * eg^2) / (td^2)
      sdnum <- ifelse(varnum > 0, sqrt(varnum), -1)
      
      if (mm == 1) {
        results[[jf]]$nse1 <- sdnum
        results[[jf]]$rne1 <- varg / (nuse * varnum)
      } else if (mm == 2) {
        results[[jf]]$nse2 <- sdnum
        results[[jf]]$rne2 <- varg / (nuse * varnum)
      } else if (mm == 3) {
        results[[jf]]$nse3 <- sdnum
        results[[jf]]$rne3 <- varg / (nuse * varnum)
      }
    }
  }
  
  return(results)
}


