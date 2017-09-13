Impute_man <- function(M_j0, D, Nw) {
  d <- dim(M_j0)[1]
  j0 <- log(dim(M_j0)[3], 2)
  tM_j1 <- array(dim = c(d, d, 2^j0))
  # P0 <- M_j0[, , 2^(j0 - 1)]
  # M_j0.log <- sapply(1:2^j0, function(i) Logm(P0, M_j0[, , i]), simplify = "array")
  if (D == 0) {
    for (i in 1:2^(j0 - 1)) {
      P0 <- M_j0[, , 2 * i - 1]
      tM_j1[, , 2 * i - 1] <- Expm(P0, (1/4) * Logm(P0, M_j0[, , 2 * i]))
      tM_j1[, , 2 * i] <- Expm(P0, (5/4) * Logm(P0, M_j0[, , 2 * i]))
    }
  } else {
    for (i in 1:2^j0) {
      if (i == 1) {
        tM_j1[, , i] <- Expm(M_j0[, , i], (1/4) * Logm(M_j0[, , i], M_j0[, , i + 1]))
      } else if (i == 2^j0) {
        tM_j1[, , i] <- Expm(M_j0[, , i], (-1/4) * Logm(M_j0[, , i], M_j0[, , i - 1]))
      } else {
        nbrs0 <- abs((0:2^j0) - (0:2^j0)[i]) <= D
        Di <- min(sum(nbrs0[1:(i - 1)]), sum(nbrs0[(i + 1):2^j0]))
        # nbrs <- abs((0:2^j0) - (0:2^j0)[i]) <= Di
        tM_j1[, , i] <- KarchMean(M = M_j0[, , (i - Di):(i + Di)], w = Nw[[Di + 1]])
      }
    }
  }
  return(tM_j1)
}

Impute2D <- function(M_j0, D){

  j0 <- log(dim(M_j0)[3], 2)
  d <- dim(M_j0)[1]

  if(identical(j0, 0)){

    tM_j1 <- array(M_j0[, , 1, 1], dim = c(d, d, 2, 2))

  } else {

    tM_j1 <- array(dim = c(d, d, 2^(j0 + 1), 2^(j0 + 1)))
    grid_k <- expand.grid(1:2^j0, 1:2^j0)

    for(i in 1:nrow(grid_k)){

      D_k <- c(min(grid_k[i, 1] - 1, D[1]), min(grid_k[i, 2] - 1, D[2]))
      nbrs <- list(x = max(grid_k[i, 1] - D[1], 1):min(grid_k[i, 1] + D[1], 2^j0),
                   y = max(grid_k[i, 2] - D[2], 1):min(grid_k[i, 2] + D[2], 2^j0))
      N_k <- 2 * D_k + 1

      if((length(nbrs$x) == 2 * D[1] + 1) & (length(nbrs$y) == 2 * D[2] + 1) & all(D <= 4)){

        ## Away from the boundary:
        ## Predict with available weights
        pred2D <- function(l1, l2){
          M <- array(c(M_j0[, , nbrs$x, nbrs$y]), dim = c(d, d, length(nbrs$x) * length(nbrs$y)))
          w <- c(pdSpecEst:::W_2D[[5 * D_k[2] + (D_k[1] + 1)]][l1, l2, , ])
          Mean_new <- M_j0[, , grid_k[i, 1], grid_k[i, 2]]
          Mean <- diag(d)
          m <- 0
          while((pdDist(Mean_new, Mean) > 1E-6) & (m < 10)){
            Mean <- Mean_new
            Mean_new <- Expm(Mean, apply(sapply(1:dim(M)[3], function(k) w[k] *
                                                  Logm(Mean, M[,,k]), simplify = "array"), c(1, 2), sum))
            m <- m+1
          }
          return(Mean_new)
        }

        tM_j1[, , (2 * grid_k[i, 1] - 1):(2 * grid_k[i, 1]), (2 * grid_k[i, 2] - 1):(2 * grid_k[i, 2])] <-
          array(c(mapply(pred2D, rep(c(1,2), times = 2), rep(c(1,2), each = 2), SIMPLIFY = "array")), dim = c(d,d,2,2))

      } else {

        ## At the boundary:
        ## Predict with Neville's algorithm
        grid_l <- expand.grid(nbrs$x, nbrs$y)
        M.bar <- array(c(mapply(function(i1,i2) KarchMean(array(c(M_j0[, , grid_l[1, 1]:i1, grid_l[1, 2]:i2]),
                                dim = c(d, d, length(grid_l[1, 1]:i1) * length(grid_l[1, 2]:i2)))), grid_l$Var1,
                                grid_l$Var2, SIMPLIFY="array")), dim = c(d, d, length(nbrs$x), length(nbrs$y)))
        M.bar_pred <- pdNeville2D(P = M.bar, ti = nbrs, grid = list(x = grid_k[i, 1] - 1 + c(0, 0.5, 1),
                                                                    y = grid_k[i, 2] - 1 + c(0, 0.5, 1)))

        ## Predict M_{j+1,2k_1+1,2k_2+1}
        C0 <- KarchMean(array(c(M.bar_pred[, , 3, 2], M.bar_pred[, , 2, 2]), dim=c(d, d, 2)), w = c(N_k[1] + 1, -N_k[1]))
        M.bar_min0 <- KarchMean(array(c(M.bar_pred[, , 2, 3], C0), dim=c(d, d, 2)),
                                w = c((N_k[1] * (N_k[2] + 1)) / ((N_k[1] + 1) * (N_k[2] + 1) - 1), N_k[2] / ((N_k[1] + 1) * (N_k[2] + 1) - 1)))
        tM_j1[, , 2 * grid_k[i, 1], 2 * grid_k[i, 2]] <- KarchMean(array(c(M.bar[, , which(nbrs$x == grid_k[i, 1]),
                                which(nbrs$y == grid_k[i, 2])], M.bar_min0), dim = c(d, d, 2)),
                                w = c((N_k[1] + 1) * N_k[2] + N_k[1] + 1, -((N_k[1] + 1) * N_k[2] + N_k[1])))

        ## Predict M_{j+1,2k_1,2k_2+1}
        C1 <- KarchMean(array(c(M.bar_pred[, , 1, 3], M.bar_pred[, , 1, 2]), dim = c(d, d, 2)), w = c(N_k[2] + 1, -N_k[2]))
        M.bar_min1 <- KarchMean(array(c(M.bar_pred[, , 3, 2], C1), dim=c(d, d, 2)),
                                w = c((N_k[1] + 1) * N_k[2] / ((N_k[1] + 1) * (N_k[2] + 1) - 2), (N_k[1] - 1) / ((N_k[1] + 1) * (N_k[2] + 1) - 2)))
        M.block1 <- KarchMean(array(c(M.bar[, , which(nbrs$x == grid_k[i, 1]), which(nbrs$y == grid_k[i, 2])], M.bar_min1),
                                dim = c(d, d, 2)), w = c(((N_k[1] + 1) * (N_k[2] + 1) - 2) / 2 + 1, -((N_k[1] + 1) * (N_k[2] + 1) - 2) / 2))
        tM_j1[, , 2 * grid_k[i, 1] - 1, 2 * grid_k[i, 2]] <- KarchMean(array(c(M.block1,
                                tM_j1[, , 2 * grid_k[i, 1], 2 * grid_k[i, 2]]), dim = c(d, d, 2)), w = c(2, -1))

        ## Predict M_{j+1,2k_1+1,2k_2}
        C2 <- KarchMean(array(c(M.bar_pred[, , 3, 1], M.bar_pred[, , 2, 1]), dim = c(d, d, 2)), w = c(N_k[1] + 1, -N_k[1]))
        M.bar_min2 <- KarchMean(array(c(M.bar_pred[, , 2, 3], C2), dim = c(d, d, 2)),
                                w = c((N_k[2] + 1) * N_k[1] / ((N_k[2] + 1) * (N_k[1] + 1) - 2), (N_k[2] - 1) / ((N_k[2] + 1) * (N_k[1] + 1) - 2)))
        M.block2 <- KarchMean(array(c(M.bar[, , which(nbrs$x == grid_k[i, 1]), which(nbrs$y == grid_k[i, 2])], M.bar_min2),
                                dim = c(d, d, 2)), w = c(((N_k[1] + 1) * (N_k[2] + 1) - 2) / 2 + 1, -((N_k[1] + 1) * (N_k[2] + 1) - 2) / 2))
        tM_j1[, , 2 * grid_k[i, 1], 2 * grid_k[i, 2] - 1] <- KarchMean(array(c(M.block2,
                                tM_j1[, , 2 * grid_k[i, 1], 2 * grid_k[i, 2]]), dim = c(d, d, 2)), w = c(2, -1))

        ## Predict M_{j+1,2k_1,2k_2}
        tM_j1[, , 2 * grid_k[i, 1] - 1, 2 * grid_k[i, 2] - 1] <- Expm(M_j0[, , grid_k[i, 1], grid_k[i, 2]],
                                -(Logm(M_j0[, , grid_k[i, 1], grid_k[i, 2]], tM_j1[, , 2 * grid_k[i, 1], 2 * grid_k[i, 2] - 1]) +
                                  Logm(M_j0[, , grid_k[i, 1], grid_k[i, 2]], tM_j1[, , 2 * grid_k[i, 1] - 1, 2 * grid_k[i, 2]]) +
                                  Logm(M_j0[, , grid_k[i, 1], grid_k[i, 2]], tM_j1[, , 2 * grid_k[i, 1], 2 * grid_k[i, 2]])))
      }
    }
  }

  return(tM_j1)
}

Impute2D_multicore <- function(M_j0, D, cores = 1){

  j0 <- log(dim(M_j0)[3], 2)
  d <- dim(M_j0)[1]

  if(identical(j0, 0)){

    tM_j1 <- array(M_j0[, , 1, 1], dim = c(d, d, 2, 2))

  } else {

    tM_j1 <- array(dim = c(d, d, 2^(j0 + 1), 2^(j0 + 1)))
    grid_k <- expand.grid(1:2^j0, 1:2^j0)

    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)

    tM1_par <- foreach::foreach(i=1:nrow(grid_k), .packages = "pdSpecEst") %dopar% {

      D_k <- c(min(grid_k[i, 1] - 1, D[1]), min(grid_k[i, 2] - 1, D[2]))
      nbrs <- list(x = max(grid_k[i, 1] - D[1], 1):min(grid_k[i, 1] + D[1], 2^j0),
                   y = max(grid_k[i, 2] - D[2], 1):min(grid_k[i, 2] + D[2], 2^j0))
      N_k <- 2 * D_k + 1

      if((length(nbrs$x) == 2 * D[1] + 1) & (length(nbrs$y) == 2 * D[2] + 1) & all(D <= 4)){

        ## Away from the boundary:
        ## Predict with available weights
        pred2D <- function(l1, l2){
          M <- array(c(M_j0[, , nbrs$x, nbrs$y]), dim = c(d, d, length(nbrs$x) * length(nbrs$y)))
          w <- c(W_2D[[5 * D_k[2] + (D_k[1] + 1)]][l1, l2, , ])
          Mean_new <- M_j0[, , grid_k[i, 1], grid_k[i, 2]]
          Mean <- diag(d)
          m <- 0
          while((pdDist(Mean_new, Mean) > 1E-6) & (m < 10)){
            Mean <- Mean_new
            Mean_new <- Expm(Mean, apply(sapply(1:dim(M)[3], function(k) w[k] *
                                                  Logm(Mean, M[,,k]), simplify = "array"), c(1, 2), sum))
            m <- m+1
          }
          return(Mean_new)
        }

        tM1 <- array(c(mapply(pred2D, rep(c(1,2), times = 2), rep(c(1,2), each = 2), SIMPLIFY = "array")), dim = c(d,d,2,2))

      } else {

        ## At the boundary:
        ## Predict with Neville's algorithm
        tM1 <- array(dim = c(d, d, 2, 2))
        grid_l <- expand.grid(nbrs$x, nbrs$y)
        M.bar <- array(c(mapply(function(i1,i2) KarchMean(array(c(M_j0[, , grid_l[1, 1]:i1, grid_l[1, 2]:i2]),
                                                                dim = c(d, d, length(grid_l[1, 1]:i1) * length(grid_l[1, 2]:i2)))), grid_l$Var1,
                                grid_l$Var2, SIMPLIFY="array")), dim = c(d, d, length(nbrs$x), length(nbrs$y)))
        M.bar_pred <- pdNeville2D(P = M.bar, ti = nbrs, grid = list(x = grid_k[i, 1] - 1 + c(0, 0.5, 1),
                                                                    y = grid_k[i, 2] - 1 + c(0, 0.5, 1)))

        ## Predict M_{j+1,2k_1+1,2k_2+1}
        C0 <- KarchMean(array(c(M.bar_pred[, , 3, 2], M.bar_pred[, , 2, 2]), dim=c(d, d, 2)), w = c(N_k[1] + 1, -N_k[1]))
        M.bar_min0 <- KarchMean(array(c(M.bar_pred[, , 2, 3], C0), dim=c(d, d, 2)),
                                w = c((N_k[1] * (N_k[2] + 1)) / ((N_k[1] + 1) * (N_k[2] + 1) - 1), N_k[2] / ((N_k[1] + 1) * (N_k[2] + 1) - 1)))
        tM1[, , 2, 2] <- KarchMean(array(c(M.bar[, , which(nbrs$x == grid_k[i, 1]),
                                which(nbrs$y == grid_k[i, 2])], M.bar_min0), dim = c(d, d, 2)),
                                w = c((N_k[1] + 1) * N_k[2] + N_k[1] + 1, -((N_k[1] + 1) * N_k[2] + N_k[1])))

        ## Predict M_{j+1,2k_1,2k_2+1}
        C1 <- KarchMean(array(c(M.bar_pred[, , 1, 3], M.bar_pred[, , 1, 2]), dim = c(d, d, 2)), w = c(N_k[2] + 1, -N_k[2]))
        M.bar_min1 <- KarchMean(array(c(M.bar_pred[, , 3, 2], C1), dim=c(d, d, 2)),
                                w = c((N_k[1] + 1) * N_k[2] / ((N_k[1] + 1) * (N_k[2] + 1) - 2), (N_k[1] - 1) / ((N_k[1] + 1) * (N_k[2] + 1) - 2)))
        M.block1 <- KarchMean(array(c(M.bar[, , which(nbrs$x == grid_k[i, 1]), which(nbrs$y == grid_k[i, 2])], M.bar_min1),
                                    dim = c(d, d, 2)), w = c(((N_k[1] + 1) * (N_k[2] + 1) - 2) / 2 + 1, -((N_k[1] + 1) * (N_k[2] + 1) - 2) / 2))
        tM1[, , 1, 2] <- KarchMean(array(c(M.block1, tM1[, , 2, 2]), dim = c(d, d, 2)), w = c(2, -1))

        ## Predict M_{j+1,2k_1+1,2k_2}
        C2 <- KarchMean(array(c(M.bar_pred[, , 3, 1], M.bar_pred[, , 2, 1]), dim = c(d, d, 2)), w = c(N_k[1] + 1, -N_k[1]))
        M.bar_min2 <- KarchMean(array(c(M.bar_pred[, , 2, 3], C2), dim = c(d, d, 2)),
                                w = c((N_k[2] + 1) * N_k[1] / ((N_k[2] + 1) * (N_k[1] + 1) - 2), (N_k[2] - 1) / ((N_k[2] + 1) * (N_k[1] + 1) - 2)))
        M.block2 <- KarchMean(array(c(M.bar[, , which(nbrs$x == grid_k[i, 1]), which(nbrs$y == grid_k[i, 2])], M.bar_min2),
                                    dim = c(d, d, 2)), w = c(((N_k[1] + 1) * (N_k[2] + 1) - 2) / 2 + 1, -((N_k[1] + 1) * (N_k[2] + 1) - 2) / 2))
        tM1[, , 2, 1] <- KarchMean(array(c(M.block2, tM1[, , 2, 2]), dim = c(d, d, 2)), w = c(2, -1))

        ## Predict M_{j+1,2k_1,2k_2}
        tM1[, , 1, 1] <- Expm(M_j0[, , grid_k[i, 1], grid_k[i, 2]],
                                -(Logm(M_j0[, , grid_k[i, 1], grid_k[i, 2]], tM1[, , 2, 1]) +
                                  Logm(M_j0[, , grid_k[i, 1], grid_k[i, 2]], tM1[, , 1, 2]) +
                                  Logm(M_j0[, , grid_k[i, 1], grid_k[i, 2]], tM1[, , 2, 2])))
      }
      return(tM1)
    }

    parallel::stopCluster(cl)

    for(l in 1:nrow(grid_k)){
      tM_j1[, , (2 * grid_k[l, 1] - 1):(2 * grid_k[l, 1]), (2 * grid_k[l, 2] - 1):(2 * grid_k[l, 2])] <- tM1_par[[l]]
    }

  }

  return(tM_j1)
}
