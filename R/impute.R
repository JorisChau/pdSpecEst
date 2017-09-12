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
        # tM_j1[, , i] <- Expm(P0, apply(array(rep(Nw[[Di + 1]], each = d^2),
        # dim = c(d, d, 2 * Di + 1)) * M_j0.log[, , (i - Di):(i + Di)], c(1, 2), sum))
      }
    }
  }
  return(tM_j1)
}

#' 2D AI refinement scheme
#'
#'
#'@export
Impute2D <- function(M_j0, D, cores = NULL, tol = NULL){

  # #'@importFrom parallel makeCluster
  # #'@importFrom doParallel registerDoParallel
  # #'@importFrom foreach foreach
  # #'@importFrom parallel stopCluster
  # #'@importFrom foreach "%dopar%"

  j0 <- log(dim(M_j0)[3], 2)
  d <- dim(M_j0)[1]

  if(identical(j0, 0)){

    tM_j1 <- array(M_j0[, , 1, 1], dim = c(d, d, 2, 2))

  } else {

    tM_j1 <- array(dim = c(d, d, 2^(j0 + 1), 2^(j0 + 1)))
    grid_k <- expand.grid(1:2^j0, 1:2^j0)
    grid_l <- expand.grid(1:(2 * D[1] + 1),1:(2 * D[2] + 1))

    # if(!is.null(cores)){
    #
    #   cl <- parallel::makeCluster(cores)
    #   doParallel::registerDoParallel(cl)
    #
    #   tM1 <- foreach::foreach(i=1:nrow(grid_k), .packages = "pdSpecEst") %dopar% {
    #
    #             tM_j1 <- array(dim = c(d, d, 2, 2))
    #
    #             ## Predict 2D Neville's algorithm
    #             N_k <- c(2 * min(grid_k[i, 1] - 1, D[1]) + 1, 2 * min(grid_k[i, 2] - 1, D[2]) + 1)
    #             nbrs <- list(x = max(grid_k[i, 1] - D[1], 1):min(grid_k[i, 1] + D[1], 2^j0),
    #                          y = max(grid_k[i, 2] - D[2], 1):min(grid_k[i, 2] + D[2], 2^j0))
    #             grid_l <- expand.grid(nbrs$x, nbrs$y)
    #             M.bar <- array(c(mapply(function(i1,i2) KarchMean(array(c(M_j0[, , grid_l[1, 1]:i1, grid_l[1, 2]:i2]),
    #                            dim = c(d, d, length(grid_l[1, 1]:i1) * length(grid_l[1, 2]:i2))), tol = tol), grid_l$Var1,
    #                            grid_l$Var2, SIMPLIFY="array")), dim = c(d, d, length(nbrs$x), length(nbrs$y)))
    #             M.bar_pred <- pdNeville2D(P = M.bar, ti = nbrs, grid = list(x = grid_k[i, 1] - 1 + c(0, 0.5, 1),
    #                                       y = grid_k[i, 2] - 1 + c(0, 0.5, 1)))
    #
    #             ## Predict M_{j+1,2k_1+1,2k_2+1}
    #             C0 <- KarchMean(array(c(M.bar_pred[, , 3, 2], M.bar_pred[, , 2, 2]), dim=c(d, d, 2)), w = c(N_k[1] + 1, -N_k[1]))
    #             M.bar_min0 <- KarchMean(array(c(M.bar_pred[, , 2, 3], C0), dim=c(d, d, 2)),
    #                                     w = c((N_k[1] * (N_k[2] + 1)) / ((N_k[1] + 1) * (N_k[2] + 1) - 1), N_k[2] / ((N_k[1] + 1) * (N_k[2] + 1) - 1)))
    #             tM_j1[, , 2, 2] <- KarchMean(array(c(M.bar[, , which(nbrs$x == grid_k[i, 1]),
    #                                     which(nbrs$y == grid_k[i, 2])], M.bar_min0), dim = c(d, d, 2)),
    #                                     w = c((N_k[1] + 1) * N_k[2] + N_k[1] + 1, -((N_k[1] + 1) * N_k[2] + N_k[1])))
    #
    #             ## Predict M_{j+1,2k_1,2k_2+1}
    #             C1 <- KarchMean(array(c(M.bar_pred[, , 1, 3], M.bar_pred[, , 1, 2]), dim = c(d, d, 2)), w = c(N_k[2] + 1, -N_k[2]))
    #             M.bar_min1 <- KarchMean(array(c(M.bar_pred[, , 3, 2], C1), dim=c(d, d, 2)),
    #                                     w = c((N_k[1] + 1) * N_k[2] / ((N_k[1] + 1) * (N_k[2] + 1) - 2), (N_k[1] - 1) / ((N_k[1] + 1) * (N_k[2] + 1) - 2)))
    #             M.block1 <- KarchMean(array(c(M.bar[, , which(nbrs$x == grid_k[i, 1]), which(nbrs$y == grid_k[i, 2])], M.bar_min1),
    #                                   dim = c(d, d, 2)), w = c(((N_k[1] + 1) * (N_k[2] + 1) - 2) / 2 + 1, -((N_k[1] + 1) * (N_k[2] + 1) - 2) / 2))
    #             tM_j1[, , 1, 2] <- KarchMean(array(c(M.block1, tM_j1[, , 2, 2]), dim = c(d, d, 2)), w = c(2, -1))
    #
    #             ## Predict M_{j+1,2k_1+1,2k_2}
    #             C2 <- KarchMean(array(c(M.bar_pred[, , 3, 1], M.bar_pred[, , 2, 1]), dim = c(d, d, 2)), w = c(N_k[1] + 1, -N_k[1]))
    #             M.bar_min2 <- KarchMean(array(c(M.bar_pred[, , 2, 3], C2), dim = c(d, d, 2)),
    #                                     w = c((N_k[2] + 1) * N_k[1] / ((N_k[2] + 1) * (N_k[1] + 1) - 2), (N_k[2] - 1) / ((N_k[2] + 1) * (N_k[1] + 1) - 2)))
    #             M.block2 <- KarchMean(array(c(M.bar[, , which(nbrs$x == grid_k[i, 1]), which(nbrs$y == grid_k[i, 2])], M.bar_min2),
    #                                   dim = c(d, d, 2)), w = c(((N_k[1] + 1) * (N_k[2] + 1) - 2) / 2 + 1, -((N_k[1] + 1) * (N_k[2] + 1) - 2) / 2))
    #             tM_j1[, , 2, 1] <- KarchMean(array(c(M.block2, tM_j1[, , 2, 2]), dim = c(d, d, 2)), w = c(2, -1))
    #
    #             ## Predict M_{j+1,2k_1,2k_2}
    #             tM_j1[, , 1, 1] <- Expm(M_j0[, , grid_k[i, 1], grid_k[i, 2]],
    #                                     -(Logm(M_j0[, , grid_k[i, 1], grid_k[i, 2]], tM_j1[, , 2, 1]) +
    #                                       Logm(M_j0[, , grid_k[i, 1], grid_k[i, 2]], tM_j1[, , 1, 2]) +
    #                                       Logm(M_j0[, , grid_k[i, 1], grid_k[i, 2]], tM_j1[, , 2, 2])))
    #
    #             return(tM_j1)
    #   }
    #
    #   parallel::stopCluster(cl)
    #
    #   for(l in 1:nrow(grid_k)){
    #     tM_j1[, , (2 * grid_k[l, 1] - 1):(2 * grid_k[l, 1]), (2 * grid_k[l, 2] - 1):(2 * grid_k[l, 2])] <- tM1[[l]]
    #   }
    #
    # } else{

      for(i in 1:nrow(grid_k)){

        ## Predict 2D Neville's algorithm
        N_k <- c(2 * min(grid_k[i, 1] - 1, D[1]) + 1, 2 * min(grid_k[i, 2] - 1, D[2]) + 1)
        nbrs <- list(x = max(grid_k[i, 1] - D[1], 1):min(grid_k[i, 1] + D[1], 2^j0),
                     y = max(grid_k[i, 2] - D[2], 1):min(grid_k[i, 2] + D[2], 2^j0))
        grid_l <- expand.grid(nbrs$x, nbrs$y)
        M.bar <- array(c(mapply(function(i1,i2) KarchMean(array(c(M_j0[, , grid_l[1, 1]:i1, grid_l[1, 2]:i2]),
                                                                dim = c(d, d, length(grid_l[1, 1]:i1) * length(grid_l[1, 2]:i2))), tol = tol), grid_l$Var1,
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

    # }

  }

  return(tM_j1)

}

