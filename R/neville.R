#' Neville's algorithm on the manifold
#'
#' @export
pdNeville <- function(P, ti, grid){

  n <- dim(P)[3] - 1
  d <- dim(P)[1]
  geo <- function(Pi, t.range, t)  Expm(Pi[, , 1], (t - t.range[1]) / (t.range[2] - t.range[1]) * Logm(Pi[, , 1], Pi[, , 2]))

  Li <- sapply(1:length(grid), function(i2) sapply(1:n, function(i1) geo(P[, , c(i1, i1 + 1)],
                                                                         ti[c(i1, i1 + 1)], grid[i2]), simplify = "array"), simplify = "array")

  if(n > 1){
    for(r in 2:n){
      Li_new <- sapply(1:length(grid), function(i2) sapply(1:(n-r+1), function(i1) geo(Li[, , c(i1, i1 + 1), i2],
                                                                                       ti[c(i1, i1 + r)], grid[i2]), simplify = "array"), simplify = "array")
      Li <- Li_new
    }
  }

  return(Li[,,1,])
}

#' 2D Neville's algorithm on the manifold
#'
#' @export
pdNeville2D <- function(P, ti, grid){

  d <- dim(P)[1]
  if(missing(ti)){
    ti <- grid
  }

  if(identical(dim(P)[3], as.integer(1)) & identical(dim(P)[4], as.integer(1))){

    P_new <- array(P[, , 1, 1], dim = c(d, d, length(grid$x), length(grid$y)))

  } else if(identical(dim(P)[3], as.integer(1)) | identical(dim(P)[4], as.integer(1))){

    if(length(ti$y) > 1){
      P_i <- pdNeville(P[, , 1, ], ti$y, grid$y)
      P_new <- aperm(sapply(1:length(grid$x), function(i) P_i, simplify = "array"), c(1, 2, 4, 3))
    } else if(length(ti$x) > 1){
      P_j <- pdNeville(P[, , , 1], ti$x, grid$x)
      P_new <- sapply(1:length(grid$y), function(j) P_j, simplify = "array")
    }

  } else {

    P_i <- sapply(1:dim(P)[3], function(i) pdNeville(P[, , i, ], ti$y, grid$y), simplify = "array")
    P_new <- sapply(1:length(grid$y), function(i) pdNeville(P_i[, , i, ], ti$x, grid$x), simplify = "array")
  }

  return(P_new)

}
