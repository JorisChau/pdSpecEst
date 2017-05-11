#' Depth-rank hypothesis testing for HPD matrices
#'
#' @importFrom stats pnorm
#' @importFrom stats pchisq
#' @export
pdRankTests <- function(samples, sample.sizes, depth = c('gdd', 'zonoid'),
                        test = c('rank.sum', 'krusk.wall', 'signed.rank', 'bartels',
                                 'rank.sum.fun')){
  if(missing(depth)){
    depth <- 'gdd'
  }
  test <- match.arg(test, c('rank.sum', 'krusk.wall', 'signed.rank', 'bartels', 'rank.sum.fun'))
  depth <- match.arg(depth, c('gdd', 'zonoid'))
  n <- sample.sizes
  ddim <- dim(samples)

  if(test == 'rank.sum'){
    if (!isTRUE((length(n) == 2) & (length(ddim) == 3) & (ddim[3] == sum(n))))  {
      stop("Incorrect input lenghts for arguments: 'samples' and/or 'sample.sizes',
                  consult the function documentation for the requested inputs.")
    }
    if(depth == 'gdd'){
      dd <- pdDepth(X = samples, method = 'gdd')
    } else if(depth == 'zonoid'){
      dd <- sapply(1:sum(n), function(i) pdDepth(y = samples[,,i], X = samples, method = 'zonoid'))
    }
    T1 <- (sum(rank(dd, ties.method = 'random')[1:n[1]]) - n[1]*(sum(n)+1)/2) / sqrt(n[1]*n[2]*(sum(n)+1) / 12)

    output <- list(p.value = 2 * stats::pnorm(abs(T1), lower.tail = F), statistic = T1,
                   null.distr = 'Standard normal distribution')
  }

  if(test == 'rank.sum.fun'){
    if (!isTRUE((length(n) == 2) & (length(ddim) == 4) & (ddim[4] == sum(n)))) {
      stop("Incorrect input lengths for arguments: 'samples' and/or 'sample.sizes',
             consult the function documentation for the requested inputs.")
    }
      dd <- sapply(1:sum(n), function(i) pdDepth(y = samples[,,,i], X = samples,
                                                 method = ifelse(depth == 'gdd', 'gdd', 'zonoid')))

      T1 <- (sum(rank(dd, ties.method = 'random')[1:n[1]]) - n[1]*(sum(n)+1)/2) / sqrt(n[1]*n[2]*(sum(n)+1) / 12)

      output <- list(p.value = 2 * stats::pnorm(abs(T1), lower.tail = F), statistic = T1,
                     null.distr = 'Standard normal distribution')
  }

  if(test == 'krusk.wall'){
    N <- sum(n)
    if (!isTRUE((length(ddim) == 3) & (ddim[3] == N)))  {
      stop("Incorrect input lenghts for arguments: 'samples' and/or 'sample.sizes',
           consult the function documentation for the requested inputs.")
    }

    if(depth == 'gdd'){
      dd <- pdDepth(X = samples, method = 'gdd')
    } else if(depth == 'zonoid'){
      dd <- sapply(1:sum(n), function(i) pdDepth(y = samples[,,i], X = samples, method = 'zonoid'))
    }
    R_bar <- unname(unlist(lapply(split(rank(dd, ties.method = 'random'), f = rep(c(1,2,3), times = n)), mean)))
    T2 <- 12 / (N * (N+1)) * sum(n * (R_bar - (N+1)/2)^2)

    output <- list(p.value = min(stats::pchisq(T2, df = 2, lower.tail = T), pchisq(T2, df = 2, lower.tail = F)),
                statistic = T2, null.distr = "Chi-squared distribution (df = 2)")
  }

  if(test == 'signed.rank'){
    if (!isTRUE((length(ddim) == 3) & (ddim[3] == 2*n)))  {
      stop("Incorrect input lenghts for arguments: 'samples' and/or 'sample.sizes',
           consult the function documentation for the requested inputs.")
    }
    d <- dim(samples[,,1])[1]
    ast <- function(A,B) t(Conj(A)) %*% B %*% A
    diff <- sapply(1:n, function(i) Re(sum(diag(Logm(diag(d), ast(iSqrt(samples[,,n+i]), samples[,,i]))))))

    test <- stats::wilcox.test(x = diff, y = rep(0, n), mu = 0, paired = T)
    output <- list(p.value = test$p.value, statistic = test$statistic, null.distr = test$method)
  }

  if(test == 'bartels'){
    if (!isTRUE((length(ddim) == 3) & (ddim[3] == n)))  {
      stop("Incorrect input lenghts for arguments: 'samples' and/or 'sample.sizes',
           consult the function documentation for the requested inputs.")
    }
    if(depth == 'gdd'){
      dd <- pdDepth(X = samples, method = 'gdd')
    } else if(depth == 'zonoid'){
      dd <- sapply(1:n, function(i) pdDepth(y = samples[,,i], X = samples, method = 'zonoid'))
    }

    T4 <- sum(diff(rank(dd, ties.method = 'random'))^2) / (n*(n^2-1)/12)
    sigma <- sqrt(4*(n-2)*(5*n^2 - 2*n - 9) / (5*n*(n+1)*(n-1)^2))

    output <- list(p.value = 2 * pnorm(abs((T4-2)/sigma), lower.tail = F), statistic = (T4-2)/sigma,
                   null.distr = 'Standard normal distribution')
  }
  return(output)
}

