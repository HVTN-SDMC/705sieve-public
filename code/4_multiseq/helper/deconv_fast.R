# ChatGPT modification of Efron's package to make it a little faster

deconv_fast <- function(tau, X, y, Q, P, n = 40,
                   family = c("Poisson", "Normal", "Binomial"),
                   ignoreZero = TRUE, deltaAt = NULL, c0 = 1,
                   scale = TRUE, pDegree = 5, aStart = 1.0, ...) {
  
  family <- match.arg(family)
  if (!is.null(deltaAt) && (family != "Normal")) {
    stop("deltaAt applies only to Normal.")
  }
  
  # ---------- Build P, Q, y (vectorized) ----------
  if (missing(Q) && missing(P)) {
    m <- length(tau)
    
    if (family == "Poisson") {
      supportOfX <- if (ignoreZero) seq_len(n) else seq.int(0, n - 1)
      # P: n x m
      P <- outer(supportOfX, tau, function(x, lam) stats::dpois(x, lam))
      if (ignoreZero) P <- sweep(P, 2, 1 - exp(-tau), "/")
      if (missing(y)) {
        # fast histogram of counts
        y <- tabulate(match(X, supportOfX), nbins = length(supportOfX))
      }
      Q1 <- splines::ns(tau, pDegree)
      if (scale) {
        Q1 <- scale(Q1, center = TRUE, scale = FALSE)
        Q1 <- sweep(Q1, 2, sqrt(colSums(Q1 * Q1)), "/")
      }
      Q <- cbind(1.0, Q1) # include intercept
      Q <- sweep(Q, 2, sqrt(colSums(Q * Q)), "/") # normalize columns
      
    } else if (family == "Normal") {
      # bins for z-scores
      r <- round(range(X), 1)
      xBin <- seq(from = r[1], to = r[2], length.out = n)
      lo <- xBin[-length(xBin)]
      hi <- xBin[-1]
      
      # P: (n-1) x m with Pr(bin | theta)
      P <- outer(lo, tau, stats::pnorm) - outer(hi, tau, stats::pnorm)
      
      intervals <- findInterval(X, xBin, all.inside = TRUE)
      y <- tabulate(pmax.int(1L, pmin.int(intervals, n - 1L)), nbins = n - 1L)
      
      Q1 <- splines::ns(tau, pDegree)
      if (scale) {
        Q1 <- scale(Q1, center = TRUE, scale = FALSE)
        Q1 <- sweep(Q1, 2, sqrt(colSums(Q1 * Q1)), "/")
      }
      if (!is.null(deltaAt)) {
        I0 <- as.numeric(abs(tau - deltaAt) < 1e-10)
        Q <- cbind(I0, Q1)
      } else {
        Q <- Q1
      }
      
    } else { # Binomial
      if (is.data.frame(X)) X <- as.matrix(X)
      if (!is.matrix(X) || ncol(X) != 2) {
        stop("For family='Binomial', X must be two columns: (n_i, s_i).")
      }
      storage.mode(X) <- "double"  # ensure numeric
      
      # Build P: nObs x m, each column j uses prob = tau[j]
      # (vapply is fast and clean; avoids outer’s overhead)
      P <- vapply(tau,
                  function(w) stats::dbinom(X[, 2], size = X[, 1], prob = w),
                  numeric(nrow(X)))
      y <- 1  # sentinel as before
      
      # Q as before
      Q <- splines::ns(tau, pDegree)
      if (scale) {
        Q <- scale(Q, center = TRUE, scale = FALSE)
        Q <- sweep(Q, 2, sqrt(colSums(Q * Q)), "/")
      }
      
    }
  } else {
    if (!missing(X) || missing(y) || missing(P) || missing(Q)) {
      stop("If supplying P/Q, pass (P, Q, y) together and omit X.")
    }
  }
  
  p <- ncol(Q)
  if (length(aStart) == 1L) aStart <- rep(aStart, p)
  if (length(aStart) != p) stop("aStart has wrong length.")
  
  softmax <- function(z) {
    zmax <- max(z); ez <- exp(z - zmax); ez / sum(ez)
  }
  
  # unified y vector (Binomial uses y=1 sentinel => yVec = 1s)
  yVec_make <- function(f, y, nrowP) {
    if (length(y) == 1L && y == 1) rep(1.0, nrowP) else as.numeric(y)
  }
  
  # ---------- Objective with fast gradient ----------
  loglik <- function(a) {
    Qa <- as.vector(Q %*% a)
    g  <- softmax(Qa)              # m
    f  <- as.vector(P %*% g)       # n
    f  <- pmax(f, 1e-300)          # guard
    yv <- yVec_make(f, y, nrow(P)) # n
    
    val <- -sum(yv * log(f)) + c0 * sqrt(sum(a * a))
    
    # Gradient: - t(Q) %*% ( g * (t(P) %*% (y/f) - sum(y)) ) + penalty'
    u   <- as.vector(t(P) %*% (yv / f))  # m
    sy  <- sum(yv)
    v   <- g * (u - sy)                  # m
    grad <- - as.vector(crossprod(Q, v)) + (c0 * a / sqrt(sum(a * a)))
    
    attr(val, "gradient") <- grad
    val
  }
  
  # ---------- Optimize (still nlm; you can switch to optim/nlminb if you like) ----------
  opt <- stats::nlm(f = loglik, p = aStart, gradtol = 1e-10, ...)
  mle <- opt$estimate
  
  # ---------- Stats (do heavy work once, but avoid the big D) ----------
  statsFunction <- function(a) {
    softmax <- function(z) { zmax <- max(z); ez <- exp(z - zmax); ez/sum(ez) }
    
    g  <- softmax(as.vector(Q %*% a))      # m
    G  <- cumsum(g)
    f  <- as.vector(P %*% g)               # n
    f  <- pmax(f, 1e-300)
    
    # match original logic
    yHat <- if (length(y) == 1L && y == 1) rep(1.0, length(f)) else sum(y) * f
    
    # --- core matrices ---
    Pt <- sweep(P, 1, f, "/")              # n x m   (row-wise scale)
    # W = g * (t(P/f)) - g * 1^T   -> m x n
    W  <- sweep(t(Pt), 1, g, "*") - g %o% rep(1.0, nrow(P))   # m x n
    
    # Fisher info pieces
    qw  <- crossprod(Q, W)                  # p x n   = t(Q) %*% W
    # (yHat * t(W)) %*% Q  -> n x p
    ywq <- (sweep(t(W), 1, yHat, "*")) %*% Q
    I1  <- qw %*% ywq                       # p x p
    
    # penalty terms
    aa       <- sqrt(sum(a * a))
    sDot     <- c0 * a / aa
    sDotDot  <- (c0 / aa) * (diag(length(a)) - (a %o% a) / (aa * aa))
    
    # R value
    Rval <- sum(diag(sDotDot)) / max(1e-300, sum(diag(I1)))
    
    # (I1 + sDotDot)^(-1) via Cholesky
    A  <- I1 + sDotDot
    R  <- tryCatch(chol(A), error = function(e) NULL)
    if (is.null(R)) { A <- A + diag(1e-8, nrow(A)); R <- chol(A) }
    I2 <- chol2inv(R)
    
    bias   <- as.vector(-I2 %*% sDot)
    Cov    <- I2 %*% (I1 %*% t(I2))
    
    Dq     <- (diag(g) - (g %o% g)) %*% Q
    bias.g <- as.vector(Dq %*% bias)
    Cov.g  <- Dq %*% Cov %*% t(Dq)
    se.g   <- sqrt(pmax(0, diag(Cov.g)))
    
    # Cov(G) via prefix sums (no giant D)
    CovG <- apply(Cov.g, 2, cumsum)
    CovG <- t(apply(CovG, 1, cumsum))
    se.G <- sqrt(pmax(0, diag(CovG)))
    
    mat <- cbind(theta = tau, g = g, SE.g = se.g, G = G, SE.G = se.G, Bias.g = bias.g)
    list(S = Rval, cov = Cov, cov.g = Cov.g, mat = mat)
  }
  stats <- statsFunction(mle)
  
  if (family == "Poisson" && ignoreZero) {
    mat <- stats$mat
    tg  <- mat[, "g"] / (1 - exp(-tau))
    stats$mat <- cbind(mat, tg = tg / sum(tg))
  }
  
  list(
    mle = mle,
    Q = Q,
    P = P,
    S = stats$S,
    cov = stats$cov,
    cov.g = stats$cov.g,
    stats = stats$mat,
    loglik = loglik,
    statsFunction = statsFunction
  )
}
