
#Linear behavior near zero (controlled by W)
#Logarithmic asymptotics for large positive values (controlled by M)
#Smooth handling of negative values (controlled by A)
#Numerically stable but slower due to uniroot per element



logicle_transform <- function(x, T = 262144, min_val = NULL, W = 0.5, M = 4.5, A = 0, tol = 1e-12, max_iter = 50) {
  
  ln10 <- log(10)
  
  # --- Solve for internal parameters ---
  # Following Moore/Parks/Robinson formulation
  
  b <- (M + A) * ln10
  w <- W / (M + A)
  
  # Solve for d numerically (core step often hidden in libraries)
  f_d <- function(d) {
    2 * (log(d) - log(b)) + w * (b + d)
  }
  
  d <- uniroot(f_d, c(1e-6, b))$root
  
  a <- T / (exp(b) - exp(-d))
  c <- a * exp(b * w) - a * exp(-d * w)
  
  # --- Newton solver ---
  y <- rep(0, length(x))
  
  # initial guess: hybrid linear/log
  y[x > 0] <- log10(1 + x[x > 0]) / M
  y[x < 0] <- -log10(1 + abs(x[x < 0])) / M
  
  for (iter in seq_len(max_iter)) {
    
    exp_by <- exp(b * y)
    exp_dy <- exp(-d * y)
    
    f  <- a * exp_by - a * exp_dy - c - x
    df <- a * b * exp_by + a * d * exp_dy
    
    step <- f / df
    y <- y - step
    
    if (max(abs(step)) < tol) break
  }
  
  # --- Normalize to display scale [-A, M] ---
  
  if(!is.null(min_val)){
      y <- ifelse(y < min_val, min_val, y)
  }
  y
}