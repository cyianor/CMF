devtools::load_all()

ps <- c(30, 50, 60)
inds <- matrix(
  c(
    1L, 2L,
    1L, 3L
  ),
  nrow = 2,
  ncol = 2,
  byrow = TRUE
)

K <- 4

alpha <- matrix(
  c(
    1, 1000,     0.1,    1,
    1, 1000,     0.5, 1000,
    1, 1000, 10000.0,    1
  ),
  nrow = length(ps),
  ncol = K,
  byrow = TRUE
)

tau <- c(5, 8)

us <- lapply(seq_along(ps), function(i) {
  do.call(
    cbind, lapply(alpha[i, ], \(alpha) rnorm(ps[i], 0, sd = 1 / sqrt(alpha)))
  )
})

xs_true <- list(
  us[[1]] %*% t(us[[2]]),
  us[[1]] %*% t(us[[3]])
)

noise <- purrr::map2(xs_true, tau, ~ {
  matrix(
    rnorm(prod(dim(.x)), sd = 1 / sqrt(.y)),
    nrow = nrow(.x),
    ncol = ncol(.x)
  )
})

xs <- purrr::map2(xs_true, noise, ~ .x + .y)

# Convert the data into the right format
triplets <- lapply(xs, matrix_to_triplets)

# Missing entries correspond to missing rows in the triple representation
# so they can be removed from training data by simply taking a subset
# of the rows.
train <- list()
test <- list()
keep_for_training <- 0.75 * sapply(xs, \(x) prod(dim(x)))
for (m in seq_along(triplets)) {
  subset <- sample(nrow(triplets[[m]]), keep_for_training[m])
  train[[m]] <- triplets[[m]][subset, ]
  test[[m]] <- triplets[[m]][-subset, ]
}

# Learn the model with the correct likelihoods
K <- 4
likelihood <- c("gaussian", "gaussian")
opts <- getCMFopts()
opts$iter.max <- 10000 # Less iterations for faster computation
opts$grad.iter <- 1
opts$grad.reg <- 0.7
opts$verbose <- 3
opts$init.tau <- 10
opts$init.alpha <- 5
opts$prior.alpha_0 <- 1
opts$prior.beta_0 <- 1e-4
model <- CMF(train, inds, K, likelihood, ps, test = test, opts = opts)

# Check the predictions
print(test[[1]][1:3, ])
print(model$out[[1]][1:3, ])
