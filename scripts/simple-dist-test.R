devtools::load_all()

# Create data for a circular setup with three matrices and three
# object sets of varying sizes.

set.seed(2837483)

ps <- c(100, 100, 100)
inds <- matrix(
  c(
    1, 2,
    1, 3,
    2, 3
  ),
  nrow = 3,
  ncol = 2,
  byrow = TRUE
)

K <- 4

alpha <- matrix(
  c(
    1000, 1000,    1,    1,
    1000,    1,    1, 1000,
    1000,    1, 1000,    1
  ),
  nrow = length(ps),
  ncol = K,
  byrow = TRUE
)

tau <- 0.5 # For Gaussian case, not applicable for other distributions

us <- lapply(seq_along(ps), function(i) {
  do.call(
    cbind, lapply(alpha[i, ], \(a) rnorm(ps[i], 0, sd = 1 / sqrt(a)))
  )
})

xs_true <- list(
  us[[1]] %*% t(us[[2]]),
  us[[1]] %*% t(us[[3]]),
  us[[2]] %*% t(us[[3]])
)

ys_true <- list(
  xs_true[[1]],
  matrix(
    sapply(1 / (1 + exp(-xs_true[[2]])), \(p) rbinom(1, 1, p)),
    nrow = ps[1],
    ncol = ps[3]
  ),
  matrix(
    sapply(log(1 + exp(xs_true[[3]])), \(r) rpois(1, r)),
    nrow = ps[2],
    ncol = ps[3]
  )
)

# For Gaussian case only
noise <- matrix(
  rnorm(prod(dim(ys_true[[1]])), sd = 1 / sqrt(tau)),
  nrow = ps[1],
  ncol = ps[2]
)

xs <- list(ys_true[[1]] + noise, ys_true[[2]], ys_true[[3]])

# Convert the data into the right format
triplets <- lapply(xs, matrix_to_triplets)
triplets_true <- lapply(xs_true, matrix_to_triplets)

# Missing entries correspond to missing rows in the triple representation
# so they can be removed from training data by simply taking a subset
# of the rows.
train <- list()
test <- list()

train_true <- list()
test_true <- list()

keep_for_training <- 0.75 * sapply(xs, \(x) prod(dim(x)))
for (m in seq_along(triplets)) {
  subset <- sample(nrow(triplets[[m]]), keep_for_training[m])
  train[[m]] <- triplets[[m]][subset, ]
  test[[m]] <- triplets[[m]][-subset, ]

  train_true[[m]] <- triplets_true[[m]][subset, ]
  test_true[[m]] <- triplets_true[[m]][-subset, ]
}

# Learn the model with the correct likelihoods
K <- 4
likelihood <- c("gaussian", "bernoulli", "poisson")
opts <- getCMFopts()
opts$iter.max <- 20000 # Less iterations for faster computation
opts$grad.iter <- 5
opts$grad.reg <- 0.7
opts$verbose <- 3
opts$prior.alpha_0 <- 1
opts$prior.beta_0 <- 1e-3
model <- CMF(train, inds, K, likelihood, ps, test = test, opts = opts)

plot(model$cost)

# Check the predictions
print(test[[1]][1:3, ])
print(model$out[[1]][1:3, ])

print(test[[2]][1:20, ])
print(model$out[[2]][1:20, ])
