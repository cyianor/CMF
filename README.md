# Collective Matrix Factorization

Collective matrix factorization (CMF) finds joint low-rank representations
for a collection of matrices with shared row or column entities. This code
learns a variational Bayesian approximation for CMF, supporting multiple
likelihood potentials and missing data, while identifying both factors shared
by multiple matrices and factors private for each matrix.

For further details on the method see
[Klami et al. (2014)](https://arxiv.org/abs/1312.5921). The package can also
be used to learn Bayesian canonical correlation analysis (CCA) and
group factor analysis (GFA) models, both of which are special cases of CMF.
This is likely to be useful for people looking for CCA and GFA solutions
supporting missing data and non-Gaussian likelihoods.

See [Klami et al. (2013)](http://www.jmlr.org/papers/v14/klami13a.html) and
[Virtanen et al. (2012)](http://proceedings.mlr.press/v22/virtanen12.html)
for details on Bayesian CCA and GFA, respectively.

- **Original package authors**: Arto Klami and Lauri Väre
- **Maintainer**: Felix Held