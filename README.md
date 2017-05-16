# HurdleNormal
Estimation and sampling for multi-dimensional Hurdle models on a Normal density with applications to single-cell co-expression

A Hurdle model is a modification of a distribution at zero. 
This package provides routines to estimate and sample from multivariate Hurdle models on a Normal density. 
These are distributions that are conditionally Normal, but with singularities along the coordinate axes, so generalize a univariate zero-inflated distribution.  
The main functionality is to estimate graphical models using group-lasso penalized neighborhood selection.

## Vignette
A vignette with an example is located
[here](http://htmlpreview.github.io/?https://github.com/amcdavid/HurdleNormal/blob/master/inst/doc/singlecell-networks.html)
or with `vignette('singlecell-networks')`.

## Citation
[Graphical Models for Zero-Inflated Single Cell Gene Expression](https://arxiv.org/abs/1610.05857)
