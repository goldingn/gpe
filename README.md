# Gaussian process everything

gpe is an R package implementing a Gaussian process (GP) statistical modelling framework. The aims of the project are to provide *some* of the functionality of fully featured GP toolboxes available in other languages, such as [GPy](https://github.com/SheffieldML/GPy) for python and [GPstuff](http://becs.aalto.fi/en/research/bayes/gpstuff/) for MATLAB. Mostly though it's about communicating the power of the GP approach to an applied statistical audience.

## Why?

Latent Gaussian process models constructed using compositional kernels - that is kernels created by adding, multiplying and convolving a series of basis kernels - can take a vast array of different functional forms, including many widely used statistical models as special cases. These include generalised linear (mixed, fixed or random) models; generalised additive models; geostatistical models; multivariate response models including joint species distribution models and (I think) model based ordination-type approaches. Under a GP approach, components of each of these models can be mixed and matched and extended in weird and wonderful directions.

The aim of the package (and a planned accompanying demonstration repo) is to demonstrate how these models are all related and stimulate research into what can be done with the GP framework in an applied setting. It's certainly hoped that the package will developed into a fully-featured toolbox for all your GP needs, but that's a long way off for now. In the mean time you'd be better off using one of the existing packages mentioned above, or one of the many other packages that I haven't mentioned.

## Installation

You can install this package directly from github using the wonderful devtools package as follows:

```r
library(devtools)
install_github('goldingn/gpe')
```

## Progress

### Kernels

There are a small handful of basis kernels currently implemented. These are easy to construct and visualise, like this:

```r
# Create an rbf kernel which acts on some variable named temperature
k1 <- rbf('temperature')
# look at the parameters
summary(k1)
# plot covariance
plot(k1)
# look at some GPs drawn from this kernel
demoKernel(k1)
```

Crucially though, these kernels may be combined to create new functional forms:

```r
# create a linear kernel
k2 <- lin('temperature')
# change a parameter
k2 <- setParameters(k2, sigma = 0.5)
# plot draws from it
demoKernel(k2)

# add the two together
k3 <- k1 + k2
# and visualise them
demoKernel(k3)

# multiply this by a periodic kernel
k4 <- k3 * per('temperature')
# visualise
plot(k4)
demoKernel(k4)
```

### GPs

There is now an interface for fitting GP models to Gaussian and some non-Gaussian (Poisson and Bernoulli at the moment) data, and a function to make predictions from these models, though functions to summarise their output and learn the kernel parameters are yet to be implemented.

## License

The package is distributed under the MIT license, which I think means you can do pretty much anything with it, but see the LICENSE file for details
