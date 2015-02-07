# Gaussian process everything

gpe is an R package implementing a Gaussian process (GP) statistical modelling framework. The aims of the project are to provide *some* of the functionality of fully featured GP toolboxes available in other languages, such as [GPy](https://github.com/SheffieldML/GPy) for python and [GPstuff](http://becs.aalto.fi/en/research/bayes/gpstuff/) for MATLAB.

Latent Gaussian process models constructed using compositional kernels - that is kernels created by adding, multiplying and convolving a series of basis kernels - can take a vast array of different functional forms, including many widely used statistical models as special cases. These include generalised linear (mixed, fixed or random) models; generalised additive models; geostatistical models; multivariate response models including joint species distribution models and (I think) model based ordination-type approaches. Under a GP approach, components of each of these models can be mixed and matched and extended in weird and wonderful different directions.

The aim of the package (and a planned accompanying demonstration repo) is to demonstrate how these models are all related and stimulate research into what can be done with these mdoels in an applied setting. It's certainly hoped that the package will developed into a fully-featured toolbox for all your GP needs, but that's a long way off for now. In the mean time you'd be better off using one of the existing packages mentioned above, or one of the many other packages that I haven't mentioned

## Progress

Not very much at the moment

![serious](https://pbs.twimg.com/media/B8XypnjIIAI333B.png)
(image via [@tom_speak on twitter](https://twitter.com/tom_speak/status/560120527311106050))


## Installation

You can install this package directly from github using the wonderful devtools package as follows:

```r
library(devtools)
install_github('goldingn/gpe')
```

but there's nothing here yet, so why would you want to do that?


# License

The package is distributed under the MIT license, which I think means you can do pretty much anything with it, but see the LICENSE file for details
