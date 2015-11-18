# covarmat_class
# class for evaluated covariance matrices

# To do:
#   fix the colour scheme in image

#' @name covarmat
#' @rdname covarmat
#'
#' @title gpe covariance matrix class
#' 
#' @description Generic functions associated with the covarmat class.
#' 
#' @param x an object of class \code{covarmat}, or an object to be tested as one.
#' 
#' @details \code{is.covarmat} returns a logical indicating whether the 
#' object is a gpe \code{covarmat} object and \code{image.covarmat} provides 
#' a nice plot of the values in the covariance matrix; in the style of, and 
#' depending on \code{\link{image.plot}}. \code{palPuRd} uses the PuRd color
#' palette from \code{\link{RColorBrewer}} to produce a non-divergent colour
#' palette for visualising these matrices.
#' 
#' @examples
#'  
#' # construct a kernel with one feature
#' k1 <- rbf('a')
#' 
#' # create a fake dataframe
#' df <- data.frame(a = sort(runif(50, -3, 3)))
#'  
#' # evaluate the kernel
#' K <- k1(df)
#' 
NULL

# covariance matrix constructor
as.covarmat <- function (x) {
  class(x) <- c('covarmat')
  return (x)
}

#' @rdname covarmat
#' @export
#' @examples
#' 
#' # is it a covarmat object? 
#' is.covarmat(K)
#'  
is.covarmat <- function (x) {
  'covarmat' %in% class(x)
}

#' @rdname covarmat
#' 
#' @param axes whether to add axes to the image plot of \code{x}.
#' 
#' @param legend whether to add a legend to the image plot of \code{x}.
#'
#' @param col colour table used to plot the image
#' 
#' @param \dots other arguments to be passed to lower-level functions, such as 
#' image.
#' 
#' @export
#' @importFrom fields image.plot
#' 
#' @examples
#' 
#' # visualise it 
#' image(K)
#' 
#' # create a new dataframe
#' df2 <- data.frame(a = sort(runif(100, -3, 3)))
#' 
#' # visualise the covariance matrix between the two dataframes
#' image(k1(df, df2))
#' 
image.covarmat <- function (x,
                            axes = TRUE,
                            legend = TRUE,
                            col = palPuRd(100),
                            ...) {
  # define a plot function for covariance matrices
  
  # set nlevel to the same length as col
  # so that matrix and legend agree
  nlevel <- length(col)
  
  covar <- x
  
  # dimension
  nx <- ncol(covar)
  ny <- nrow(covar)
  
  # flip the matrix round (image does this for some reason)
  z <- as.matrix(t(x[ny:1, ]))
  
  # axis numbers
  x <- 1:nx
  y <- 1:ny
  
  # use image plot, make the legend plot fail silently
  tmp <- tryCatch(fields::image.plot(x = x,
                                     y = y,
                                     z = z,
                                     useRaster = TRUE,
                                     ylab = '',
                                     xlab = '',
                                     axes = FALSE,
                                     ...,
                                     legend.shrink = NA,
                                     nlevel = nlevel,
                                     col = col),
                  error = function(e) NULL)
  
  if (axes) {
    
    # get the axis points at which to print matrix indices
    by <- min(nx, ny) %/% 5
    at_x <- seq(0, nx, by = by)
    at_y <- seq(0, ny, by = by)
    
    # add axes
    axis(3,
         at_x[-1],
         at_x[-1],
         tick = FALSE,
         cex.axis = 0.8,
         line = -0.5)
    
    axis(2,
         at_y[-length(at_y)] + 1,
         rev(at_y[-1]),
         tick = FALSE,
         cex.axis = 0.8,
         line = -0.5,
         las = 2)
    
    # axis labels
    mtext('K_j',
          side = 3,
          line = 2,
          at = nx / 2)
    
    mtext('K_i',
          side = 2,
          line = 2,
          at = ny / 2,
          las = 2)
    
  }
  
  if (legend) {
    
    # check there is more than one unique covariance
    if (length(unique(as.vector(z))) > 1) {
      # add the legend
      fields::image.plot(x = x,
                         y = y,
                         z = z,
                         useRaster = TRUE,
                         ylab = '',
                         xlab = '',
                         axes = FALSE,
                         ...,
                         legend.only = TRUE,
                         add = TRUE,
                         nlevel = nlevel,
                         col = col)
      
    } else {
      # otherwise, get the value
      val <- mean(z, na.rm = TRUE)
      
      # and add a single colour legend
      fields::image.plot(covar,
                         legend.only = TRUE,
                         col = rep(col[1], 8),
                         smallplot = c(0.85, 0.9, 0.5, 0.55),
                         lab.breaks = rep(c('', val, ''),
                                          c(4, 1, 4)),
                         axis.args = list(tcl = 0))

    }
    
  }
}



# colour palettes
# purple to red, skipping the lightest tones

#' @rdname covarmat
#' 
#' @param n length of the continuous colour vector required
#' 
#' @export
#' @examples
#' 
#' # visualise the covariance matrix between these dataframes using
#' # the purple-red colour scheme (actually the default) but with fewer levels
#' image(k1(df, df2), col = palPuRd(10))
#' 
palPuRd <- function (n) {
  cols <- RColorBrewer::brewer.pal(9, 'PuRd')[-(1:2)]
  pal <- colorRampPalette(colors = cols)
  ans <- pal(n)
  return (ans)
}