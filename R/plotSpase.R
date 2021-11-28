#' Main function for plotting ASE fits from spase fit
#'
#' Takes in the result from spase and plots raw data and smoothed allele
#' probability surface
#'
#' @description This function plots stuff for spase results.
#'
#' @param matrix1 a matrix of counts where  the rows are genes and the columns
#' are cells for allele 1. This one is the one that gets its probability
#' modeled. Must  have  row  names and column names to identify genes and cells.
#'
#' @param matrix2 a matrix of counts where the rows are genes and the columns
#' are cells for allele 2.
#'
#' @param covariates data frame that matches pixel IDs to spatial coordinates
#' and other covariates. First column assumed to be pixel IDs, next two columns
#' assumed to be x,y coordinates. Any additional columns assumed to be
#' covariates that will be included in the baseline model.
#'
#' @param spasefit the resulting data frame from a call to scase
#'
#' @param coords data frame of coordinates to evaluate the model on. Default is
#' NULL which will plot coordinates along the range of the x1 and x2 coordinate
#' values. Highly suggest supplying a custom set of coordinates since the
#' default probably doesn't work well for all data sets.
#'
#' @param extrapolate.width numeric which defines the neighborhood around each
#' coordinate from coords to look for an input data point. If there is no
#' point in the neighborhood,
#'
#' @param crosshairs boolean default of TRUE indicating plotting the cross-
#' hairs on the 2D spatial plot and the associated confidence intervals of
#' those slices in a separate plot. Set to FALSE to turn off.
#'
#' @param size.scale boolean default of TRUE indicating the raw data to be
#' plotted with points that change size depending on the total UMI count at
#' that pixel. Set to FALSE to turn off.
#'
#' @param point.size numeric to set the point size, default is 0.5. Must be
#' used with setting size.scale to FALSE.
#'
#' @param point.outline boolean to determine whether or not to draw points on
#' raw data with a black outline. default is TRUE. Must be used with
#' setting size.scaleto FALSE.
#'
#' @param theme string to set what theme to use for plotting smooth 2D
#' probability functions, default is the default ggplot theme, despite ugly
#' gray background, since the scale uses white as a color. If you're not
#' worried about that you can set to 'classic' to remove gray background
#' and grid lines.
#'
#' @param cross.x1 numeric indicating which x1 value to make crosshairs at.
#' default is NULL which, if crosshairs is TRUE, defaults to the middle of the
#' x1 range.
#'
#' @param cross.x2 numeric indicating which x2 value to make crosshairs at.
#' default is NULL which, if crosshairs is TRUE, defaults to the middle of the
#' x2 range.
#'
#' @param genes a vector of genes given as strings to plot. default is NULL
#' which will plot the first gene in the spasefit results dataframe.
#'
#' @param save the filename prefix to save if want to save instead of visualize.
#' Default is null meaning no file will be saved.
#'
#'
#' @return Multiple ggplot2 plots or, instead, each one gets saved to a different
#' file starting with the prefix defined in the save parameter.
#'
#' @import ggplot2
#' @importFrom mgcv smoothCon s Predict.matrix
#' @importFrom reshape2 melt
#' @importFrom RColorBrewer brewer.pal
#' @import dplyr
#' @import latex2exp
#'
#' @export

plotSpase <- function(matrix1, matrix2, covariates, spasefit, coords=NULL,
                      extrapolate.width=0.1, crosshairs=TRUE, size.scale=TRUE,
                      point.size = 0.5, point.outline = TRUE, theme='default',
                      cross.x1=NULL, cross.x2=NULL, genes=NULL, save=NULL) {
  if(!(is(matrix1, 'matrix') | is(matrix1, 'sparseMatrix')) &
     !(is(matrix2, 'matrix') | is(matrix2, 'sparseMatrix'))) {
    stop('input matrix1 and matrix2 must be a matrix or sparseMatrix')
  }
  if (ncol(matrix1)!=ncol(matrix2)) {
    stop('matrices dont have same number of columns')
  }
  if (nrow(matrix1) != nrow(matrix2)) {
    stop('matrices dont have same number of rows')
  }
  if (is.null(rownames(matrix1))) {
    warning('input matrix1 doesnt have rownames')
  }
  if (is.null(colnames(matrix1))) {
    warning('input matrix1 doesnt have colnames')
  }
  num.cov <- ncol(covariates)
  mydim <- 2
  if (num.cov<2) {
    stop('covariates should contain 1) pixel names, 2) and possibly 3), spatial coordinates')
  } else if (num.cov==2) {
    message(paste('found 2 columns in covariates; going to assume that first column is pixel names and second column is 1D coordinates'))
    colnames(covariates) <- c('pixel','x1')
    covariates.norm <- covariates
    covariates.norm$x1 <- (covariates$x1-spasefit$coord.mean)/spasefit$coord.sd
    sm <- smoothCon(s(x1,k=df,fx=T,bs='tp'),data=covariates.norm)[[1]]
    mydim <- 1
  } else if (num.cov==3) {
    message(paste('found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates'))
    colnames(covariates)[1:3] <- c('pixel','x1','x2')
    covariates.norm <- covariates
    covariates.norm[,2:3] <- as.data.frame(sweep(covariates[,2:3],2,spasefit$coord.mean,'-')/spasefit$coord.sd)
    sm <- smoothCon(s(x1,x2,k=spasefit$df,fx=T,bs='tp'),data=covariates.norm)[[1]]

  } else if (num.cov>3) {
    colnames(covariates)[1:3] <- c('pixel','x1','x2')
    covariates.norm <- covariates
    covariates.norm[,2:3] <- as.data.frame(sweep(covariates[,2:3],2,spasefit$coord.mean,'-')/spasefit$coord.sd)
    sm <- smoothCon(s(x1,x2,k=spasefit$df,fx=T,bs='tp'),data=covariates.norm)[[1]]
    baseline.cov <- colnames(covariates)[-c(1:3)]
    message(paste('using', paste(baseline.cov, collapse=','),
                  'as baseline covariates'))
    baseline.cov.classes <- sapply(covariates, 'class')[-c(1:3)]
  }
  if (is.null(genes)) {
    genes <- spasefit$result$gene[1]
  }
  zscore.palette <- rev(brewer.pal(n=7, 'RdBu'))
  for (i in 1:length(genes)) {
    print(genes[i])
    if (is.null(spasefit$fits[genes[i]][[1]])) {
      warning(paste(genes[i], 'is not present in the results, probably did not pass minimum threshold of UMI or pixels'))
      next
    }
    y <- matrix1[genes[i],]
    total <- y + matrix2[genes[i],]
    y <- y[total > 0]; total <- total[total > 0]
    n <- sum(total)
    covari <- left_join(data.frame(pixel=names(y)), covariates.norm, by='pixel')
    # accounting for pixels w/o coordinates:
    if (mydim > 1) {
      remove.idx <- which(is.na(covari$x1)|is.na(covari$x2))
    } else {
      remove.idx <- which(is.na(covari$x1))
    }
    if (length(remove.idx)>0) {
      warning(paste(length(remove.idx), 'pixels did not have matching coordinates'))
      y <- y[-remove.idx]; total<-total[-remove.idx]
      covari <- covari[-remove.idx,]
    }
    covari.non.norm <- left_join(data.frame(pixel=covari$pixel), covariates, by='pixel')
    if (is.null(coords)) {
      xcoords <- seq(min(covari$x1), max(covari$x1), length.out=50)
      ycoords <- seq(min(covari$x2), max(covari$x2), length.out=50)
      pred.coords<- expand.grid(xcoords, ycoords)
    } else {
      pred.coords <- coords
    }
    df <- spasefit$df
    if (mydim > 1) {
      colnames(pred.coords) <- c('x1','x2')

    } else {
      colnames(pred.coords) <- c('x1')
    }
    smooth.model <- spasefit$fits[[genes[i]]]
    smooth.coef <- coef(smooth.model)
    pred.coords.norm <- sweep(pred.coords, 2, spasefit$coord.mean, '-')
    pred.coords.norm <- pred.coords.norm/spasefit$coord.sd
    Xp <- Predict.matrix(sm, pred.coords.norm)
    Xp <- cbind(Xp[,(df-2):df], Xp[,1:(df-3)])
    Xp[,2:df] <- sweep(Xp[,2:df], 2, spasefit$spline.mean, '-')
    Xp[,2:df] <- sweep(Xp[,2:df], 2, spasefit$spline.sd, '/')
    pred.coords$fit <- Xp%*%smooth.coef
    pred.coords$se <- sqrt(diag(Xp%*%vcov(smooth.model)%*%t(Xp)))
    pred.coords$z <- (pred.coords$fit/pred.coords$se)
    pred.coords <- pred.coords %>%
      mutate(zBin = case_when(
        (z < 1.5) & (z > -1.5) ~ '-1.5 < z < 1.5',
        (z >= 1.5 ) & (z < 2) ~ '1.5 <= z < 2',
        (z <= -1.5 ) & (z > -2) ~ '-2 < z <= -1.5',
        (z >= 2 ) & (z < 3) ~ '2 <= z < 3',
        (z <= -2 ) & (z > -3) ~ '-3 < z <= -2',
        z >=3 ~ 'z >= 3',
        z <= -3 ~ 'z <= -3'
      ))  %>%
      mutate(fillColor = case_when(
        zBin == 'z <= -3' ~ zscore.palette[1],
        zBin == '-3 < z <= -2' ~ zscore.palette[2],
        zBin == '-2 < z <= -1.5' ~ zscore.palette[3],
        zBin == '-1.5 < z < 1.5' ~ zscore.palette[4],
        zBin == '1.5 <= z < 2' ~ zscore.palette[5],
        zBin == '2 <= z < 3' ~ zscore.palette[6],
        zBin == 'z >= 3' ~ zscore.palette[7]
      ))
    pred.coords$zBin <- factor(pred.coords$zBin,
                               levels = c('z <= -3', '-3 < z <= -2',
                                          '-2 < z <= -1.5', '-1.5 < z < 1.5',
                                          '1.5 <= z < 2', '2 <= z < 3',
                                          'z >= 3' ))
    pr <- ggplot(cbind(data.frame(y=y, total=total), covari.non.norm),
                 aes(x=x1,y=x2))
    if (size.scale) {
      pr <- pr +
        geom_point(aes(fill=y/total,size=total), shape=21, color='black') +
        scale_fill_gradient2(name='y/total', low='blue', mid='white', high='red',
                             midpoint=0.5, limits=c(0,1))
    } else {
      if (point.outline) {
        pr <- pr +
          geom_point(aes(fill=y/total), size=point.size, shape=21, color='black') +
          scale_fill_gradient2(name='y/total', low='blue', mid='white', high='red',
                               midpoint=0.5, limits=c(0,1))
      } else {
        pr <- pr +
          geom_point(aes(color=y/total), size=point.size) +
          scale_color_gradient2(name='y/total', low='blue', mid='white', high='red',
                               midpoint=0.5, limits=c(0,1))
      }
    }
    pr <- pr +
      theme_classic() +
      ggtitle(paste(genes[i],'raw'))
    ps <- ggplot(pred.coords, aes(x=x1, y=x2)) +
      geom_raster(aes(fill=exp(fit)/(1+exp(fit)))) +
      scale_fill_gradient2(name='fitted p', low='blue', mid='white', high='red',
                           midpoint=0.5, limits=c(0,1)) +
      ggtitle(paste(genes[i],'smooth'))
    pover <- ps +
      geom_point(data=cbind(data.frame(y=y, total=total), covari.non.norm),
                 aes(fill=y/total, size=total),shape=21,color='black') +
      scale_size_continuous(limits=c(1,60), breaks=c(1,20,40,60), range=c(1,3))
    pz <- ggplot(pred.coords, aes(x=x1, y=x2, fill=fillColor)) +
      geom_raster(interpolate=F) +
      scale_fill_identity() +
      ggtitle(paste(genes[i], 'z score surface'))
    if (theme=='classic') {
      ps <- ps + theme_classic()
      pover <- pover + theme_classic()
    }
    if (crosshairs) {
      if (is.null(cross.x1)) {
        cross.x1 <- round(mean(pred.coords$x1),2)
      }
      if (is.null(cross.x2)) {
        cross.x2 <- round(mean(pred.coords$x2),2)
      }
      ps <- ps +
        geom_hline(yintercept = cross.x1, lty='dashed', color='black') +
        geom_vline(xintercept = cross.x2, lty='dashed', color='black')
      pover <- pover +
        geom_hline(yintercept = cross.x1, lty='dashed', color='black') +
        geom_vline(xintercept = cross.x2, lty='dashed', color='black')
    }
    if (is.null(save)) {
      print(pr)
      print(ps)
      print(pover)
      print(pz)
    } else {
      savename <- paste0(save,'-raw-',genes[i],'.png')
      ggsave(file=savename, plot=pr, height=3, width=4)
      message(paste('saved', savename))
      savename <- paste0(save,'-smooth-',genes[i],'.png')
      ggsave(file=savename, plot=ps, height=3, width=4)
      message(paste('saved', savename))
      savename <- paste0(save,'-overlay-',genes[i],'.png')
      ggsave(file=savename, plot=pover, height=3, width=4)
      message(paste('saved', savename))
      savename <- paste0(save,'-zscore-',genes[i],'.png')
      ggsave(file=savename, plot=pz, height=3, width=3)
      message(paste('saved', savename))
    }
    if (crosshairs) {
      x1.pred <- data.frame(x1=cross.x1,
                            x2=seq(min(pred.coords$x2),max(pred.coords$x2),
                                   length.out=50))
      x2.pred <- data.frame(x1=seq(min(pred.coords$x1),max(pred.coords$x1),
                                   length.out=50),
                            x2=cross.x2)
      x1.pred.norm <- sweep(x1.pred, 2, spasefit$coord.mean, '-')/spasefit$coord.sd
      x2.pred.norm <- sweep(x2.pred, 2, spasefit$coord.mean, '-')/spasefit$coord.sd
      Xp <- Predict.matrix(sm, x1.pred.norm)
      Xp <- cbind(Xp[,(df-2):df], Xp[,1:(df-3)])
      Xp[,2:df] <- sweep(Xp[,2:df], 2, spasefit$spline.mean, '-')
      Xp[,2:df] <- sweep(Xp[,2:df], 2, spasefit$spline.sd, '/')
      se <- sqrt(diag(Xp%*%vcov(smooth.model)%*%t(Xp)))
      x1.pred$fit <- Xp%*%smooth.coef
      x1.pred$ci.low <- x1.pred$fit-qnorm(.975)*se
      x1.pred$ci.high <- x1.pred$fit+qnorm(.975)*se
      Xp <- Predict.matrix(sm, x2.pred.norm)
      Xp <- cbind(Xp[,(df-2):df], Xp[,1:(df-3)])
      Xp[,2:df] <- sweep(Xp[,2:df], 2, spasefit$spline.mean, '-')
      Xp[,2:df] <- sweep(Xp[,2:df], 2, spasefit$spline.sd, '/')
      se <- sqrt(diag(Xp%*%vcov(smooth.model)%*%t(Xp)))
      x2.pred$fit <- Xp%*%smooth.coef
      x2.pred$ci.low <- x2.pred$fit-qnorm(.975)*se
      x2.pred$ci.high <- x2.pred$fit+qnorm(.975)*se
      px1 <- ggplot(x1.pred, aes(x=x2, y=exp(fit)/(1+exp(fit)))) +
        geom_line() +
        geom_ribbon(aes(ymin=exp(ci.low)/(1+exp(ci.low)),
                        ymax=exp(ci.high)/(1+exp(ci.high))),
                    lty='blank',alpha=0.25) +
        geom_hline(yintercept=0.5, lty='dashed') +
        ylim(c(0,1)) +
        xlab('x2 coordinate') +
        ylab('estimated p') +
        theme_classic() +
        ggtitle(paste(genes[i], 'x1 =', cross.x1, 'slice'))
      px2 <- ggplot(x2.pred, aes(x=x1, y=exp(fit)/(1+exp(fit)))) +
        geom_line() +
        geom_ribbon(aes(ymin=exp(ci.low)/(1+exp(ci.low)),
                        ymax=exp(ci.high)/(1+exp(ci.high))),
                    lty='blank',alpha=0.25) +
        geom_hline(yintercept=0.5, lty='dashed') +
        ylim(c(0,1)) +
        xlab('x1 coordinate') +
        ylab('estimated p') +
        theme_classic() +
        ggtitle(paste(genes[i], 'x2 =', cross.x2, 'slice'))
      if (is.null(save)) {
        print(px1)
        print(px2)
      } else {
        savename <- paste0(save,'-x1-',genes[i],'.png')
        ggsave(file=savename, plot=px1, height=2, width=2)
        message(paste('saved', savename))
        savename <- paste0(save,'-x2-',genes[i],'.png')
        ggsave(file=savename, plot=px2, height=2, width=2)
        message(paste('saved', savename))
      }
    }
  }
}
