#' Main function for plotting ASE fits in single-cell data
#'
#' Takes in the result from scase and plots desired genes with beta distribution
#' overlaid on the histogram of allele1 proportion, confidence interval plot
#' for the desired genes, and the overall distribution of p and phi for all
#' genes
#'
#' @description This function plots stuff for scase results.
#'
#' @param matrix1 a matrix of counts where  the rows are genes and the columns
#' are cells for allele 1. This one is the one that gets its probability
#' modeled. Must  have  row  names and column names to identify genes and cells.
#'
#' @param matrix2 a matrix of counts where the rows are genes and the columns
#' are cells for allele 2.
#'
#' @param scasefit the resulting data frame from a call to scase
#'
#' @param genes a vector of genes given as strings to plot. default is NULL
#' which will plot the first gene in the scasefit results dataframe.
#'
#'
#' @param save the filename prefix to save if want to save instead of visualize.
#' Default is null meaning no file will be saved.
#'
#' @param breaks how to split up the total UMI counts for coloring. there is
#' a default setting which should work for most data sets.
#'
#'
#' @return Multiple ggplot2 plots or, instead, each one gets saved to a different
#' file starting with the prefix defined in the save parameter.
#'
#' @import ggplot2
#' @importFrom ggExtra ggMarginal
#' @importFrom ggrepel geom_label_repel
#' @importFrom reshape2 melt
#' @import dplyr
#' @import viridis
#' @import latex2exp
#'
#' @export

plotScase <- function(matrix1, matrix2, scasefit, genes=NULL, save=NULL,
                      breaks=c(1,2,5,10,50,100,250,500,1000,2000)) {

  if(!is(matrix1, 'matrix') | !is(matrix2, 'matrix')) {
    stop('input matrix1 and matrix2 must be a matrix')
  }
  if (ncol(matrix1)!=ncol(matrix2)) {
    stop('matrices dont have same number of columns')
  }
  if (nrow(matrix1) != nrow(matrix2)) {
    stop('matrices dont have same number of rows')
  }
  if (is.null(rownames(matrix1))) {
    stop('input matrix1 doesnt have rownames; need rownames to be gene names')
  }
  if (is.null(colnames(matrix1))) {
    warning('input matrix1 doesnt have colnames')
  }

  matrix1[is.na(matrix2)] <- NA
  matrix2[is.na(matrix1)] <- NA
  if (is.null(genes)) {
    idx <- 1
  } else {
    idx <- which(rownames(matrix1)%in%genes)
  }
  if (length(idx)==1) {
    df.hist <- data.frame(gene = rownames(matrix1)[idx], cell = colnames(matrix1),
                          y = matrix1[idx,], total = matrix1[idx,]+matrix2[idx,]) %>%
      filter(total>0)
  } else {
    df.hist <- melt(matrix1[idx,], varnames=c('gene','cell'), value.name='y')
    df.hist <- left_join(
      df.hist,
      melt(matrix1[idx,]+matrix2[idx,],  varnames=c('gene','cell'),
           value.name='total')
    )  %>%
      filter(total>0)
  }

  df.pdf <- scasefit
  df.pdf$p <- exp(df.pdf$logit.p)/(1+exp(df.pdf$logit.p))
  df.pdf$p.upper <- exp(df.pdf$logit.p+2*df.pdf$logit.p.sd)/
    (1+exp(df.pdf$logit.p+2*df.pdf$logit.p.sd))
  df.pdf$p.lower <- exp(df.pdf$logit.p-2*df.pdf$logit.p.sd)/
    (1+exp(df.pdf$logit.p-2*df.pdf$logit.p.sd))
  df.pdf$a1 <- df.pdf$p*(1-df.pdf$phi)/df.pdf$phi
  df.pdf$a2 <- (1-df.pdf$p)*(1-df.pdf$phi)/df.pdf$phi
  df.pdf$var.p <- df.pdf$a1*df.pdf$a2/(df.pdf$a1+df.pdf$a2)^2/
    (df.pdf$a1+df.pdf$a2+1)

  mybreaks <- quantile(df.hist$total,probs=seq(0,1,0.2))

  df.hist$totalBin <- cut(df.hist$total,breaks=mybreaks,include.lowest=T)
  # color palette
  colornames <- levels(df.hist$totalBin)
  mycolors <- viridis(length(colornames))

  # first plot the raw data and overlaid beta pdf for each gene
  for (i in 1:length(genes)) {
    pidx <- which(df.pdf$gene == genes[i])
    a <- df.pdf$a1[pidx]
    b <- df.pdf$a2[pidx]
    pdfvals <- data.frame(x = seq(0.01,0.99,0.01),
                          density = dbeta(x=seq(0.01,0.99,0.01), shape1=a, shape2=b))
    df.hist.sub <- df.hist[which(df.hist$gene==genes[i]),]
    colors.subset <- mycolors[which(colornames%in%df.hist.sub$totalBin)]
    p <- ggplot(df.hist.sub, aes(x = y/total)) +
      geom_histogram(aes(fill = totalBin))
   scale.factor <- layer_scales(p)$y$get_limits()[2]/max(pdfvals$density)
    pdfvals$scaled.density <- pdfvals$density*scale.factor
    p <- p +
      geom_line(data = pdfvals, aes(x=x, y=scaled.density)) +
      scale_y_continuous(sec.axis = sec_axis(~ . /scale.factor, name='fitted density')) +
      scale_fill_manual(name = 'total UMI/cell', values = colors.subset) +
      theme_classic() +
      theme(axis.title = element_text(size=10), axis.text = element_text(size=9),
            plot.title = element_text(size=12)) +
      ggtitle(genes[i])
    if (is.null(save)) {
      print(p)
    } else {
      savename <- paste0(save,'-',genes[i],'.png')
      ggsave(file=savename, plot=p, height=2, width=4.5)
      message(paste('saved', savename))
    }
  }
  # next plot the MLE and confidence intervals for all genes
  pci <- ggplot(df.pdf[which(df.pdf$gene %in% genes),], aes(x = as.factor(gene))) +
    geom_point(aes(y = p)) +
    geom_errorbar(aes(ymin = p.lower, ymax=p.upper)) +
    geom_hline(yintercept=0.5, lty='dashed', color = 'grey') +
    ylim(c(0,1)) +
    theme_classic() +
    theme(axis.title = element_text(size=10), axis.text = element_text(size=9),
          plot.title = element_text(size=12),
          axis.text.x = element_text(size=9, angle=45, hjust=1)) +
    xlab('gene') +
    ylab('estimated p')
  if (is.null(save)) {
    print(pci)
  } else {
    savename <- paste0(save,'-confint.png')
    ggsave(file=savename, plot=pci, height=3, width=2)
    message(paste('saved', savename))
  }
  # finally, plot p and var(p) for all genes
  mybreaks <- quantile(df.pdf$totalUMI/df.pdf$totalCells,probs=seq(0,1,0.2))
  df.pdf$totalBin <- cut(df.pdf$totalUMI/df.pdf$totalCells,
                         breaks=mybreaks,
                         include.lowest=T)
  mycolors <- plasma(length(unique(df.pdf$totalBin)))
  pidx <- which(df.pdf$gene %in% genes)
  labelthese <- df.pdf[pidx,]
  pphi <- ggplot(df.pdf, aes(x = p, y = var.p)) +
    geom_point(aes(color = totalBin), alpha=1,size=0.5) +
    scale_color_manual(name='total UMI/gene/cell', values=mycolors) +
    annotate('point', df.pdf$p[pidx], df.pdf$var.p[pidx], color='black') +
    geom_label_repel(data = labelthese, aes(label=gene), fill='white',
               size=3, label.padding=0.05, min.segment.length=0) +
    theme_classic() +
    theme(axis.title = element_text(size=10), axis.text = element_text(size=9),
          plot.title = element_text(size=12)) +
    xlab(TeX(r'($\hat{p}$)')) +
    ylab('variance')
  pphi <- ggMarginal(pphi, type = 'histogram')
  if (is.null(save)) {
    print(pphi)
  } else {
    savename <- paste0(save,'-all.png')
    ggsave(file=savename, plot=pphi, height=3, width=5)
    message(paste('saved', savename))
  }
}
