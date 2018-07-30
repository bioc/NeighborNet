#' Neighbor Net: An approach to infer putative disease-specific mechanisms using neighboring gene networks.
#'
#'
#' @param de a vector including the differentially expressed genes; \code{de} must use the same id's as \code{ref} and the genes in the \code{listofgenes}
#' @param ref the reference vector for all genes in the analysis
#' @param listofgenes a list representing the neighbor networks associated to each gene; the name of the list must be the same as genes in the \code{de}
#' @param threshold a threshold of choosing significant neighbor networks (default is 0.1)
#' @param minsize minimum size of the neighbor networks that should be considered in the analysis (default is 2)
#'
#' @details
#'
#' See details in the cited articles.
#'
#' @return
#'
#' An object of class \code{\link{graphNEL}}.
#'
#' @author
#'
#' Sahar Ansari and Sorin Draghici
#'
#' @references
#'
#' Sahar Ansari, Michele Donato, Nafiseh Saberian, Sorin Draghici; An approach to infer putative disease-specific mechanisms using neighboring gene networks, Bioinformatics, Volume 33, Issue 13, 1 July 2017, Pages 1987–1994
#'
#' @examples
#'
#' # load multiple colorectal cancer study (public data available in GEO
#' # ID: GSE4183, GSE9348, GSE21510, GSE32323, GSEl8671)
#' # These files contains the tables, produced by the limma package with
#' # added gene information.
#' # The table contains the expression fold change and signficance of each
#' # probe set comparing colorectal cancer disease and normal.
#' load(system.file("extdata/dataColorectal4183.RData", package = "NeighborNet"))
#' load(system.file("extdata/dataColorectal9348.RData", package = "NeighborNet"))
#' load(system.file("extdata/dataColorectal21510.RData", package = "NeighborNet"))
#' load(system.file("extdata/dataColorectal32323.RData", package = "NeighborNet"))
#' load(system.file("extdata/dataColorectal8671.RData", package = "NeighborNet"))
#' head(dataColorectal4183)
#'
#' load(system.file("extdata/listofgenes.RData", package = "NeighborNet"))
#' head(listofgenes)
#'
#' # select differentially expressed genes for each data set at p-value below 1%
#' # and absolute value for more than 1.5 and save their entrez ID in a vector de1 to de5
#' pvThreshold <- 0.01
#' foldThreshold <- 1.5
#' de1 <- dataColorectal4183$EntrezID [
#'   dataColorectal4183$adj.P.Val < pvThreshold &
#'   abs(dataColorectal4183$logFC) > foldThreshold]
#' de2 <- dataColorectal9348$EntrezID [
#'   dataColorectal9348$adj.P.Val < pvThreshold &
#'   abs(dataColorectal9348$logFC) > foldThreshold]
#' de3 <- dataColorectal21510$EntrezID [
#'   dataColorectal21510$adj.P.Val < pvThreshold &
#'   abs(dataColorectal21510$logFC) > foldThreshold]
#' de4 <- dataColorectal32323$EntrezID [
#'   dataColorectal32323$adj.P.Val < pvThreshold &
#'   abs(dataColorectal32323$logFC) > foldThreshold]
#' de5 <- dataColorectal8671$EntrezID [
#'   dataColorectal8671$adj.P.Val < pvThreshold &
#'   abs(dataColorectal8671$logFC) > foldThreshold]
#' all <- unique( c(dataColorectal4183$EntrezID, dataColorectal9348$EntrezID,
#'   dataColorectal21510$EntrezID, dataColorectal32323$EntrezID,
#'   dataColorectal8671$EntrezID))
#' de <- unique( c(de1,de2,de3,de4,de5))
#'
#' sig_net <- neighborNet (de, all, listofgenes)
#'
#' @import graph
#' @importFrom stats fisher.test p.adjust
#' @importFrom methods as
#' @export
neighborNet <- function(de,ref,listofgenes,threshold=0.1,minsize=2){

    ## remove any duplication and null values from the de and ref vector
    DE <- de
    ALL <- ref
    DE <- unique(DE)
    ALL <- unique(ALL)
    DE <- DE [!is.na(DE)]
    ALL <- ALL [!is.na(ALL)]

    ## remove neighbor networks with size lower than minsize (default=2)
    listofgenesComplete <- listofgenes
    len <- unlist(lapply(listofgenes,function(x) length(x)))
    listofgenes <- listofgenes [  (len > minsize)]

    ## calculate a hypergoemetric p-value for each neighbor network
    DE_inNetworks <- unlist(lapply(listofgenes,function(x) sum(x %in% DE)))
    genes_in_background <- unique( c( unlist(listofgenes),names(listofgenes) ))
    totalDE <- sum( DE %in% genes_in_background)
    total_inTest <- sum (genes_in_background %in% ALL)
    commonGenes_inNetworks <- unlist(lapply(listofgenes,function(x) sum(x %in% ALL)))
    pv <- unlist(lapply(c(1:length(listofgenes)), function(x){
      testor <- rbind(c(DE_inNetworks[[x]],(commonGenes_inNetworks[x]-DE_inNetworks[[x]])),
                      c(totalDE-DE_inNetworks[[x]],(total_inTest - commonGenes_inNetworks[x]-totalDE+DE_inNetworks[[x]])))
      fisher.test(testor,alternative = "greater")$p.value
    }))
    names(pv) <- names(listofgenes)

    ### correct for multiple comparison
    pvFDR <- p.adjust(unlist(pv),"fdr")
    ### select significant neighbor networks and the genes included
    sigNets  <- pvFDR [ pvFDR < threshold ]
    genes_in_sigNets <- unique(c( unlist(listofgenes[names(sigNets)]), names(sigNets) ))


    ####Connect identified nodes
    if ( length(sigNets) > 0 ){
      adj_sigGenes <- matrix(0,nrow=length(genes_in_sigNets),ncol=length(genes_in_sigNets))
      rownames(adj_sigGenes) <- colnames(adj_sigGenes) <- genes_in_sigNets
      for(i in 1:length( names(sigNets) )){
        adj_sigGenes[ names(sigNets)[i], listofgenes[[ names(sigNets)[i] ]] ] <- 1
      }

      ######### make the graphs  undirected
      for( i in 1:nrow(adj_sigGenes)){
        for(j in 1:i){
          if ( adj_sigGenes[i,j] == 1 | adj_sigGenes[j,i] == 1){adj_sigGenes[i,j] <- adj_sigGenes[j,i] <- 1 }
        }
      }

      ######### remove self edges
      for(i in 1:dim(adj_sigGenes)[1] ){
        adj_sigGenes[ i, i ] <- 0
      }

      ########## construct the final network
      # library(graph)
        ghList <- as(adj_sigGenes,"graphNEL")
    }

    ### return NULL if there is no significant neighbor networks
    if (length(sigNets) == 0 ) { ghList <- NULL }

    ghList

}
