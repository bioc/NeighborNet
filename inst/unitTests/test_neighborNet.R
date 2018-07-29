test_neighbor.net <- function() {

  library("NeighborNet")
  library("graph")

	load(system.file("extdata/dataColorectal4183.RData", package = "NeighborNet"))
	load(system.file("extdata/dataColorectal9348.RData", package = "NeighborNet"))
	load(system.file("extdata/dataColorectal21510.RData", package = "NeighborNet"))
	load(system.file("extdata/dataColorectal32323.RData", package = "NeighborNet"))
	load(system.file("extdata/dataColorectal8671.RData", package = "NeighborNet"))

  load(system.file("extdata/listofgenes.RData", package = "NeighborNet"))

  pvThreshold <- 0.01
  foldThreshold <- 1.5
  de1 <- dataColorectal4183$EntrezID [
    dataColorectal4183$adj.P.Val < pvThreshold &
    abs(dataColorectal4183$logFC) > foldThreshold
  ]
  de2 <- dataColorectal9348$EntrezID [
    dataColorectal9348$adj.P.Val < pvThreshold &
    abs(dataColorectal9348$logFC) > foldThreshold
  ]
  de3 <- dataColorectal21510$EntrezID [
    dataColorectal21510$adj.P.Val < pvThreshold &
    abs(dataColorectal21510$logFC) > foldThreshold
  ]
  de4 <- dataColorectal32323$EntrezID [
    dataColorectal32323$adj.P.Val < pvThreshold &
    abs(dataColorectal32323$logFC) > foldThreshold
  ]
  de5 <- dataColorectal8671$EntrezID [
    dataColorectal8671$adj.P.Val < pvThreshold &
    abs(dataColorectal8671$logFC) > foldThreshold
  ]

  de <- unique( c(de1,de2,de3,de4,de5))

  ref <- unique( c(
    dataColorectal4183$EntrezID,
    dataColorectal9348$EntrezID,
    dataColorectal21510$EntrezID,
    dataColorectal32323$EntrezID,
    dataColorectal8671$EntrezID
  ))

  sig_genes <- neighbor.net(de = de, ref = ref, listofgenes=listofgenes)

  checkEquals(numNodes(sig_genes), 144)
  checkEquals(numEdges(sig_genes), 251)

}
