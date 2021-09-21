getBowtiePaths <- function(){
	bowtieIndicesPath    <- file.path(.getCrisprIndicesPath(), "bowtie")	
	bowtieIndexHuman     <- file.path(bowtieIndicesPath, "GRCh38/GRCh38")
	bowtieIndexHumanMaj  <- file.path(bowtieIndicesPath, "GRCh38_major/GRCh38_major")
	bowtieIndexHumanMin  <- file.path(bowtieIndicesPath, "GRCh38_minor/GRCh38_minor")
	bowtieIndexMouse     <- file.path(bowtieIndicesPath, "GRCm38/GRCm38")
    bowtieIndexRat       <- file.path(bowtieIndicesPath, "Rnor6/Rnor6")
    bowtieIndexFruitfly  <- file.path(bowtieIndicesPath, "BDGP6/BDGP6")
	bowtieIndexHumanTx   <- file.path(bowtieIndicesPath, "gencode.v34.human/gencode.v34.human")
	bowtieIndexMouseTx   <- file.path(bowtieIndicesPath, "gencode.vM25.mouse/gencode.vM25.mouse")
	bowtieIndexHumanExon <- file.path(bowtieIndicesPath, "igis4.human/exons.igis4.human")
	bowtieIndexMouseExon <- file.path(bowtieIndicesPath, "igis4.mouse/exons.igis4.mouse")
	list(bowtieIndexHuman=bowtieIndexHuman, 
	   bowtieIndexHumanMaj=bowtieIndexHumanMaj,
	   bowtieIndexHumanMin=bowtieIndexHumanMin,
	   bowtieIndexMouse=bowtieIndexMouse,
       bowtieIndexRat=bowtieIndexRat,
       bowtieIndexFruitfly=bowtieIndexFruitfly,
	   bowtieIndexHumanTx=bowtieIndexHumanTx,
	   bowtieIndexMouseTx=bowtieIndexMouseTx,
	   bowtieIndexHumanExon=bowtieIndexHumanExon,
	   bowtieIndexMouseExon=bowtieIndexMouseExon
	)
}
