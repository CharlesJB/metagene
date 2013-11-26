# Created by Charles Joly Beauparlant
# 2013-11-25

annotateBED <- function(BEDfile) {
	library(ChIPpeakAnno)
	library(org.Mm.eg.db)
	data(TSS.mouse.NCBIM37)
	myPeaks <- BED2RangedData(BEDfile)
	annotatedPeaks <- annotatePeakInBatch(myPeaks, AnnotationData=TSS.mouse.NCBIM37)
	IDs_df <- as.data.frame(addGeneIDs(annotatedPeaks, "org.Mm.eg.db", IDs2Add="ensembltrans"))
	return(as.data.frame(IDs_df))
}
