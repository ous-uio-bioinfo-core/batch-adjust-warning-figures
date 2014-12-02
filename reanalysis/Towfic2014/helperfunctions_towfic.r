
# Downloads the non-normalized data from GSE40566
# Using the covariate annotation from GSE61901
loadtowfic = function(downloaddata=TRUE)
{
	geoaccession="GSE40566"
  rawfn = "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40566&format=file&file=GSE40566%5Fnon%5Fnormalized%2Etxt%2Egz"	
	
	#createsampleannotation() # ad.hoc function to create sampleannotation.txt based on GSE61901 and GSE40566
	sampleannotation = read.table("sampleannotation.txt", sep="\t",
                                header=TRUE, stringsAsFactors=FALSE)
	
	if(downloaddata)
	{
		temp = tempfile()    
		download.file(url=rawfn, destfile=temp, mode = "wb")
		datamatrix_raw = read.table(temp, sep="\t", header=TRUE, 
														 stringsAsFactors=FALSE, skip=0, strip.white=TRUE, fill=TRUE)
		unlink(temp)
	}else{
		datamatrix_raw = read.table(paste("not_in_github/",geoaccession,"_non-normalized.txt", sep=""), sep="\t", header=TRUE, 
														 stringsAsFactors=FALSE, skip=0, strip.white=TRUE, fill=TRUE)
	}
	
	rownames(datamatrix_raw) =  datamatrix_raw[,1]
	datamatrix_raw = datamatrix_raw[,-1]
	datamatrix_raw = as.matrix(datamatrix_raw)

	sampleannotation = sampleannotation[ match(sampleannotation$title, colnames(datamatrix_raw) ) , ]
	return(list(sampleannotation=sampleannotation, data=datamatrix_raw))
}	

# ad hoc function to fetch the covariate labels which are described in GSE61901 but not in GSE40566
# GSE40566 has a "relation" tag that connects to the sample in GSE61901 which has covaraite label.
# I dont know a method of getting the sample annotation from GEO whitout downloading the whole matrix (including the GPL).
createsampleannotation = function()
{
	require("GEOquery")
	GSE61901eset=getGEO("GSE61901")[[1]]
	GSE40566eset=getGEO("GSE40566")[[1]]
	
	#t(pData(GSE61901eset)[1,])
	GSE61901sa = pData(GSE61901eset)[, c("characteristics_ch1", "characteristics_ch1.3", "characteristics_ch1.2")]
	names(GSE61901sa) = c("batch", "covariate", "array_strip_address")
	GSE61901sa[,1]=gsub("batch: ", "", GSE61901sa[,1] )
	GSE61901sa[,2]=gsub("batchcovar: ", "", GSE61901sa[,2] )
	GSE61901sa[,3]=gsub("array_address: ", "", GSE61901sa[,3] )
	GSE61901sa[,"array_hyb_address"]=gsub("_[12]", "", GSE61901sa[,3] )
	# some renaming for conveniance
	
	GSE61901sa$covariate = make.names(GSE61901sa$covariate)
	GSE61901sa$covariate[GSE61901sa$covariate=="GA.DP"] = "DP"
	GSE61901sa$covariate[GSE61901sa$covariate=="GA.Q"] = "N"
	GSE61901sa$covariate[GSE61901sa$covariate=="Medium"] = "M"
	GSE61901sa$covariate[GSE61901sa$covariate=="GA.RS"] = "RS"
	
	#t(pData(GSE40566eset)[1,])
	GSE40566sa  = pData(GSE40566eset)[, c("title", "relation")]
	GSE40566sa$title = make.names(GSE40566sa$title)
	GSE40566sa[,"relation"]=gsub("Reanalyzed by: ", "", GSE40566sa[,"relation"] )
	GSE40566sa = cbind(GSE40566sa, GSE61901sa[GSE40566sa$relation,])
	write.table(GSE40566sa, file="sampleannotation.txt", sep="\t", col.names=NA, quote=FALSE)
}





