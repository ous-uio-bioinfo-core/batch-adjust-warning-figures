

# Towfic et al.  For the correct GEO accession GSE61901
# load data and sample annotation and some formatting
loadtowfic = function(downloaddata=TRUE)
{
	geoaccession="GSE61901"
	
	if(downloaddata)
	{
		geoseries=getGEO(geoaccession)
		eset = geoseries[[1]]
	}else{
		eset=getGEO(filename="not_in_github/GSE61901_series_matrix.txt.gz")
	}
	# will only use sample annotation from this eset. The data will be taken from the non-normalized matrix file.
	
	# getting the sample annotation, some renaming
	sampleannotation = pData(eset)[, c("characteristics_ch1", "characteristics_ch1.3", "characteristics_ch1.2")]
	names(sampleannotation) = c("batch", "covariate", "array_strip_address")
	sampleannotation[,1]=gsub("batch: ", "", sampleannotation[,1] )
	sampleannotation[,2]=gsub("batchcovar: ", "", sampleannotation[,2] )
	sampleannotation[,3]=gsub("array_address: ", "", sampleannotation[,3] )
	sampleannotation[,"array_hyb_address"]=gsub("_[12]", "", sampleannotation[,3] )
	
	# some renaming for conveniance
	sampleannotation$covariate = make.names(sampleannotation$covariate)
	sampleannotation$covariate[sampleannotation$covariate=="GA.DP"] = "DP"
	sampleannotation$covariate[sampleannotation$covariate=="GA.Q"] = "N"
	sampleannotation$covariate[sampleannotation$covariate=="Medium"] = "M"
	sampleannotation$covariate[sampleannotation$covariate=="GA.RS"] = "RS"
	rownames(sampleannotation) = sampleannotation[, "array_strip_address"]
		
	
	if(downloaddata)
	{
		temp = tempfile()
		download.file(url="http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE61901&format=file&file=GSE61901%5Fnon%2Dnormalized%2Etxt%2Egz",
									destfile=temp, mode = "wb")
		orgrawtable = read.table(temp, sep="\t", header=TRUE, 
														 stringsAsFactors=FALSE, skip=5, strip.white=TRUE, fill=TRUE)
		unlink(temp)
	}else{
		orgrawtable = read.table("not_in_github/GSE61901_non-normalized.txt", sep="\t", header=TRUE, 
														 stringsAsFactors=FALSE, skip=5, strip.white=TRUE, fill=TRUE)
	}
	
	#to matrix, taking out the empty detection p-val columns.
	datamatrix_raw = orgrawtable
	rawdataprobeids = datamatrix_raw[,1]
	datamatrix_raw = datamatrix_raw[,(1:214)*2]
	datamatrix_raw = as.matrix(datamatrix_raw)
	mode(datamatrix_raw) = "numeric"
	# table(is.na(datamatrix_raw))
	
	#ordercheck
	#table(colnames(datamatrix_raw)==paste( sampleannotation[, "array_strip_address"], sep="")) 
	colnames(datamatrix_raw)=sampleannotation[, "array_strip_address"]
	
	# probeannotation not used
	# probeid check. Different order from the probeannotation, but all there.
	# probeannotation is not used further.
	#gplannotation =  getGEO(annotation(eset))
	#probeannotation = data.frame(Table(dataTable(gplannotation)), stringsAsFactors=FALSE)
	#probeannotation[] <- lapply(probeannotation, as.character)
	#table(rawdataprobeids==probeannotation$ID)
	#table(rawdataprobeids %in% probeannotation$ID)
	#table(probeannotation$ID %in% rawdataprobeids)
	
	rownames(datamatrix_raw) = rawdataprobeids
	rm(orgrawtable)
	return(list(sampleannotation=sampleannotation, data=datamatrix_raw))
}	

	
