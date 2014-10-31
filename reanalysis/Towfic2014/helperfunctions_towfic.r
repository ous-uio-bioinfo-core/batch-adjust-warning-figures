

# Towfic et al.  For the correct GEO accession GSE61901
# load data and probe annotation and some formatting
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
	duplicatesampleannotation = pData(eset)[, c("characteristics_ch1", "characteristics_ch1.3", "characteristics_ch1.2")]
	names(duplicatesampleannotation) = c("batch", "covariate", "array_strip_address")
	duplicatesampleannotation[,1]=gsub("batch: ", "", duplicatesampleannotation[,1] )
	duplicatesampleannotation[,2]=gsub("batchcovar: ", "", duplicatesampleannotation[,2] )
	duplicatesampleannotation[,3]=gsub("array_address: ", "", duplicatesampleannotation[,3] )
	duplicatesampleannotation[,"array_hyb_address"]=gsub("_[12]", "", duplicatesampleannotation[,3] )
	
	# some renaming for conveniance
	duplicatesampleannotation$covariate = make.names(duplicatesampleannotation$covariate)
	duplicatesampleannotation$covariate[duplicatesampleannotation$covariate=="GA.DP"] = "DP"
	duplicatesampleannotation$covariate[duplicatesampleannotation$covariate=="GA.Q"] = "N"
	duplicatesampleannotation$covariate[duplicatesampleannotation$covariate=="Medium"] = "M"
	duplicatesampleannotation$covariate[duplicatesampleannotation$covariate=="GA.RS"] = "RS"
	rownames(duplicatesampleannotation) = duplicatesampleannotation[, "array_strip_address"]
	
	# Creating a sample annotation for the data set where replicates are combined.
	sampleannotation=duplicatesampleannotation[(1:(nrow(duplicatesampleannotation)/2))*2,c(1,2,4)]
	rownames(sampleannotation) = sampleannotation[, "array_hyb_address"]
	
	# Currently takes about ? minutes to read in.
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
	#table(colnames(datamatrix_raw)==paste( duplicatesampleannotation[, "array_strip_address"], sep="")) 
	colnames(datamatrix_raw)=duplicatesampleannotation[, "array_strip_address"]
	
	# probeid check. Different order from the probeannotation, but all there.
	# probeannotation is not used further.
	#gplannotation =  getGEO(annotation(eset))
	#probeannotation = data.frame(Table(dataTable(gplannotation)), stringsAsFactors=FALSE)
	#probeannotation[] <- lapply(probeannotation, as.character)
	#table(rawdataprobeids==probeannotation$ID)
	#table(rawdataprobeids %in% probeannotation$ID)
	#table(probeannotation$ID %in% rawdataprobeids)
	
	rownames(datamatrix_raw) = rawdataprobeids
	
	return(list(duplicatesampleannotation=duplicatesampleannotation, data=datamatrix_raw, probes=probeannotation))
	
	rm(orgrawtable)
}

# Towfic et al.  For the first geo accession, GSE40566
# load data and probe annotation and some formatting
loadtowfic = old_function(downloaddata=TRUE)
{
  
  sampleannotation = read.table("data/sampleannotation.csv", sep="\t",
                                header=TRUE,  stringsAsFactors=FALSE)
  sampleannotation$code = make.names(sampleannotation$code)
  sampleannotation$chip = as.character(sampleannotation$chip)
  dimnames(sampleannotation)[[1]] = sampleannotation$code
  # take out 3 samples that are not assign to a geoaccession. Failed QC?
  sampleannotation = sampleannotation[!is.na(sampleannotation$geoaccession),] 
  
  
  if(downloaddata)
  {
    temp = tempfile()
    download.file(url="http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40566&format=file&file=GSE40566%5Fnon%5Fnormalized%2Etxt%2Egz",
                  destfile=temp, mode = "wb")
    rawdata = read.table(temp, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    unlink(temp)
    
    temp = tempfile()
    download.file(url="http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40566&format=file",
                  destfile=temp, mode = "wb")
    tardirtemp = file.path(tempdir(), "annottmp")
  	dir.create(tardirtemp)
    untar(temp, exdir = tardirtemp,tar="internal")
    rawannotation = read.table(paste(tardirtemp ,"/GPL6887_MouseWG-6_V2_0_R3_11278593_A.txt.gz", sep=""), 
                               sep="\t", header=TRUE, stringsAsFactors=FALSE, 
                               skip=8, comment.char="", quote="", fill=TRUE)
    
    unlink(temp)
    unlink(tardirtemp, recursive=TRUE)
  }else{
    
    # if download did not work change downloaddata to FALSE  
    # download data from the GEO deposit
    # unpack and place files in a folder named "not_in_github"
    rawdata = read.table("not_in_github/GSE40566_non_normalized.txt", 
                         sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    # the probe annotation file found inside GSE40566_RAW.tar
    rawannotation = read.table("not_in_github/GPL6887_MouseWG-6_V2_0_R3_11278593_A.txt", 
                               sep="\t", header=TRUE, stringsAsFactors=FALSE, 
                               skip=8, comment.char="", quote="", fill=TRUE)
  }  
  
  tmp =  rawannotation$Species=="Mus musculus"
  experimentalannot = rawannotation[tmp,]
  experimentalannot$Array_Address_Id = as.numeric(experimentalannot$Array_Address_Id)
  controlannot = rawannotation[!tmp,]
  dimnames(controlannot)[[2]] = rawannotation[rawannotation[,2]=="Array_Address_Id",]
  controlannot$Array_Address_Id = suppressWarnings(as.numeric(controlannot$Array_Address_Id))
  controlannot = controlannot[!is.na(controlannot$Array_Address_Id),]
  controlannot=controlannot[,1:6]
  probeannotation = merge(experimentalannot, controlannot, all=TRUE )
  #dim(probeannotation)
  rm(tmp, experimentalannot, controlannot)
  
  probeannotation = probeannotation[!duplicated(probeannotation$Array_Address_Id),]
  probeannotation = probeannotation[probeannotation$Array_Address_Id %in% rawdata$ID_REF, ] # 
  probeannotation$Symbol=tolower(probeannotation$Symbol)
  dimnames(probeannotation)[[1]] = probeannotation$Probe_Id
  #dim(probeannotation)
  
  #sort and filter probe and data similar.
  datamatrix_raw = as.matrix(rawdata[,-1])
  datamatrix_raw = datamatrix_raw[match( probeannotation$Array_Address_Id , rawdata$ID_REF), ]
  dimnames(datamatrix_raw)[[1]] = probeannotation$Probe_Id
  #dim(datamatrix_raw)
  #dim(probeannotation)
  
  #and match data to samples.
  #table(sampleannotation$code %in% dimnames(datamatrix_raw)[[2]])# check
  #table(dimnames(datamatrix_raw)[[2]] %in% sampleannotation$code)# check
  datamatrix_raw = datamatrix_raw[, match(sampleannotation$code , dimnames(datamatrix_raw)[[2]])]
  
  return(list(sampleannotation=sampleannotation, data=datamatrix_raw, probes=probeannotation))
}