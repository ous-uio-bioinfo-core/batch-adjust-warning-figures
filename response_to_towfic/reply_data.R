###
### Reply.R: subsoutines for reading in data
###

# Set DataDir as working directory if not already defined
if (!exists('DataDir')) DataDir='.'

# Download if not already exists
ensure.download <- function(url,destfile,mode='wb') {
  if (!file.exists(destfile)) {
    download.file(url=url,destfile=destfile,mode=mode)
  }  
}

# Data processing logging
Log.print <- function(x) {
  Log <- attr(x,'log')
  for (tag in names(Log)) cat(tag,': ',Log[[tag]],"\n",sep='')
}
Log.name <- function(x,name='data') paste(c(name,attr(x,'log')$tags),collapse='_')

# Read data
ensure.download(url = 'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE61901&format=file&file=GSE61901%5Fnon%2Dnormalized%2Etxt%2Egz',
                destfile = file.path(DataDir,'GSE61901_non-normalized.txt.gz'),mode = 'wb')
if (!exists('data.all')) {
  data.raw <- read.delim(file=file.path(DataDir,'GSE61901_non-normalized.txt.gz'),header=T,skip=5,row.names=1,as.is=T,check.names=F)
  data.raw <- data.raw[,!(colnames(data.raw) %in% c('Detection Pval'))]
  data.all <- normalizeQuantiles(data.raw)
  data.all.avg=apply(data.all,1,mean)
  data.all.sd=apply(data.all,1,sd)
}

# Read labels
getLabels <- function() {
  file <- file.path(DataDir,'GSE61901_series_matrix.txt.gz')
  ensure.download(url='ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE61nnn/GSE61901/matrix/GSE61901_series_matrix.txt.gz',destfile = file,mode = 'wb')
  data <- read.delim(file=file,header=F,skip=32)
  SAMPLE <- as.character(unlist(data[2,-1]));
  CLASS <- as.character(unlist(data[26,-1]));
  SLOT <- as.character(unlist(data[12,-1]));
  SLOT <- gsub(pattern = 'array_address: ',replacement = '',x = SLOT)
  BATCH <- as.character(unlist(data[10,-1]));
  BATCH <- gsub(pattern = 'batch: ','B',BATCH)
  data.frame(SAMPLE=factor(SAMPLE),CLASS=factor(CLASS),SLOT=factor(SLOT),BATCH=factor(BATCH))
}
if (!exists('labels.all')) labels.all <- getLabels()

# Select data
select_data <- function(sample=NULL,avg.min=0,sd.min=0,subset=NULL,subset.tags='subset') {
  if (is.null(subset)) {
    data.subset <- which(data.all.avg>=avg.min & data.all.sd>=sd.min)
    data <- data.all[data.subset,]
    attr(data,'log')$tags <- c(paste0('A',avg.min),paste0('S',sd.min))
    attr(data,'log')$data = paste("Included data with avg>=",avg.min," and sd>=",sd.min)
    attr(data,'log')$avg.min = avg.min
    attr(data,'log')$sd.min = sd.min
  } else {
    data <- data.all[subset,]
    attr(data,'log')$tags <- c(subset.tags)
    attr(data,'log')$data = "Custom data subset"
  }
  if (!is.null(sample) && sample<nrow(data)) {
    cat("Subsample: ",sample)
    data <- data[sample.int(nrow(data),sample),] # Subsample probes
    attr(data,'log')$sample = sample
    attr(data,'log')$tags = c(attr(data,'log')$tags,paste0('n',sample))
  }
  return(data)
}

# Store results
get.tag.file <- function(x,...) paste0(Log.name(x,...),'.rds')
store.results <- function(filename=NULL,file=NULL,dir=DataDir) {
  store <- new.env()
  store$data.dim <- dim(data)
  store$result_summary <- result_summary
  store$P <- P
  attr(store,'log') <- attr(data,'Log')
  store$log <- attr(data,'log')
  attr(store$P,'log') <- attr(data,'log')
  if (is.function(filename)) filename=filename(data)
  if (!is.null(filename)) file=file.path(dir,filename)
  if (!is.null(file)) {
    cat('Saving data to',file,"\n")
    saveRDS(store,file=file)
    cat("Done!\n")
  }
}

