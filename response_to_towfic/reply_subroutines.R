###
### Reply subroutines
###

# Write to file
write.batch.sample.distribution <- function() {
  write.csv(as.matrix(table(unique_samples$labels.BATCH,unique_samples$labels.CLASS)),file='/workspace/copaxone/stage2/rna/ur_02112016_nygaard_data_balance/sample_table.csv',quote=F)
}

# Permute samples (replicates in positions 2j-1 & 2j for j=1,...)
#    rows = TRUE for permutation per row, FALSE for all rows permutet together
permute_dp_q_samples <- function(data,labels,rows=FALSE,seed=NULL) {
  permdata = data
  if (!is.null(seed)) {
    set.seed(seed)
    attr(permdata,'log')$seed=seed
  }
  used=labels$CLASS %in% c('GA DP','GA Q')
  for(b in unique(labels$BATCH)) {  
    sel = labels$BATCH==b & used
    if(any(sel)) {
      if (rows) permdata[, sel] = t(apply(data[, sel], 1 , permute_dp_q_row))  # randomise each row separately
      else permdata[, sel] = permute_dp_q_all(data[, sel])  # randomise all rows simultaneously
    }
  }
  attr(permdata,'log')$randomise=TRUE
  attr(permdata,'log')$randomise.rows=rows
  tag=paste0(if (rows) 'RndPerRow' else 'RndSmp',seed)
  attr(permdata,'log')$tags <- c(attr(permdata,'log')$tags,tag)
  return(permdata)
}

permute_dp_q_row_alt = function(rowdata) {
  u=sample(length(rowdata)/2)
  ret = rowdata
  ret[c(TRUE,FALSE)] = rowdata[2*u-1]
  ret[c(FALSE,TRUE)] = rowdata[2*u]
  return(ret)
}
permute_dp_q_row = function(rowdata)
{
  a=sample((1:length(rowdata))[c(TRUE,FALSE)])
  b=a+1
  ret = rowdata
  ret[c(TRUE,FALSE)] = rowdata[a]
  ret[c(FALSE,TRUE)] = rowdata[b]
  return(ret)
}

permute_dp_q_all = function(rowdata) {
  u=sample(length(rowdata)/2)
  ret = rowdata
  ret[,c(TRUE,FALSE)] = rowdata[,2*u-1]
  ret[,c(FALSE,TRUE)] = rowdata[,2*u]
  return(ret)
}

# BH-correction to get Q-values
N.sign = function(P,method='BH',limit=0.05) if (sum(!is.na(P))) length(which(p.adjust(P,method=method) < limit)) else NA

# Return *-string to represent how small a number is
simstar = function(x,f=10,max=5,symbol='*') {
  if (x>1 || x<0) {
    return('')
  } else if (x<f**(-max)) {
    n=max
  } else {
    n=ceiling(-log(x)/log(f))
  }
  return(cat(rep(symbol,n),sep=''))
}

# Histogram for pairs plot
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
