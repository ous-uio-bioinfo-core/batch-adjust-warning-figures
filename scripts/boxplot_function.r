
adhoc.cex=1
adhoc.legendcex=adhoc.cex
adhoc.pch = c(1, 2, 3, 4)


adhoc.usecolor = function(usecolor=TRUE) {
  if (usecolor) {
    adhoc.palette <<- c("darkred", "darkblue", "darkgreen")
    adhoc.palette.line <<- c("red", "blue", "seagreen"); # adjustcolor(adhoc.palette[g], linealpha)
    adhoc.palette.fill <<- c("lightpink", "lightblue", "seagreen1"); # adjustcolor(adhoc.palette[g],boxalpha))
    adhoc.batch.colour<<- "darkgray"; #adjustcolor("black", 0.3)
  } else {
    adhoc.palette <<-  c("black","black","black","black")
    adhoc.palette.line <<- c("darkgray", "darkgray", "darkgray"); # adjustcolor(adhoc.palette[g], linealpha)
    adhoc.palette.fill <<- c("lightgray", "lightgray", "lightgray"); # adjustcolor(adhoc.palette[g],boxalpha))
    adhoc.batch.colour <<- "darkgray"; #adjustcolor("black", 0.3)
  }
}
adhoc.usecolor();

plot_one_gene = function(y, group, batch=NULL, ylim=NULL, main="", estimatemethod="none", lwd=1, boxlabel='CI',leftmargin=7, bbh=NA, bblo=NA)
{  
	
	
  # Boxplots for CI etc.  
  xboxplots=round(length(y)* 1.1)    
  boxwidth=round(length(y)/20)
  boxseparation = (boxwidth * 1.5)
 
  
  a = order(group)
  if(!is.null(batch))
  {
    a = order(batch, group)
    batch=batch[a]
  }
  group=group[a]
  
  y=y[a]  
  
  if(is.null(ylim))
  {
    ymax = round(max(y)+2, 1)
    ymin = round(min(y)-2, 1)
    ylim=c(ymin,ymax)
  }
  #xlim = c(1, round(xboxplots + (boxseparation) * length(unique(group)))  )
	xlim = c(1-leftmargin, length(y))
  
  # measurements
  plot(y, ylim=ylim, xlim=xlim, col=adhoc.palette[as.factor(group)],  pch=adhoc.pch[as.factor(group)], main=main,  ylab=NA, xlab=NA, lwd=lwd, cex.main=adhoc.legendcex, xaxt="n", yaxt="n")
  
  # batch separator
  if(!is.null(batch))    
  {
    for(b in unique(batch))
    {
      batchnames=unique(batch)
      #x = match(batchnames, batch)
      #abline(v = x[-1]-0.5, lty=3)#
    
      xleft = match(b, batch)
      xright = length(batch) - match(b, batch[length(batch):1]) + 1
      batchmean = mean(y[batch==b])
      ybottom = batchmean - bblo[as.numeric(b)]
      ytop = ybottom + bbh    
      rect(xleft, ybottom, xright, ytop, lty=3, border=adhoc.batch.colour)
    
    # Decide where to put batch label, up of down. furthest from spot 1.
      if( (y[batch==b][1]-ybottom)  >  (ytop - y[batch==b][1]) ) # bottom
      {
        labpos=3
        laby=ybottom
      }else{ # top
        labpos=1
        laby=ytop
      }
      text(labels=paste("             Batch", b, sep=""), y=laby, x=xleft, pos = labpos , offset=0.3,
        cex=adhoc.legendcex, col=adhoc.batch.colour)
    }
  }
  
  
  groupnames=factor(unique(group))
  for(g in groupnames)
  {
    g=factor(g,levels=levels(groupnames))
    
    if(estimatemethod=="lsmeans")
    {
      fit_lm=lm(y ~ group+batch)
      means_lm=lsmeans(fit_lm,~group)
      anovaest = summary(means_lm)      
      m =anovaest[g , "lsmean"]      
      ybottom = anovaest[g , "lower.CL"]      
      ytop = anovaest[g , "upper.CL"]      
    } else if(estimatemethod=="CI") {
      m = mean(y[group==g])
      ci = t.test(y[group==g])[["conf.int"]]
      ybottom = ci[1]
      ytop = ci[2]
    }
    
    
    xextra = 0.3
    xleft = match(g, group) - xextra
    xright = length(group) - match(g, group[length(group):1]) + 1 + xextra
    
    linealpha=0.2
    boxalpha=0.1
    rect(xleft, ybottom, xright, ytop, border=adhoc.palette.line[g], lwd=lwd, 
         lty="solid", density=-1, col=adhoc.palette.fill[g])
    lines(c(xleft,xright), c(m,m), col=adhoc.palette.line[g], lwd=lwd)
    
    if(g==1)
      text(labels=boxlabel, x=xleft, y=m, pos=2, offset = 0.3)
    

  # measurements: display
  points(y, col=adhoc.palette[as.factor(group)],  pch=adhoc.pch[as.factor(group)], lwd=lwd)
  } 
  
}

estimatesboxesonly = function(y, group, batch, ylim=NULL, main="", lwd=1)
{
  
  groupnames=factor(unique(group))
  plot( as.numeric(unique(group)),  xlim=c(1, length(unique(group))+1) , ylim=ylim, type="n", xaxt="n", yaxt="n", main=main, cex.main=adhoc.legendcex, xlab=NA, ylab=NA)
  for(g in groupnames)
  {
    g=factor(g,levels=levels(groupnames))
    
    
      fit_lm=lm(y ~ group+batch)
      means_lm=lsmeans(fit_lm,~group)
      anovaest = summary(means_lm)      
      m =anovaest[g , "lsmean"]      
      ybottom = anovaest[g , "lower.CL"]      
      ytop = anovaest[g , "upper.CL"]      
    

    xleft = as.numeric(g)
    xright = as.numeric(g) + 0.9
    
    linealpha=0.2
    boxalpha=0.1
    rect(xleft, ybottom, xright, ytop, border=adhoc.palette.line[g], lwd=lwd, 
         lty="solid", density=-1, col=adhoc.palette.fill[g])
    lines(c(xleft,xright), c(m,m), col=adhoc.palette.line[g], lwd=lwd)
    #points(x=(xleft+xright)/2, y=ytop+( (ylim[2]-ylim[1])/30), pch=adhoc.pch[g], cex=adhoc.legendcex, lwd=lwd, col=adhoc.palette[g])

  }
}
