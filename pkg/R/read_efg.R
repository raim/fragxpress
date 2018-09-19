

read_efg <- function(type=c('aati'), path=getwd(), pattern, 
                     blank, ladder, verb=1) {
  if ( type=='aati')
    fgs <-read_aati(path=path, pattern='.*.Electropherogram.csv$', verb=verb)

  return(fgs)
}


read_aati <- function(path=getwd(), pattern='.*.Electropherogram.csv$', 
                      blank, ladder, verb=1) {
  
  if ( missing(blank) ) blank <- NULL

  fgf <- list.files(path=path,
                    pattern=pattern,
                    recursive=TRUE, include.dirs=TRUE)
  if ( verb>0 )
    cat(paste("found ", length(fgf), "electropherogram file(s)\n"))
  
  fgs <- NULL
  for ( fg in fgf ) {
    
    fgdat <- read.csv(file.path(path, fg), stringsAsFactors=FALSE)
    
    ## parse run time
    time <- strptime(sub(sub("\\.\\*","",pattern),"",basename(fg)), format="%Y %m %d %HH %MM")
        
    if ( verb>0 ) 
      cat(paste0("parsing ", fg, ", ", ncol(fgdat), " columns, ", nrow(fgdat), " rows\n"))
    
    for ( j in 2:ncol(fgdat) ) {
      
      ## parse well & sample info
      wll <- sub("\\.\\..*","",colnames(fgdat)[j])
      smp <- sub(".*\\.\\.","",colnames(fgdat)[j])
      
      ## skip blank samples
      if ( smp%in%blank ) next
    
      
      ## replicate?
      rep <- ifelse(smp%in%fgs$sample,
                    tail(fgs$rep[fgs$sample==smp],1)+1,
                    1)
      ## add to dataframe
      fgs <- rbind(fgs, cbind.data.frame(size=fgdat[,1],
                                         signal=fgdat[,j],
                                         well=wll,
                                         sample=smp, 
                                         rep=rep,
                                         time=time))
    }
  }
  return(fgs)
}

plot_efg <- function(fgs, sample, xlim=c(1,2.5e3), ylim=c(0,1.1),
                     xlab="fragment size, bp", ylab="norm. signal", ...) {
  
  reps <- unique(fgs$rep[fgs$sample==sample])
  dates <- unlist(sapply(reps, function(x) as.character(fgs$time[fgs$rep==x])[1]))
  
  plot(1, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
  
  for ( i in 1:length(reps) ) {
    idx <- fgs$sample==sample & fgs$rep==reps[i]
    midx <- idx & (fgs[,1] > 10 & fgs[,1] < 5000) # normalize only in 1:5000
    lines(fgs[idx,1], fgs[idx,2]/max(fgs[midx,2]),col=i)
    
  }
  legend("topright", legend=dates, col=1:length(reps), lty=1)
}
