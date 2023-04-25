#Build vector 'host' indicating in which host each node of the ctree is found 
#Added ttree (multitype) as input so as to include host types
.computeMultitypeHost = function(ctree, ttree)  {
  fathers <- rep(NA, nrow(ctree)) 
  fathers[ctree[ ,2] + 1] <- 1:nrow(ctree) 
  fathers[ctree[ ,3] + 1] <- 1:nrow(ctree) 
  fathers <- fathers[-1] 
  host <- rep(0, nrow(ctree)) 
  nsam <- sum(ctree[ ,2] == 0&ctree[ ,3] == 0) 
  for (i in 1:nsam) { 
    j <- i 
    while (1)  { 
      if (host[j]>0) warning('Warning: two leaves in same host')
      host[j] <- i 
      j <- fathers[j] 
      if (ctree[j,3] == 0) break 
    } 
  } 
  if (nsam==1) return(host)

  dispo=nsam+1
  f <- which( ctree[,3] == 0 & ctree[,2]>0 & (1:nrow(ctree))<nrow(ctree) ) #transmission events other than root
  for (i in 1:length(f)) { 
    j <- f[i] 
    tocol <- c() 
    while (ctree[fathers[j],3]>0&&fathers[j]<nrow(ctree))  { 
      tocol <- c(tocol,j) 
      j <- fathers[j] 
    } 
    if (host[j]==0) {host[j]=dispo;dispo=dispo+1}
    host[tocol] <- host[j] 
  } 
  # now that hosts have been computed, add types on from the ttree
  types <- rep(0, length(host))
  # type for individual 0 (root node infector) is "0" as a placeholder
  for (i in 1:length(host[-1])){
    types[i] <- ttree[host[i],4]
  }
  return( unname(cbind(host, types)) )
} 
