#' Plot a transmission tree in an economic format
#' @param ttree Transmission tree
#' @param showLabels Boolean for whether or not to show the labels
#' @param hostTypeNames Vector of names for host types, ordered.
#' @param showMissingLinks Option for how to show missing links: (0) as dots, (1) as several gray levels, (2) as a single gray level.
#' @param maxTime Maximum value of time to show on x axis
#' @param cex Expansion factor
#' @return Returns invisibly the first parameter
#' @examples 
#' plotTTree2(extractTTree(simulateOutbreak()))
#' @export
plotTTree2 = function(ttree,showLabels=TRUE, hostTypeNames = NA,showMissingLinks=0,maxTime=NA,cex=1) {
  nam=ttree$nam
  ttree=ttree$ttree
  multitype <- ifelse( ncol(ttree) == 4, yes = T, no = F )
  if ( multitype ) ttree <- cbind(ttree[,1:3], rep(1,nrow(ttree)), ttree[,4]) 
  else ttree <- cbind(ttree,rep(1,nrow(ttree))) 
  if (showMissingLinks>0) {
    i=which(is.na(ttree[,2]))[1]
    while (i<nrow(ttree)) {
      w=which(ttree[,3]==i)
      if (length(w)==1) {
        ttree[w,3]=ttree[i,3]
        ttree[w,4]=ttree[w,4]+ttree[i,4]
        ttree=ttree[-i,]
        ttree[which(ttree[,3]>i),3]=ttree[which(ttree[,3]>i),3]-1
      } else i=i+1
    }
  }
  
  if (showMissingLinks==2) {
    ttree[which(ttree[,4]>=2),4]=2
  }
  
  #Determine ys 
  n=nrow(ttree)
  ys <- rep(0, n)
  scale <- rep(1,n)
  todo=c(which(ttree[,3]==0))
  while (length(todo)>0) {
    f=which(ttree[,3]==todo[1])
    o=rank(-ttree[f,1])
    f[o]=f
    for (i in f) {ys[i]=ys[todo[1]]+scale[todo[1]]*which(f==i)/(length(f)+1);scale[i]=scale[todo[1]]/(length(f)+1);todo=c(todo,i)}
    todo=todo[-1]
  }
  ys=rank(ys)
  
  
  #Do the plot
  oldpar <- par('yaxt','bty')
  on.exit(par(oldpar))
  par(yaxt='n',bty='n')
  #mi=min(ttree[,2])#,ttree[,1])
  mi=min(ttree[which(!is.na(ttree[,1])),1])
  ma=max(ttree[which(!is.na(ttree[,1])),1])
  if (!is.na(maxTime)) ma=maxTime
  plot(c(),c(),xlim=c(mi,ma),ylim=c(0,n+1),xlab='',ylab='')
  pal <- grDevices::gray.colors(max(ttree[,4]))
  for (i in 1:n) {
    if (ttree[i,3]!=0) {
      dircol <- pal[ttree[i,4]] 
      #arrows(ttree[i,1],ys[ttree[i,3]],ttree[i,1],ys[i],length=0)
      arrows(ttree[ttree[i,3],1],ys[ttree[i,3]],ttree[i,1],ys[i],length=0,col=dircol)
    }
    #ma=max(ttree[i,1],ttree[which(ttree[,3]==i),1])
    #mi=min(ttree[i,1],ttree[i,1],ttree[which(ttree[,3]==i),1])
    # if (showLabels && !is.na(ttree[i,2])) text(ttree[i,1],ys[i],nam[i],pos=4,cex=cex) THIS LINE GIVES LABELS ONLY TO OBSERVED HOSTS
    if (showLabels) text(ttree[i,1],ys[i],i,pos=4,cex=cex) ## I ADDED THIS LINE TO LABEL ALL HOSTS TO CHECK THINGS
    #lines(c(mi,ma),c(ys[i],ys[i]))
  }
  if ( multitype ) type_colors <- grDevices::rainbow(max(ttree[,5]))
  else type_colors <- "black"
  for (i in 1:n) {
    #points(ttree[i,1],ys[i],pch=21,bg='black',cex=0.2+0.2*(!is.na(ttree[i,2])))
    points(ttree[i,1], ys[i],
           pch = ifelse(is.na(ttree[i,2]),21,16),
           bg = "white",
           col = type_colors[ttree[i,5]],
           cex = cex)
  }
  if ( length(pal)>2 ) legend('topleft',legend = 0:(length(pal)-1),col = pal,lty=1,cex=cex,title='Missing links')
  if ( is.na(hostTypeNames) ) hostTypeNames <- paste( "Type", 1:length(type_colors ) )
  if ( multitype ) legend('topleft', legend = hostTypeNames, fill = type_colors, title = "Host Type")
  return(invisible(ttree))
}
plotTTree2(tree)
