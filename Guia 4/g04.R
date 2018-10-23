

require(igraph)

#preprocess Data
if(FALSE){
 (load("G04/dataRNAseqVSN.Rdata"))
   
  #Me focalizo en una Temp y encuentro los genes que mas varian
  ttemp   <- 12   #Temperatura del tratamiento
  Nvar    <- 500  #nro de genes mas variable a clusterizar
  
 
  # En los datosd originales cada punto aparece por replicado
  # la siguiente linea promedia replicas
  y      <- t(apply(ycpm.vsn[,sampleNames],1,function(x){return(apply(matrix(x,nrow=2),2,mean))}))
  yhmean <- apply(y,1,function(s){1/mean(1/s)})
  y      <- y[!rownames(y)%in%names(which(yhmean<0)),]
 
  # estos son los nombres de las muestras
  sampleNames<-as.character(pheno[pheno[,"temperature"]==12 & pheno[,"replicate"]=="A","sample"])
 
  #en geneX me quedo solo con los Nvar genes mas variables
  ymad <- apply(y,1,mad)
  iok  <- which(ymad > quantile(ymad,1-Nvar/nrow(y))) 
  geneX <- y[iok,]
  colnames(geneX)<-sampleNames
  
  metadata <- pheno[pheno[,"temperature"]==12 & pheno[,"replicate"]=="A",]
  
  if(FALSE){
   save(geneX,metadata,file="dataG04.Rdata")
  }
}else{
 (load("G04/dataG04.Rdata"))
}
 
 metadata 

##red de coexpresio

 #voy a armar una red con enlaces que correspondan a correlaciones de pearson mayores que simMin
 simMin <- 0.95
 
 #defino matriz de similaridad  0 < s_ij < 1
 sim <- (1+cor(t(geneX)))/2
 sim <- abs(cor(t(geneX)))
 
 # Asi estan distribuidos los valores
 a <- sim[upper.tri(sim)]
 plot(density(a),xlab="correlation similarity")
 abline(v=simMin,col="gray")

 # pongo a cero todos los valores de correlacion < simMin
 # y a 1 los que superen dicho umbral
 ssim             <- sim
 ssim[ssim<simMin]<- 0
 ssim[ssim>0]     <- 1
 
 # armo el grafo a partir de la matriz de adyacencia
 gc  <- graph_from_adjacency_matrix(ssim,mode="undirected",diag=FALSE)
 
 # Deteccion de comunidades por fastgreedy
 fgc <- fastgreedy.community(gc)
 
 # Deteccion de comunidades por infomap
 ifm <- infomap.community(gc)
 
 #Visualizo ambas soluciones
 # parametros de visualizacion
 V(gc)$size        <- 4
 V(gc)$frame.color <- "white"
 E(gc)$curved      <- 0.2
 lo <- layout_nicely(gc)
 
 layout(matrix(1:2,1,2))
 
 # primero fastgreedy
 # coloreo nodos segun pertenencia a clusters dandole color a los 10 clusters mas grandes
 ccolor <- c(rainbow(10),"#CCCCCCCC")[pmin(11,membership(fgc))]
 V(gc)$color       <- ccolor
 
 # calculo modularidad de la particion
 mod  <- modularity(gc,membership(fgc))
 saux <- paste("fastgreedy modularity:",signif(mod,2))
 plot(gc, vertex.label=NA,layout=lo,main=saux)
 
 # lo mismo para infomap
 # coloreo nodos segun pertenencia a clusters
 ccolor <- c(rainbow(10),"#CCCCCCCC")[pmin(11,membership(ifm))]
 V(gc)$color       <- ccolor
 # calculo modularidad de la particion
 mod  <- modularity(gc,membership(ifm))
 saux <- paste("infomap modularity:",signif(mod,2))
 plot(gc, vertex.label=NA,layout=lo,main=saux)
 
 
  
## Clustering Jerarquico 
 
 
 #Matriz de distancia
 d<-1-(1+cor(t(geneX)))/2
 
 #Clustering Jerarquico
 method1 <- c("single","complete","average")[2]
 h<-hclust(d=as.dist(d),method=method1)
 
 layout(matrix(1:2,1,2))
 plot(h)
 
 # Uso modularidad sobre grafo para definir un corte
 res<-c()
 for(numClusters in 2:50){
   ct <- cutree(h,k=numClusters)
   res <- c(res,modularity(gc,ct))
 }
 plot(2:50,res,typ="b",xlab="num clusters",ylab="modularity")
 
 
## Perfil de expresion 
 #elijo una de estas tres particiones
 label <-  membership(fgc)
 label <-  membership(ifm)
 label <-  ct
 
 #grafico perfil de los 6 clusters mas grandes de la particion
 layout(matrix(1:6,2,3,byrow=TRUE))
 par(mar=c(2,2,2,2))

 #para visualizar mejor ESTANDARIZO el perfil de expresion (lo centro y normalizo segun la varianza)
 # z es la matriz ahora de perfiles de expresion estandarizados
 z <- t(apply(geneX,1,function(x){(x-mean(x))/sd(x)}))
 
 #ordeno las clusters por tamanio. 
 #io. tiene los nombres de los clusters segun tamanio decreciente
 io <- order(table(label),decreasing=TRUE)
 

 for(ilabel in io[1:6]){
   ia<-which(label==ilabel)
   aux<-length(ia)

   #ploteo perfiles del cluster
   hpi <- metadata[,"hpi"] 
   #matplot(hpi ,t(z[ia,1:12]),typ="b",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),ylim=range(t(z[ia,])))
   matplot(hpi ,t(z[ia,1:12]),typ="b",ylab="",xlab="",col=rgb(0.2,0.2,0.2),ylim=range(t(z[ia,])))
   mtext(paste("Cluster",ilabel," - ",aux),cex=0.7)

   #calculo y agrego perfil medio del cluster
   xx<-apply(z[ia,],2,mean)
 
   #lines(hpi,xx[1:12],lwd=7,col=rainbow(10,alpha=0.4)[ilabel])
   lines(hpi,xx[1:12],lwd=7,col=rainbow(10,alpha=1)[ilabel])
   
 }
  
if(FALSE){ 
 nncol <- min(6,length(ua))
 npanel<- signif(length(ua)/nncol,0)*nncol
 par(mar=c(2,2,2,2))
  if(bplot1) layout(matrix(1:npanel,ncol=nncol,byrow=TRUE))
  mgenex<-c()
  z <- t(apply(geneX,1,function(x){(x-mean(x))/sd(x)}))
  for(i in 1:min(length(ua),npanel)){
   ia<-which(ct==ua[i])
   aux<-length(ia)
   xx<-apply(z[ia,],2,mean)
 
   mgenex<-rbind(mgenex,apply(z[ia,],2,mean))


   matplot(1:12 ,t(z[ia,1:12]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),xlim=c(1,12),ylim=range(t(z[ia,])))
    
   mtext(paste("preCluster",i," - ",aux),cex=0.7)
   lines(1:12,xx[1:12],lwd=7,col="gray")
   
   if(i%%npanel==0) readline("Sacale una foto y Press <enter>\n")

  }
  
 
 
}

