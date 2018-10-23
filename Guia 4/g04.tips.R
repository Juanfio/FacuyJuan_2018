


#matriz de correlacion entre filas de una matriz A
cor(t(A))

#defino matriz de similaridad  0 < s_ij < 1
sim <- (1+cor(t(A)))/2
 
# Como obtener los elementos de la matriz triangular superior
a <- A[upper.tri(A)]
 
# Como graficar la densidad de probabilidad de una ristra de valores 
plot(density(a))


# Binarizar una matriz de numeros reales de acuerdo a si superan un valor umbral o no
 A[A<umbral]<- 0
 A[A>umbral]<- 1
 
 ## Armado de grafo y particiones
 # armo el grafo a partir de la matriz de adyacencia
 gc  <- graph_from_adjacency_matrix(A,mode="undirected",diag=FALSE)
 
 # Deteccion de comunidades por fastgreedy
 fgc <- fastgreedy.community(gc)
 
 # Deteccion de comunidades por infomap
 ifm <- infomap.community(gc)
 
 # calculo modularidad de la particion fgc
 mod  <- modularity(gc,membership(fgc)) 
 
 ## Visualizacion de grafo
 # parametros de visualizacion
 V(gc)$size        <- 4
 V(gc)$frame.color <- "white"
 E(gc)$curved      <- 0.2
 lo <- layout_nicely(gc)
 
 # Para fijar el color de cada nodo coloreo nodos segun pertenencia a clusters dandole color a los 10 clusters mas grandes
 ccolor <- c(rainbow(10),"#CCCCCCCC")[pmin(11,membership(fgc))]
 V(gc)$color       <- ccolor
  
 saux <- paste("fastgreedy modularity:",signif(mod,2))
 plot(gc, vertex.label=NA,layout=lo,main=saux)
 
 
 
 ## Clustering Jerarquico 
 #Matriz de distancia
 d<-1-(1+cor(t(geneX)))/2
 
 #Clustering Jerarquico
 method1 <- c("single","complete","average")[2]
 h<-hclust(d=as.dist(d),method=method1)
 
 #ploteo dendrograma
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
   matplot(hpi ,t(z[ia,1:12]),typ="b",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),ylim=range(t(z[ia,])))
   mtext(paste("preCluster",i," - ",aux),cex=0.7)

   #calculo y agrego perfil medio del cluster
   xx<-apply(z[ia,],2,mean)
 
   lines(hpi,xx[1:12],lwd=7,col=rainbow(10,alpha=0.4)[ilabel])
   
 }

