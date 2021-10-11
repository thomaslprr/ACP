showIndivFun <- function(data,x=data[,1],y=data[,2],label=TRUE){
  library(ggplot2)
  if(label){
    ggplot(data, aes(x=x, y=y)) +
      geom_label(
        label=rownames(data), 
        nudge_x = 0.25, nudge_y = 0.25, 
        check_overlap = T
      )+ geom_hline(yintercept=0,color="black")+ geom_vline(xintercept = 0,color="black")
  }else{
    ggplot(data, aes(x=x, y=y)) +
      geom_point() + # Show dots
      geom_text(
        label=rownames(data), 
        nudge_x = 0.25, nudge_y = 0.25, 
        check_overlap = T
      )+ geom_hline(yintercept=0,color="black")+ geom_vline(xintercept = 0,color="black")
  }
}

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

showCorrelFun <- function(data,xvalue=1,yvalue=2,labelname=rownames(data)){
  
  dat <- circleFun(c(0,0),2,npoints = 700)
  
  #affichage du cercle de corrélation
  library(ggplot2)
  library(ggrepel)
  variable.graphe <- ggplot(data, aes(data[,xvalue], data[,yvalue])) +
    geom_hline(yintercept=0,color="black")+
    geom_vline(xintercept = 0,color="black")+
    geom_point(color = "blue", size = 3)+
    geom_point(aes(x=x, y=y), data=dat,color="black")
  
  variable.graphe + geom_segment(aes(x = 0, y = 0, xend = data[,xvalue], yend = data[,yvalue]),
                                 arrow = arrow(length = unit(0.5, "cm")))+
    geom_label_repel(aes(label = labelname),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'blue') +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank())+
    expand_limits(x=c(-1,1), y=c(-1, 1))+theme(aspect.ratio=1)
}

acpFun <- function(dataVal){
  nblignes <- dim(dataVal)[1]
  nbcols <- dim(dataVal)[2]
 
  #centrée réduite
  X <- scale(dataVal,center = TRUE, scale=TRUE)*sqrt(nblignes/(nblignes-1)); 
  
  #matrice de variance/covariance
  S = t(X)%*%X * (1/nblignes)
 
  #matrice de corrélation
  R <- t(S)%*%S * (1/nblignes) 
  
  var_cum <- round(eigen(S)$values,3)
  #affichage sous forme de graphique
  barplot(main = "Variance cumulée",
          col="blue", 
          var_cum) 
  
  nb_axis <- readline(prompt="Combien d'axe ?")
  # convert character into integer
  nb_axis <- as.integer(nb_axis)

  total <- sum(eigen(S)$values)
  #somme des valeurs des variances cumulées
  round(cumsum(eigen(S)$values*(100/total)),3)  
  
  # 4. ACP Normée diagonalisation de S
  ACP=eigen(S) 
  
  # Inertie Ig
  round(ACP$values,3) 
  
  # axes principaux
  u=ACP$vectors
  
  # composantes principales F : individus
  F=X%*%u 
  
  data <- as.data.frame(F[,1:nb_axis])
  
  
  # G : CP coord des variables
  G<-matrix(0,nbcols,nbcols)
  for(i in 1:nb_axis){
    G[,i]=sqrt(ACP$values[i])%*%ACP$vectors[,i]
  }
  
  #conversion matrice en dataframe
  dataCor=as.data.frame(G[,1:nb_axis])
  
  listOfDataframe = list(
    "datas" = dataVal,
    "indiv" = data,
    "var" = dataCor
  )
  
  return(listOfDataframe)
}
