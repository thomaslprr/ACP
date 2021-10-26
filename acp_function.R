showIndivFun <- function(data,xvalue=1,yvalue=2,label=TRUE,number=TRUE){
  library(ggplot2)
  valueToShow <- c(1:dim(data)[1])
  if(!number){
    valueToShow <- rownames(data)
  }
  base <- ggplot(data, aes(x=data[,xvalue], y=data[,yvalue])) +
    ggtitle("Représentation des individus")+
    theme(plot.title = element_text(hjust = 0.5,face="bold"))+
    labs(y=paste(c("Axe", yvalue), collapse = " "), x =paste(c("Axe", xvalue), collapse = " "))
  
  if(label){
    base +
      geom_label(
        label=valueToShow, 
        nudge_x = 0.25, nudge_y = 0.25, 
        check_overlap = T
      )+ geom_hline(yintercept=0,color="black")+ geom_vline(xintercept = 0,color="black")
  }else{
    base +
      geom_point() + # Show dots
      geom_text(
        label=valueToShow, 
        nudge_x = 0.25, nudge_y = 0.25, 
        check_overlap = T
      )+ geom_hline(yintercept=0,color="black")+ geom_vline(xintercept = 0,color="black")
    
     }
}

#Fonction trouvée sur internet
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

showCorrelFun <- function(data,xvalue=1,yvalue=2,number=TRUE){
  
  dat <- circleFun(c(0,0),2,npoints = 700)
  valueToShow <- c(1:dim(data)[1])
  if(!number){
    valueToShow <- rownames(data)
  }
  
  
  #affichage du cercle de corrélation
  library(ggplot2)
  library(ggrepel)
  
  variable.graphe <- ggplot(data, aes(data[,xvalue], data[,yvalue])) +
    ggtitle("Cercle des corrélations des variables")+
    theme(plot.title = element_text(hjust = 0.5,face="bold"))+
    geom_hline(yintercept=0,color="black")+
    geom_vline(xintercept = 0,color="black")+
    geom_point(color = "blue", size = 3)+
    geom_point(aes(x=x, y=y), data=dat,color="black")
  
  variable.graphe + geom_segment(aes(x = 0, y = 0, xend = data[,xvalue], yend = data[,yvalue]),
                                 arrow = arrow(length = unit(0.5, "cm")))+
    geom_label_repel(aes(label = valueToShow),
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
  rownames(dataCor)<- colnames(dataVal)
  
  
  # CP qualité des variables
  dataCor.cos2 <- as.data.frame((G[,1:nb_axis]^2)*100,row.names = colnames(dataVal))

  # CP contribution des variables
  dataCor.contr <- dataCor.cos2 * 100
  sums <- colSums(dataCor.cos2)
  for(i in 1:nb_axis){
    dataCor.contr[,i] <- dataCor.contr[,i] / sums[i] 
  }
  dataCor.contr <- as.data.frame(dataCor.contr,row.names = colnames(dataVal))
  
  
  listOfDataframe = list(
    "datas" = dataVal,
    "indiv" = data,
    "var.coord" = dataCor,
    "var.cos2" = dataCor.cos2,
    "var.contr" = dataCor.contr
  )
  
  return(listOfDataframe)
}
