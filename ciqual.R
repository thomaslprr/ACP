library(readr)
data_ciqual <- read_delim("Desktop/Master Informatique/Analyse des données/Projet/data_ciqual.csv"
                          ,delim = ";", escape_double = FALSE, locale = locale(encoding = "latin1"),
                          trim_ws = TRUE,show_col_types = FALSE)

dim(data_ciqual)
# 3186 76

names(data_ciqual)
# voir le noms des colonnes

str(data_ciqual)
# voir les compositions des colonnes, on remarque leur format caractère

summary(data_ciqual)
# rien de vraiment intéressant car caractère

#ISOLATION DES DONNEES FROMAGES


sapply(data_ciqual, class) 

i <- c(10: 76)
data <- data_ciqual
data[ , i] <- apply(data_ciqual[ , i], 2,            # Specify own function within apply
                    function(x) as.numeric(sub(",", ".", as.character(x), fixed = TRUE)))

length(data[data=='fromages et assimilés'])
#135 lignes sur le fromages

install.packages("dplyr")
library("dplyr")
fromages <- dplyr::filter(data, grepl('fromages et assimilés', alim_ssgrp_nom_fr))
fromages <- fromages[-1,]
fromages 
summary(fromages) #bcp plus intéressant car on a des vrais données pour chaque colonne

#TRAITEMENT DES DONNEES NA 

sum(is.na(fromages))
#Nombre de données NA

colSums(is.na(fromages))
#Nombre de données NA par colonne

length(fromages[fromages=='-'])
#Nombre de données non renseignées

fromage <- fromages[,1:61]
drops <- c("Sélénium (µg/100 g)","Potassium (mg/100 g)","Phosphore (mg/100 g)",
           "Manganèse (mg/100 g)","Magnésium (mg/100 g)","Iode (µg/100 g)","Fer (mg/100 g)",
           "Cuivre (mg/100 g)","Chlorure (mg/100 g)","Cholestérol (mg/100 g)",
           "AG 22:6 4c,7c,10c,13c,16c,19c (n-3) DHA (g/100 g)",
           "AG 20:5 5c,8c,11c,14c,17c (n-3) EPA (g/100 g)",
           "AG 20:4 5c,8c,11c,14c (n-6), arachidonique (g/100 g)",
           "AG 18:3 c9,c12,c15 (n-3), alpha-linolénique (g/100 g)",
           "AG 4:0, butyrique (g/100 g)",
           "AG 6:0, caproïque (g/100 g)",
           "AG 8:0, caprylique (g/100 g)",
           "AG 10:0, caprique (g/100 g)",
           "AG 12:0, laurique (g/100 g)",
           "AG 14:0, myristique (g/100 g)",
           "AG 16:0, palmitique (g/100 g)",
           "AG 18:0, stéarique (g/100 g)",
           "AG 18:1 9c (n-9), oléique (g/100 g)",
           "AG 18:2 9c,12c (n-6), linoléique (g/100 g)",
           "Fructose (g/100 g)",
           "Galactose (g/100 g)",
           "Glucose (g/100 g)",
           "Lactose (g/100 g)",
           "Maltose (g/100 g)",
           "Saccharose (g/100 g)",
           "Amidon (g/100 g)",
           "Alcool (g/100 g)",
           "alim_nom_sci"
           )

fromage <- fromage[ , !(names(fromage) %in% drops)]

colSums(is.na(fromage))

fromage <- fromage[!is.na(fromage$`Energie, Règlement UE N° 1169/2011 (kJ/100 g)`),]
fromage <- fromage[!is.na(fromage$`Calcium (mg/100 g)`),]
fromage <- fromage[!is.na(fromage$`Sodium (mg/100 g)`),]
fromage <- fromage[!is.na(fromage$`AG saturés (g/100 g)`),]
fromage <- fromage[!is.na(fromage$`Eau (g/100 g)`),]
fromage <- fromage[!is.na(fromage$`AG polyinsaturés (g/100 g)`),]
#Amidon, sucre ????
fromage <- fromage[!is.na(fromage$`Sucres (g/100 g)`),]

colSums(is.na(fromage))

names <- fromage[,1:8]
fro <- fromage[9:28]


round(colMeans(fro),2) #on en déduit que la matrice n'est pas centrée réduite
X <- scale(fro,center = TRUE, scale=TRUE)*sqrt(89/88) #centrée réduite
S = t(X)%*%X * (1/89)  #matrice de variance/covariance
R <- t(S)%*%S * (1/89) #matrice de corrélation

nb_axe <- round(eigen(S)$values*5,3)
round(cumsum(eigen(S)$values*5),3);nb_axe  #somme des valeurs des variances cumulées
#affichage sous forme de graphique
barplot(main = "Variance cumulée",
        col="blue", 
        nb_axe) 

# 4. ACP N diagonalisation de R

ACP=eigen(R) ; ACP

# inertie Ig
round(ACP$values,3) ; sum(ACP$values) ; sum(diag(R)) 

# axes principaux
u=ACP$vectors;u

# composantes principales F : individus
fr <- data.matrix(fro);fr
F=fr%*%u ; F

# plan 1-2 des individus
plot(F[,1],F[,2],xlab="Axe 1", ylab="Axe 2")


# G : CP coord des variables
G<-matrix(0,20,20)
G[,1]=sqrt(ACP$values[1])%*%ACP$vectors[,1]
G[,2]=sqrt(ACP$values[2])%*%ACP$vectors[,2]
round(G,3)
plot(G)
#s.corcircle()
s.corcircle (G[,1:2])


###version ave ADE4###

library(ade4)
acp<-dudi.pca(fro,center=TRUE,scale=TRUE,scannf=TRUE)
round(acp$eig,2)
round(cumsum(acp$eig*5),2)

#ON ANALYSE LES VARIABLES
inertie <-inertia.dudi(acp, col.inertia=TRUE) 
round(acp$co,2)

#contribution axe 1 et axe 2 
round(inertie$col.abs,2)

#qualité de représentation des variables
round(inertie$col.rel,2) 

s.corcircle(acp$co,xax=1,yax=2)

#Analyse des individus 
inertie <-inertia.dudi(acp, row.inertia=TRUE)
round(acp$li,2) 

#Contribution individu
round(inertie$row.abs,2)

#Qualité de représentation
round(inertie$row.rel,2)

#Affichage graphique
s.label(acp$li,xax=1,yax=2) 
