#installation du package ggrepel
install.packages("ggrepel")

#installation du package dplyr
install.packages("dplyr")

#installation du package funModeling
install.packages("funModeling")

#installation du package naniar
install.packages("naniar")

#chargement de la librairie de lecture csv
library(readr)
#import des données
data_ciqual <- read_delim("Desktop/Master Informatique/Analyse des données/Projet/data_ciqual.csv"
                          ,delim = ";", escape_double = FALSE, locale = locale(encoding = "latin1"),
                          trim_ws = TRUE,show_col_types = FALSE)

#dimension du dataset
dim(data_ciqual)

#nom des colonnes
names(data_ciqual)

#compositions des colonnes, on remarque leur format caractère
str(data_ciqual)

# rien de vraiment intéressant car les colonnes sont des caractères
summary(data_ciqual)



###ISOLATION DES DONNEES FROMAGES###

#copie des données originales par sécurité
data <- data_ciqual

#classes de chaque colonne
sapply(data_ciqual, class) 

#nombre de lignes sur le fromage
length(data[data=='fromages et assimilés'])

#isolation des colonnes qui n'ont pas le bon format
i <- c(10: 76)

#changement du type des colonnes de caractère à numérique 
# +
#remplacement des virgules par des points
data[ , i] <- apply(data_ciqual[ , i], 2, 
                    function(x) as.numeric(sub(",", ".", as.character(x), fixed = TRUE)))

#import de la librairie
library("dplyr")
#on garde seulement les fromages et assimilés
fromages <- dplyr::filter(data, grepl('fromages et assimilés', alim_ssgrp_nom_fr))

#on sauvegarde l'individu supplémentaire
indiv_supp <- fromages[1,]

#suppression du fromage moyen
fromages <- fromages[-1,]
fromages 
#sommaire des fromages
#bcp plus intéressant qu'avant car on a des données numériques pour chaque colonne
summary(fromages) 





### TRAITEMENT DES DONNEES NA ###
#Nombre de données NA
sum(is.na(fromages))

#Nombre de données NA par colonne
colSums(is.na(fromages))

library(funModeling)
df_status(fromages)
library(naniar)
missVar.p1 <- fromages[,1:38]
missVar.p2 <- fromages[,39:76]

gg_miss_var(missVar.p1,show_pct = TRUE)
gg_miss_var(missVar.p2,show_pct = TRUE) 

#Nombre de données non renseignées
length(fromages[fromages=='-']) 
#lors de la conversion des colonnes, les "-" ont été remplacé par des NA ce qui explique
#qu'on obtient 0 

#Elimination des colonnes "inutiles"
fromages <- fromages[,1:61]
indiv_supp <- indiv_supp[,1:61]
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

fromages <- fromages[ , !(names(fromages) %in% drops)]

indiv_supp <- indiv_supp[ ,!(names(indiv_supp) %in% drops)]

#observation des NA par colonne après le premier balayage
colSums(is.na(fromages))


#Elimination des lignes "inutiles"
fromages <- fromages[!is.na(fromages$`Energie, Règlement UE N° 1169/2011 (kJ/100 g)`),]
fromages <- fromages[!is.na(fromages$`Calcium (mg/100 g)`),]
fromages <- fromages[!is.na(fromages$`Sodium (mg/100 g)`),]
fromages <- fromages[!is.na(fromages$`AG saturés (g/100 g)`),]
fromages <- fromages[!is.na(fromages$`Eau (g/100 g)`),]
fromages <- fromages[!is.na(fromages$`AG polyinsaturés (g/100 g)`),]

#Dilemme entre suppression des lignes ou il n'y a pas de sucre renseigné
#ou des lignes ou il n'y a pas d'amidon renseigné
#Suppression des lignes où la valeur du sucre n'est pas renseigné
fromages <- fromages[!is.na(fromages$`Sucres (g/100 g)`),]

#observation des NA par colonne après le balayage par colonne puis par ligne
colSums(is.na(fromages))

#isolation des données nominatives
names <- as.data.frame(fromages[,8])
#isolation des données quantitatives
fro <- fromages[9:28]
indiv_supp <- indiv_supp[9:28]

rownames(fro) <- names[,1]

### ACP ###

acpMan <- acpFun(fro)
# plan 1-2 des individus
showIndivFun(acpMan$indiv,label=FALSE)
showIndivFun(acpMan$indiv,label=FALSE,number=FALSE)
showIndivFun(acpMan$indiv)


#affichage du cercle de corrélation
showCorrelFun(acpMan$var.coord)
showCorrelFun(acpMan$var.coord,number=FALSE)



### VERSION AVEC ADE4 ###

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
ligsup<-suprow(acp,indiv_supp)

inertie <-inertia.dudi(acp, row.inertia=TRUE)
round(acp$li,2) 

#Contribution individu
round(inertie$row.abs,2)

#Qualité de représentation
round(inertie$row.rel,2)

#Affichage graphique sans l'individu supplémentaire
s.label(acp$li,xax=1,yax=2) 
#Affichage graphique avec l'individu supplémentaire
#coordonnées des individus actifs
cl1<-acp$li[,1]
cl2<-acp$li[,2]
#coordonnées des individus supplémentaires
csup1<-ligsup$lisup[,1]
csup2<-ligsup$lisup[,2]
#le graphique "vide"
plot(cl1,cl2,type="n",main="Les individus",xlim=c(-8,8))
abline(h=0,v=0)
#on ajoute les individus actifs
text(cl1,cl2,row.names(acp$li),)
#on ajoute les individus supplémentaires
text(csup1,csup2,row.names(ligsup$lisup),col="red",cex=1.2)

#isolement de l'individu supplémentaire
plot(cl1,cl2,type="n",main="Les individus",xlim=c(-8,8))
abline(h=0,v=0)
text(csup1,csup2,row.names(ligsup$lisup),col="red",cex=1.2)
