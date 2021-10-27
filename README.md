# ACP avec R

```fromages.csv``` correspond au fichier des données utilisées pour faire fonctionner notre fonction

```data_ciqual.csv``` correspond au fichier des données de base (notre projet travaille à partir de ce fichier qui est celui d'origine sur Madoc)

```acp_function.R``` correspond au fichier avec l'ensemble des fonctions (ACP, affichages)

```ciqual.R``` correspond au fichier principal. Il contient une fonction de nettoyage des données, 
puis exécute la démonstration de notre fonction ACP et de ses usages ainsi qu'une démonstration d'ACP via ADE4

## Prérequis

Installation du package ggrepel : ```install.packages("ggrepel")``` (obligatoire pour la fonction ACP)

Installation du package dplyr : ```install.packages("dplyr")``` (obligatoire pour le nettoyage des données) 

Installation du package funModeling : ```install.packages("funModeling")``` (obligatoire pour le nettoyage des données)

Installation du package naniar : ```install.packages("naniar")``` (obligatoire pour le nettoyage des données)

Installation du package ggplot2 : ```install.packages("ggplot2")``` (obligatoire pour la fonction ACP)

Installation du package ade4 : ```install.packages("ade4")``` (obligatoire pour utiliser ade4)

## Utilisation 

### Etape 1 
Il vous faut simplement exécuter tout le code du fichier ```acp_function.R``` 

**NB : Libraries GGREPEL, GGPLOT2 sont nécessaires.**

### Etape 2
• Ouvrir le fichier ```ciqual.R``` 

• Installer les packages nécessaires *(dplyr, funModeling, naniar)*  

• Changer l'adresse source du fichier à charger dans la fonction ```getFro()``` selon là où se situe votre jeu de données

• Compiler la fonction ```getFro()```

• Exécuter le reste du script au fur et à mesure selon les besoins


## Composition du groupe
Lucas LELIEVRE | Thomas LAPIERRE
