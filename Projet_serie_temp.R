# ==============================================================================
# PROJET DE SÉRIES TEMPORELLES : MODÉLISATION ARIMA
# Membres du binôme : [Prénom Nom 1] & [Prénom Nom 2]
# Série : Indice de la Production Industrielle (IPI) - Métallurgie
# Source : INSEE (Série 010767941)
# ==============================================================================

# Chargement des bibliothèques nécessaires
if(!require(zoo)) install.packages("zoo")
if(!require(tseries)) install.packages("tseries")
library(zoo)
library(tseries)

# ==============================================================================
# PARTIE I : LES DONNÉES
# ==============================================================================

# ------------------------------------------------------------------------------
# QUESTION 1 : Représentation de la série et traitements éventuels
# ------------------------------------------------------------------------------

# IMPORTATION
# La série brute est l'Indice de la Production Industrielle (base 100 en 2021) 
# pour le secteur "Métallurgie des autres métaux non ferreux".
# Le fichier INSEE contient 3 lignes d'en-tête, nous utilisons skip=3.
data_raw <- read.csv("valeurs_mensuelles.csv", sep=";", dec=".", skip=4, header=FALSE)
colnames(data_raw) <- c("Date", "IPI", "Code")

# Conversion en objet temporel "zoo" (fréquence mensuelle)
dates <- as.yearmon(data_raw$Date, format="%Y-%m")
X <- zoo(as.numeric(data_raw$IPI), order.by=dates)
X <- na.omit(X)

# JUSTIFICATION DES TRAITEMENTS :
# CVS-CJO : La série choisie est déjà corrigée des variations saisonnières et 
#    des jours ouvrables par l'INSEE. Il n'est donc pas nécessaire d'appliquer 
#    de différenciation saisonnière (D=0) ou de modèle SARIMA.
#  TRANSFORMATION LOG : Comme la série couvre une longue période (1990-2026), 
#    on observe une légère augmentation de la volatilité avec le niveau de l'indice. 
#    Nous appliquons une transformation logarithmique Y_t = log(X_t) pour 
#    stabiliser la variance (transformation de Box-Cox avec lambda = 0).

logX <- log(X)

# ------------------------------------------------------------------------------
# QUESTION 2 : Transformation pour rendre la série stationnaire
# ------------------------------------------------------------------------------

# ANALYSE DE LA STATIONNARITÉ (Série en Log)
# Nous testons la présence d'une racine unitaire sur log(X_t).
adf_log <- adf.test(logX)
print(adf_log)
# Si p-value > 0.05 : On ne rejette pas l'hypothèse de racine unitaire (DS).

# DIFFÉRENCIATION
# Pour stationnariser la série en moyenne, nous appliquons l'opérateur de  différence première : dlogX_t = (1 - L) log(X_t) = log(X_t) - log(X_{t-1}).

dlogX <- diff(logX, differences = 1)

# VÉRIFICATION DE LA STATIONNARITÉ (Série différenciée)

# Test ADF (H0 : racine unitaire / non-stationnaire)
adf_test <- adf.test(dlogX)
print(adf_test) # On attend une p-value < 0.05 pour rejeter H0.

# Test KPSS (H0 : stationnarité) - Pour confirmer par un test inverse
kpss_test <- kpss.test(dlogX, null="Level")
print(kpss_test) # On attend une p-value > 0.05 pour ne pas rejeter H0.

# Conclusion : La série dlogX est stationnaire au second ordre, notée I(0).

# ------------------------------------------------------------------------------
# QUESTION 3 : Représentations graphiques (Avant et Après)
# ------------------------------------------------------------------------------

# On s'assure que l'affichage n'est pas divisé en plusieurs morceaux
par(mfrow=c(1,1)) 

# La série originale (Brute)
plot(X, main="Graphique 1 : Série brute de l'IPI (Jan. 1990 - Déc. 2019)", 
     ylab="IPI (Base 100)", xlab="Temps", col="royalblue", lwd=2)

# La série avec le Log
plot(logX, main="Graphique 2 : Série au logarithme", 
     ylab="Log(IPI)", xlab="Temps", col="darkorange", lwd=2)

#La série différenciée (Stationnaire)
plot(dlogX, main="Graphique 3 : Série stationnarisée (Différence du Log)", 
     ylab="Taux de croissance", xlab="Temps", col="forestgreen", lwd=2)
abline(h=mean(dlogX), col="red", lty=2, lwd=2) # Ajoute la ligne rouge de moyenne


# ==============================================================================
# PARTIE II : MODÈLES ARMA 
# ==============================================================================

# ------------------------------------------------------------------------------
#IDENTIFICATION GRAPHIQUE DES ORDRES (p, q)
# -----------------------------------------------------------------------------

# On centre la série ( y <- desaison - mean(desaison))
y <- dlogX - mean(dlogX)

par(mfrow=c(1,2)) 

# ACF et PACF purs. 

acf(as.numeric(y), 36, main="ACF de la série centrée")
pacf(as.numeric(y), 36, main="PACF de la série centrée")

par(mfrow=c(1,1))

# ------------------------------------------------------------------------------
# Estimation des paramètres 
# -----------------------------------------------------------------------------

# Graphiquement,grâce a l'ACF et le PACF, on remarque bien que q=2 et p=3
modele_AR3 <- arima(y, order=c(3,0,0))
modele_MA2 <- arima(y, order=c(0,0,2))
modele_ARMA11 <- arima(y, order=c(1,0,1))
modele_ARMA21 <- arima(y, order=c(2,0,1))

#tableau pour comparer les critères d'information (AIC et BIC)

criteres <- data.frame(
  Modele = c("AR(3)", "MA(2)", "ARMA(1,1)", "ARMA(2,1)"),
  AIC = c(AIC(modele_AR3), AIC(modele_MA2), AIC(modele_ARMA11), AIC(modele_ARMA21)),
  BIC = c(BIC(modele_AR3), BIC(modele_MA2), BIC(modele_ARMA11), BIC(modele_ARMA21))
)


print(criteres)

