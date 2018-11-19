#1. - Visualize distribution of expression values across samples. Use function boxplot (). File = TableS5.csv
#Set directory where Table with genes VS cell type is and gather data.
#setwd("C:/Users/Maksym/Desktop/UPM/Genomic Data Analysis/")
mydata <- read.csv("TableS5.csv", header = TRUE, sep = ",", dec = ".", row.names = 1) 

#Visualize it. You can´t see well because of the outliers. Most of FPKMs are almost at 0.
boxplot (mydata)
boxplot (mydata, outline=FALSE)



##2.-Compare PCA separation and read distribution
#Is there any possible relationship between expression values distribution and PCA separation?


PCA<-princomp(mydata,cor=FALSE,scores=TRUE)
PCA[1:7]
plot(PCA$loadings)
text(PCA$loadings,labels=rownames(PCA$loadings),pos=2,cex=0.5)

#You can see that S18 and S30 are clearly separated from the rest in the plane.
#You might not need many number of changes to separate the cell types.

plot(PCA$loadings, main = "Ground Tissue Sorted Cells")
text (PCA$loadings, labels =rownames(PCA$loadings), cex = 0.6)

twocomps = PCA$scores[,1:2]

#Si quieres usar dos graficas a la vez, utiliza funcion par
?
par(mfrow=c(2,1)) #dos filas con 1 columna
boxplot (mydata)
plot(PCA$loadings, main = "Ground Tissue Sorted Cells")
text (PCA$loadings, labels =rownames(PCA$loadings), cex = 0.6)

#3.- Find biomarkers for Stem Cells as those only (WOX = Stemm cells)
#expressed in WOX5 domain (threshold=1). Use function
#subset (). Obtain biomarker names and visualize
#expression patterns in a heatmap.


#Tus muestras mas variables son las 2 componentes de PCA mas separadas, pero podemos encontrar 
#biomarcadores para celulas madre como aquellos expresados solo en VOX5 (marcador de celulas madres)
#Si encuentras genes especificos de celulas madres, vas a poder separar celulas que tienen mas
#probabilidad de ser celulas madres. Define cuantas veces los biomarcadores tienen que ser 
#MAS expresadas relativamente a los otros.
#Coge los valores genes con expresion >1 en WOX5 y los demas <1

marcadores <- subset(mydata, WOX5>1 & APL<1 &	CO2<1	& COBL9	<1 & COR	<1 & E30<1 & GL2<1 &	PET111<1 &	S17	<1 & S18<1 &	S32<1 &	S4<1 &	SCR<1 &	WER	<1 & WOL <1)
head(marcadores)

boxplot(marcadores, outline=FALSE)


dim(marcadores)#These 417 genes would be our biomarkers
rownames(marcadores)


#Heatmap function needs a matrix as input
marcadoresmatrix <-as.matrix(marcadores)

library(RColorBrewer)
newcolors<-colorRampPalette(colors=c("blue","black","yellow"))(256)
heatmap(marcadoresmatrix,scale="row",col =newcolors)


##4.- Separate samples by PCA using Stem Cell Biomarkers
par(mfrow=1)
PCAmarcadores<-princomp(marcadores,cor=FALSE,scores=TRUE)
summary(PCAmarcadores)
plot(PCAmarcadores$loadings)
text(PCAmarcadores$loadings,labels=rownames(PCAmarcadores$loadings),pos=3,cex=0.5)
#We also separate SCR, it is ok because it contains stem cells
windows()
plot(PCAmarcadores$scores)# few genes mosty contribute to SCR sepration 


##5.-Select biomarkers based on PCA scores, it is complicated

#Let's check what the scores of our stem biomarkers are in original dataset
summary(marcadores)
par(mfrow=c(1,5))
plot(PCA)
plot(PCA$scores)
plot(PCA$scores[rownames(marcadores),])#bad separation
plot(PCA$scores[,c(3:4)])
plot(PCA$scores[rownames(marcadores),c(3:4)])#better separation
