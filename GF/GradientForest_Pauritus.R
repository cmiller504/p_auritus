library(gradientForest)

Efrog=read.csv("env_GF.csv",header=T) #Environmental data at sites
Gfrog=read.table("frog_grqs",header=T,sep="",na.strings="NA",fill=T) 
Gfrog = Gfrog[ , colSums(is.na(Gfrog)) == 0] #SNPdata [with no -nan values (sed -i 's/-nan/0/g' frog_snps) and cols with NAs removed]


frogsites=Efrog[,1:2]

preds <- colnames(Efrog)
specs <- colnames(Gfrog)

nSites <- dim(Gfrog)[1]
nSpecs <- dim(Gfrog)[2]

# set depth of conditional permutation
lev <- floor(log2(nSites*0.368/2))

lev

frogforest=gradientForest(cbind(Efrog,Gfrog), predictor.vars=preds, response.vars=specs, ntree=2000, transform = NULL, compact=F,nbin=100, maxLevel=0,trace=F)

#predictoroverallimportance

plot.gradientForest(frogforest,plot.type="Overall.Importance")

predictors=names(importance(frogforest))


#splitsdensityplots
plot(frogforest, plot.type="S", imp.vars=predictors, leg.posn="topright", cex.legend=0.4, cex.axis=0.6, cex.lab=0.7, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), mar=c(3.1,1.5,0.1,1)))

#predictorcumulative
plot(frogforest, plot.type="C", imp.vars=predictors, show.species=F, common.scale=T, cex.axis=0.6, cex.lab=0.7, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), mar=c(2.5,1.0,0.1,0.5), omi=c(0,0.3,0,0)))

#R2
plot(frogforest, plot.type="P", show.names=F, horizontal=F, cex.axis=1, cex.labels=0.7, line=2.5)

# read in grid of points
froggrid=read.table("input_files/random_grid_EVs_100000_below7.csv",header=T,sep=",")


#transform grid and environmental predictors
predictors <- names(importance(frogforest)[1:3])
tgrid=cbind(froggrid[,c("long","lat")], predict.gradientForest(frogforest,froggrid[,predictors]))
Trns_site <- predict(frogforest)

#pcs
PCs=prcomp(tgrid[,3:5])
sgn <- sign(PCs$rotation[1:3])
PCs$rotation <- sweep(PCs$rotation,2,sgn,"*")
PCs$x <- sweep(PCs$x,2,sgn,"*")
# set up a colour palette for the mapping
a1 <- PCs$x[,1]
a2 <- PCs$x[,2]
a3 <- PCs$x[,3]
r <- a1+a2
g <- -a2
b <- a3+a2-a1
r <- (r-min(r)) / (max(r)-min(r)) * 255
g <- (g-min(g)) / (max(g)-min(g)) * 255
b <- (b-min(b)) / (max(b)-min(b)) * 255

nvs <- dim(PCs$rotation)[3] # number of variables
vec <- c("bio19", "lat" ,"bio15") 
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
# choose a scaling factor to plot the vectors over the grid
scal <- 60
xrng <- range(PCs$x[,1], PCs$rotation[,1]/scal)*1.1
yrng <- range(PCs$x[,2], PCs$rotation[,2]/scal)*1.1
plot((PCs$x[,1:2]), xlim=xrng, ylim=yrng, pch=".", cex = 10, col=rgb(r,g,b, max = 255), asp=1)
# plot the other predictors with "+"
points(PCs$rotation[! vind,1:2]/scal, pch="+")  
# plot the chosen predictors as arrows
arrows(rep(0,lv), rep(0,lv), PCs$rotation[,1]/scal, PCs$rotation[,2]/scal, length = 0.1,lwd = 4)
jit <- 0.0015
jit = 0.003
text(PCs$rotation[vec,1]/scal+jit*sign(PCs$rotation[vec,1]), PCs$rotation[vec,2]/scal+jit*sign(PCs$rotation[vec,2]), labels = vec, cex = 2, font = 2)

# first predict the PCs for the transformed site data
PCsites <- predict(PCs,Trns_site[,predictors])
# plot all the sites as points on the biplot
points(PCsites[,1:2])

#plot these in geographic space

frog.pred <- predict(frogforest,tgrid[,predictors])
plot(tgrid[,c("long","lat")],pch=15,cex=1.0,asp=1,col=rgb(r,g,b, max=255),main="SNP turnover")
points(frogsites[,2:1])


#export map for use in ArcGIS
frogcols=rgb(r,g,b,max=255)
frogcols2=col2rgb(frogcols)
frogcols3=t(frogcols2)
gradients=cbind(tgrid[c("long","lat")],frogcols3)
write.csv(gradients,file="results/frog_gradients4arcgis.csv")


frogtotal=frogforest$species.pos.rsq #this will give you the total SNPs that had an R^2 greater than 0
frogtotal
frogaverage=sum(frogforest$result)/frogtotal #this will give you the average R^2 value of those SNPs
frogaverage
range(frogforest$result)
median(frogforest$result)
