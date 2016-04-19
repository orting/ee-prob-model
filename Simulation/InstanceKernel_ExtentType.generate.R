nBags <- 10
bagSize <- 100
nInstances <- nBags*bagSize

## Model patches as an average intensity
emphysemaPopulationProportion <- 0.3
nEmphysema <- ceiling( nInstances * emphysemaPopulationProportion )
nCentrilobular <- ceiling( nEmphysema / 2 )
nParaseptal <- nEmphysema - nCentrilobular
nHealthy <- nInstances - nEmphysema

centrilobular <- rgamma( nCentrilobular, shape=6, scale=10 ) - 1000
paraseptal <- rgamma( nParaseptal, shape=5, scale=10 ) - 1000 
healthy <- rgamma( nHealthy, shape=10, scale=10 ) - 1000;

instanceLabels <- rep(0, nInstances)
emphysemaInstances <- sample(nInstances, nEmphysema)
paraseptalInstances <- emphysemaInstances[ sample( nEmphysema, nCentrilobular ) ]
instanceLabels[ emphysemaInstances ] <- 1
instanceLabels[ paraseptalInstances ] <- 2

instances <- rep(0, nInstances)
instances[ instanceLabels == 0 ] <- healthy
instances[ instanceLabels == 1 ] <- centrilobular
instances[ instanceLabels == 2 ] <- paraseptal

K <- exp( - as.matrix( dist( instances, diag=T, upper=T, 'euclidean' ) ) )
dimnames(K) <- c()
xi <- rep(100,nBags)
xiTy <- 100
Ty <- matrix(nrow=nBags, ncol=3)
P <- rep(0,nBags)
for ( i in 1:nBags ) {
    low <- 1 + (i-1)*bagSize
    high <- i*bagSize
    for ( j in 1:3 ) {
        Ty[i,j] <- sum( instanceLabels[low:high] == (j-1) ) / bagSize
    }
    P[i] <- sum( instanceLabels[low:high] != 0 ) / bagSize
}
dump( c('nBags', 'bagSize', 'nInstances', 'K', 'P', 'Ty', 'xi', 'xiTy', 'instanceLabels'), 'InstanceKernel_ExtentType.data.R')
