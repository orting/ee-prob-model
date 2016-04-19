/* Emphysema estimation with Learning with Label Proportions (LLP)
 * We have a set of lung CT scans with estimates of emphysema extent in six
 * regions of the lung and estimate of predominant emphysema sub-type in the
 * scan.
 * We want to categorize patches extracted from the regions as one of:
 *  Healty
 *  Centrilobular emphysema
 *  Paraseptal emphysema
 *
 * The model here is based on the formulation from 
 *  Hendrik Kück and Nando de Freitas
 *  Learning about individuals from group statistics
 */
data {
  int<lower=0> nBags;
  int<lower=0> bagSize;
  int<lower=nBags*bagSize,upper=nBags*bagSize> nInstances;

  // Kernel matrix, where K[i,j] = kernel(instance i, instance j)
  //  row_vector<lower=0>[nInstances] K[nInstances];
  matrix<lower=0>[nInstances,nInstances] K;

  // Estimated extent of emphysema in region AKA proportion of patches with
  // emphysema. Set this to the midpoint of the combined interval of the raters.
  // NOTE: Must not be 0 or 1
  real<lower=0,upper=1> P[nBags];

  // Estimated predominant type of emphysema.
  //  Index 0 = Not emphysema
  //  Index 1 = Centrilobular
  //  Index 2 = Paraseptal
  // Mixed is coded as an equal  mixture of centrilobular and paraseptal.
  // Disagreement between raters is coded as a mixture of the types.
  // e.g. Two raters, one say mixed and one say centrilobular is coded as
  //  (0, 0.5/2 + 1/2, 0.5/2) = (0, 0.75, 0.25)
  // NOTE: No entry may be 0 or 1
  simplex[3] Ty[nBags]; 

  // Certainty of extent label, higher is more certain
  // This should be set so a desired percentage of the density from
  //   beta( xi*IntervalMidpoint + 1, xi*(1 - IntervalMidpoint) + 1)
  // lies inside the interval.
  real<lower=1> xi[nBags];

  // Certainty of the type labels. This is a global certainty because Ty already
  // includes uncertainty about the specific labels.
  real<lower=1> xiTy;
}

transformed data {
  // Division of instances into bags
  int<lower=0,upper=nInstances> I[nBags,2];

  for ( n in 1:nBags ) {
    I[n,1] <- 1 + (n-1)*bagSize;
    I[n,2] <- n*bagSize;
  }
}
  
parameters {
  // Weights for the kernels
  matrix<lower=0>[3,nInstances] w;

  // Predicted probability of tissue type for each instance
  //  Index 0 = Not Emphysema
  //  Index 1 = Centrilobular
  //  Index 2 = Paraseptal  
  simplex[3] z[nInstances];
}

transformed parameters {
  // Parameters for the z generating Dirichlet distribution
  matrix<lower=0>[3,nInstances] alpha;

  // lambda is the "true" proportion of emphysema patches in each bag
  real<lower=0,upper=1> lambda[nBags];

  // tau is the "true" distribution of tissue types in each bag
  vector<lower=0>[3] tau[nBags];
  
  alpha <- w*K;

  for ( n in 1:nBags ) {
    real lambdaSum; lambdaSum <- 0;
    for ( i in 1:3 ) {
      tau[n][i] <- 0;
    }
    for ( i in I[n,1]:I[n,2] ) {
      lambdaSum <- lambdaSum + sum(z[i,2:3]);
      tau[n] <- tau[n] + z[i];
    }
    lambda[n] <- lambdaSum / bagSize;
  }  
}

model {
  # prior
  for ( i in 1:nInstances ) {
    w[,i] ~ gamma(1,1);
  }

  for ( i in 1:nInstances ) {
    z[i] ~ dirichlet( alpha[,i] );
  }
  
  for ( i in 1:nBags ) {
    P[i] ~ beta(xi[i]*lambda[i] + 1.0, xi[i]*(1.0 - lambda[i]) + 1.0);
    Ty[i] ~ dirichlet( xiTy*tau[i] );
  }
}
