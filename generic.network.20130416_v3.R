# 2013 April 16
# Compare power-law, poisson, and regular network model of aging. 

rm(list=ls());

power_numbers = function( k0, kN, gamma) {
  #k0=2; kN=50; gamma=2.5
  my.k= seq(k0, kN, 1)
  my.p = my.k^(-gamma)
  my.p = my.p / (sum(my.p))
  my.Frep = round( (kN-k0+1)*my.p + 0.5)
  return( rep(k0:kN, my.Frep) )
}
powerLawLinks = power_numbers (4,100,3)

calculate.s.m = function( lifespan, bin.size=10 ){ #lifespan is simulated data
   my.data = sort( lifespan[!is.na(lifespan)] );
   my.max = max(my.data)
   my.interval.size = bin.size
   deathFreq = c(NA, NA)
   for ( i in 1:round(length(my.data) / my.interval.size))  {
     sub = my.data[ ((i-1)* my.interval.size + 1) : (i* my.interval.size) ]
     deathFreq = rbind( deathFreq, c(mean(sub, na.rm=T), my.interval.size) )
   }
   deathFreq = data.frame(deathFreq)
   deathFreq = deathFreq[-1,] #remove NA
   #deathFreq       = table( my.data )
   deathCumulative = deathFreq
   for( i in 2:length(deathFreq[,1])) {
        deathCumulative[i, 2] = deathCumulative[i-1, 2] + deathCumulative[i, 2]                
   }
   tot = length(my.data)
   s = 1 - deathCumulative[,2]/tot
   currentLive  = tot - deathCumulative[,2]        
   m =  deathFreq[,2] / currentLive; 

   #list( s=s, t=unique(my.data));
   #ret = data.frame( cbind(s, m, unique(my.data)));
   ret = data.frame( cbind(s, m, deathFreq[,1]) );
   names(ret) = c("s", "m", "t");
   ret;
}

#test
lifespan = rnorm(1000, mean=100, sd=10)
ret = calculate.s.m(lifespan)

#mm modules
#nn elements in each module
simulate.age.single.population = function(mm, nn, Npop, meanMu) {
  SystemAges = numeric(Npop);#store the lifespan for all individuals 
  BlockAges = numeric(mm) #buffer for temporary storage
  for( i in 1: Npop){  # i-th individual in the population  
   for( j in 1:mm) {#Block loop
    #mu.vec strores the constant failure rates of elements
    #We are not sure whether mu.vec should be inside of Nop loop
    #mu.vec = rlnorm(n[j], mean=0.005, sd=1) #Element constant mortality rates
    
    #In GG01 paper, mu.vec are the same. 
    mu.vec = rep(meanMu, nn[j])
    ElementAges =  rexp(nn[j], rate=mu.vec);    
    
    #maximum of elelementAges -> BlackAges
    BlockAges[j] = round( max(ElementAges), 2 );  
  }
  SystemAges[i] = min(BlockAges)  
 }
 return( SystemAges )
}


##################
NSims = 50
mm = 15;  # numOfBlocks in a system
Npop = 1E2; # numOfSystems (individuals), ie population size
meanMu = 2

#meanNN = 5; #mean number of elements per module
meanNN = mean(powerLawLinks)
medianNN = median(powerLawLinks)

nn = rep(meanNN,mm) ; #regular network using mean Powerlaw
CV_regular = numeric(NSims)
for ( sim in 1:NSims) {
  SystemAges = simulate.age.single.population(mm, nn, Npop, meanMu ); 
  CV_regular[sim] = sd(SystemAges ) / mean( SystemAges ) 
}
CV_regular

nn = rep(medianNN,mm) ; #regular network using median Powerlaw
CV_regular2 = numeric(NSims)
for ( sim in 1:NSims) {
  SystemAges = simulate.age.single.population(mm, nn, Npop, meanMu ); 
  CV_regular2[sim] = sd(SystemAges ) / mean( SystemAges ) 
}
CV_regular2

#powerlaw
CV_powerlaw = numeric(NSims)
for ( sim in 1:NSims) {
  SystemAges = simulate.age.single.population(mm, powerLawLinks, Npop, meanMu ); 
  CV_powerlaw[sim] =  sd(SystemAges ) /mean( SystemAges )
}
CV_powerlaw

#Poisson network, using powerLaw mean
nn = rpois(mm, meanNN); # Poisson network
CV_poisson = numeric(NSims)
for ( sim in 1:NSims) {
  SystemAges = simulate.age.single.population(mm, nn, Npop, meanMu ); 
  CV_poisson[sim] =  sd(SystemAges ) /mean( SystemAges )
}
CV_poisson

#Poisson network, using powerLaw median
nn = rpois(mm, medianNN); # Poisson network
CV_poisson2 = numeric(NSims)
for ( sim in 1:NSims) {
  SystemAges = simulate.age.single.population(mm, nn, Npop, meanMu ); 
  CV_poisson2[sim] =  sd(SystemAges ) /mean( SystemAges )
}
CV_poisson2

tb = cbind(CV_powerlaw, CV_regular, CV_regular2, CV_poisson, CV_poisson2)
summary(tb)


#t.test( CV_regular, CV_poisson)

######################
#tb = calculate.s.m( SystemAges )
#head(tb)
#plot( tb$s ~ tb$t)
#plot(tb$m ~ tb$t)
#plot(log10(tb$m) ~ tb$t)

#sub = tb[1:floor(length(tb[,1])/4), ]
#sub = tb[ tb$s > 0.25, ]
#with( sub, plot( m ~ t))
#with( sub, plot( s~ t ))

#sub$m[sub$m==0] = NA;  #remove zero for log operations
#plot (log(sub$m) ~ sub$t, main = paste("sd=", mysd) ) #plot the mortality ~ age

  
