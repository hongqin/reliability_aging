# 2013 April 15
# Regular network model of cellular aging 
# Ref: Gavrilov & Gavrilova, 2006 Handbook of models of for human aging

rm(list=ls());
single.component.vialbility = function(mu, t) { exp( - mu * t) }
module.viability  = function(nn,mu,t) { 1 - (1 - exp(-mu*t)^nn ) } #n elements per module
system.viability =  function(mm, nn,mu,t) {

  
}
  
mm=10; #10 modules
nn = 10; #10 node simulation
mu = 10 #constant decaying rate
#mu = rpois(n, 10); #mean of exponetial distribution rexp( n, lambda)

maxmu = max(mu);
#t = seq(0, 7*maxmu, by= maxmu/50)




