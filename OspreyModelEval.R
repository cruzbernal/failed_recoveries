###############################################################################
###                  Code created by Jen Cruz                         ########
#### This script contains methods for calculating model evaluation stats  ####
####            based on likelihood of observed and predicted data      ######
####                    given model p(data|model).                        ####
### Model estimates nest occupancy, persistence, colonization, and       #####
### success for ospreys at Voyageurs National Park,  USA (1986-2012).    #####
##############################################################################

########### clean workspace and load required packages ####################
###########################################################################

#####clean workspace to improve efficiency: ###
rm(list = ls() ) 

####### load relevant packages ###
library( dplyr ) #dataframe manipulations.
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( tidyr ) #to use spread and other df functions
library( ggplot2 ) #fancy plots.
#library( lubridate ) #easy date adjustments and calculations.
library( jagsUI ) #to run RJAGS

########## end of package loading ###########

#######################    import relevant data   

#start with results from osprey model analysis (1986-2012):
load( "D:/RProjects/Osprey/DynOccSuc_Results.RData" )

##################################################################################################
## Calculating summary stats from likelihood of data conditional on model: p(y|theta)  ###
##################################################################################################
# Define model we want to extract output from 
mr <- m1
summary( mr )
#Total number of iterations ran:
S <- length( mr$sims.list$int.phi )#

# Define function to calculate model evaluation parameters 
EvaluateModel <- function( lik_yobs, lik_sucobs ){
  # create temporary objects
  CPO_yobs_ik <- CPO_sucobs_ik <- ll_yobs_ik <- ll_sucobs_ik <- var_yobs_ik <- var_sucobs_ik <- matrix( ncol = K, nrow = M )
  CPO_yobs_k <- CPO_sucobs_k <- ll_yobs_k <- ll_sucobs_k <- var_yobs_k <- var_sucobs_k <- rep( NA, K )
  
  for( k in 1:K ){ # loop over years
    for( i in 1:M ){ #loop over sites
      # Calculate the approximation of the Conditional Predictive Ordinate (CPO) across MCMC iterations 
      #for each site, each year:
      CPO_yobs_ik[ i, k ] <- S / sum( 1 / lik_yobs[ 1:S, i, k ] ) 
      CPO_sucobs_ik[ i, k ] <- S / sum( 1 / lik_sucobs[ 1:S, i, k ] ) 
      #Calculate components of the Watanabe-Akaike Information Criterion (WAIC):
      # Based on Appendix S1 equations from Broms et al. 2016.
      #Start by summing likelihood of data across MCMC iterations, dividing by MCMC steps, then logging it
      ll_yobs_ik[ i, k ] <- log( sum( lik_yobs[ 1:S, i, k ] ) / S )
      # Now for nest success data:
      ll_sucobs_ik[ i, k ] <- log( sum( lik_sucobs[ 1:S, i, k ] ) / S )
      #calculate variance between log likelihood computed at each MCMC step for each site, each
      #year and the mean log likelihood across MCMC steps:
      var_yobs_ik[ i, k ] <-  ( 1 / ( S - 1 ) ) *
                      sum( ( log( lik_yobs[ 1:S, i, k ] ) - log( mr$mean$lik_yobs[ i, k ] ) )^2 )
      var_sucobs_ik[ i, k ] <-  ( 1 / ( S - 1 ) ) *
                    sum( ( log( lik_sucobs[ 1:S, i, k ] ) - log( mr$mean$lik_sucobs[ i, k ] ) )^2 )
    }
    # Calculate the summary statistic of individual CPO values by summing over sites that were sampled each year inside VNP:
    CPO_yobs_k[ k ] <- - sum( log( CPO_yobs_ik[ 1:M, k ] ) *
                      XM$inVNP[ 1:M ] * surv_ind[ 1:M, k ]  )
    CPO_sucobs_k[ k ] <- - sum( log( CPO_sucobs_ik[ 1:M, k ] ) *
                              XM$inVNP[ 1:M ] * surv_ind[ 1:M, k ]  )
    # sum log likelihood only across sites that were sampled each year inside VNP:
    ll_yobs_k[ k ] <- sum( ll_yobs_ik[ 1:M, k ] *
                           XM$inVNP[ 1:M ] * surv_ind[ 1:M, k ] )
    ll_sucobs_k[ k ] <- sum( ll_sucobs_ik[ 1:M, k ] *
                             XM$inVNP[ 1:M ] * surv_ind[ 1:M, k ] )
    # sum variances only across sites sampled each year inside VNP:
    var_yobs_k[ k ] <- sum( var_yobs_ik[ 1:M, k ] *
                            XM$inVNP[ 1:M ] * surv_ind[ 1:M, k ] )
    var_sucobs_k[ k ] <- sum( var_sucobs_ik[ 1:M, k ] *
                              XM$inVNP[ 1:M ] * surv_ind[ 1:M, k ] )
  }

  # Calculate measures of predictive performance (CPO) for each model aspect:
  CPOyobs <- sum( CPO_yobs_k[ 1:K ] ) 
  CPOsucobs <- sum( CPO_sucobs_k[ 1:K ] )
  # Calculate elppd of WAIC:
  elppd <- sum( ll_yobs_k[ 1:K ] ) + sum( ll_sucobs_k[ 1:K ] )
  # Calculate model complexity component of WAIC following Gelman et al. 2014 pWAIC2:
  pWAIC2 <- sum( var_yobs_k[ 1:K ] ) + sum( var_sucobs_k[ 1:K ] )
  # Calculate WAIC
  WAIC <- -2 * elppd + 2 * pWAIC2 
  # Function output:
  return( data.frame( CPOyobs, CPOsucobs, elppd, pWAIC2, WAIC ) )
}
# Run function and store:
EvalRes <- EvaluateModel( mr$sims.list$lik_yobs, mr$sims.list$lik_sucobs )
EvalRes

# Define function to calculate Model deviances to be used in Bayesian p value calculations
CalcDevs <- function( ){
  # Define initial values:
  lik_yobs <- mr$sims.list$lik_yobs
  lik_yhat <- mr$sims.list$lik_yhat
  lik_sucobs <- mr$sims.list$lik_sucobs
  lik_suchat <- mr$sims.list$lik_suchat
  #assign temporary objects
  ll_yobs_k <- ll_yhat_k <- ll_sucobs_k <- ll_suchat_k <- matrix( nrow = S, ncol = K )
  ll_yobs <- ll_yhat <- ll_sucobs <- ll_suchat <- Dev_obs <- Dev_hat <- rep( NA, S )

  for( s in 1:S ){ #loop over MCMC iterations
    for( k in 1:K ){ #loop over years
      # sum log likelihoods only for sites sampled each year inside VNP:
      ll_yobs_k[ s, k ] <- sum( log( lik_yobs[ s, 1:M, k ] ) *
                               XM$inVNP[ 1:M ] * surv_ind[ 1:M, k ]  )
      ll_yhat_k[ s, k ] <- sum( log( lik_yhat[ s, 1:M, k ] ) *
                               XM$inVNP[ 1:M ] * surv_ind[ 1:M, k ] )
      ll_sucobs_k[ s, k ] <- sum( log( lik_sucobs[ s, 1:M, k ] ) *
                                 XM$inVNP[ 1:M ] * surv_ind[ 1:M, k ] )
      ll_suchat_k[ s, k ] <- sum( log( lik_suchat[ s, 1:M, k ] ) *
                                 XM$inVNP[ 1:M ] * surv_ind[ 1:M, k ] )
    } #close K loop
    #sum log likelihoods across years
    ll_yobs[ s ] <- sum( ll_yobs_k[ s, 1:K ] ) #for observed occupancy
    ll_yhat[ s ] <- sum( ll_yhat_k[ s, 1:K ] ) #for estimated occupancy
    ll_sucobs[ s ] <- sum( ll_sucobs_k[ s, 1:K ] ) #for observed success
    ll_suchat[ s ] <- sum( ll_suchat_k[ s, 1:K ] ) #for estimated success
    #print( ll_yobs[ s ] )
    # Calculate overall deviances
    Dev_obs[ s ] <- -2 * ( ll_yobs[ s ] + ll_sucobs[ s ] )
    Dev_hat[ s ] <- -2 * ( ll_yhat[ s ] + ll_suchat[ s ] )
  } # close s loop
  
  return( data.frame( Dev_obs, Dev_hat, ll_yobs, ll_sucobs, ll_yhat, ll_suchat ) )
  
} # close function

# Run function to calculate deviances
ModDevs <- CalcDevs(  )
#head( ModDevs ); dim( ModDevs )
# Bayesian p-value define as mean number of times that Deviance of observed data was #
# greater than the deviance of predicted model observations:
Baypvalue <- mean( ModDevs$Dev_obs > ModDevs$Dev_hat )
plot( ModDevs$Dev_obs, ModDevs$Dev_hat, main = paste( "Bayesian P value = ", 
        round( Baypvalue, 3 ) ), xlab = "Deviance of observed data", 
      ylab = "Deviance of predicted data", tcl = 0.2, cex.lab = 1.5, bty = "l"  ) 
abline( 0, 1, col = 'black', lwd = 2 )


######################    END OF SCRIPT             ###########################################
##################################################################################################