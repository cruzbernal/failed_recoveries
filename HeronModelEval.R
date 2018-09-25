##############################################################################
###                  Code created by Jen Cruz                         ########
#### This script contains methods for calculating model evaluation stats  ####
####            based on likelihood of observed and predicted data      ######
####                    given model p(data|model).                        ####
###   Model estimates heronry occurrence and mean size for great blue    #####
###       herons at Voyageurs National Park, Minnesota, USA (1986-2012). #####
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
#start with results from heron model:
load( "D:/RProjects/herons/FullModelResultsV2.RData" )

##################################################################################################
## Calculate summary stats from likelihood of data conditional on model: p(y|theta)  ###
##################################################################################################
# Define model we want to extract output from 
mr <- cm10
summary( mr )
#Total number of iterations ran:
S <- length( mr$sims.list$int.omg )#

# Define function to calculate model evaluation parameters 
EvaluateModel <- function( lik_yobs, llik_nobs ){
  # create temporary objects
  CPO_yobs_ik <- CPO_nobs_ik <- ll_yobs_ik <- ll_nobs_ik <- var_yobs_ik <- var_nobs_ik <- matrix( NA, ncol = K, nrow = M )
  CPO_yobs_k <- CPO_nobs_k <- ll_yobs_k <- ll_nobs_k <- var_yobs_k <- var_nobs_k <- rep( NA, K )
  # avoid numerical overflow of exp(>710) as Inf
  llik_nobs[ llik_nobs > 700 ] <- 700
#  print( sum( as.numeric( which( llik_nobs > 709 ) ) ) )
  for( k in 1:K ){ # loop over years
    for( i in 1:M ){ #loop over sites
      if ( XM$inVNP[ i ] * surv_ind[ i, k ] == 1 ){
        # Calculate the approximation of the Conditional Predictive Ordinate (CPO) across MCMC iterations 
        #for each site, each year:
        CPO_yobs_ik[ i, k ] <- S / sum( 1 / lik_yobs[ 1:S, i, k ] ) 
        # replace division by 0 with 0 outcome
        CPO_nobs_ik[ i, k ] <- ifelse( sum( 1 / exp( llik_nobs[ 1:S, i, k ] ) )  == 0, 0, 
                                      S / sum( 1 / exp( llik_nobs[ 1:S, i, k ] ) ) )
        CPO_nobs_ik[ i, k ] <- log( CPO_nobs_ik[ i,k ] + 1e-10 )
        #Calculate components of the Watanabe-Akaike Information Criterion (WAIC):
        # Based on Appendix S1 equations from Broms et al. 2016.
        #Start by summing likelihood of data across MCMC iterations, dividing by MCMC steps, then logging it
        ll_yobs_ik[ i, k ] <- log( sum( lik_yobs[ 1:S, i, k ] ) / S )
        # Now for heronry size:
        ll_nobs_ik[ i, k ] <- log( 1e-10 + sum( exp( llik_nobs[ 1:S, i, k ] ) ) / S )
        #calculate variance between log likelihood computed at each MCMC step for each site, each
        #year and the mean log likelihood across MCMC steps:
        var_yobs_ik[ i, k ] <-  ( 1 / ( S - 1 ) ) *
                      sum( ( log( lik_yobs[ 1:S, i, k ]  +  1e-10 ) - 
                               log( mr$mean$lik_yobs[ i, k ]  +  1e-10 ) )^2 )
        var_nobs_ik[ i, k ] <-  ( 1 / ( S - 1 ) ) *
                    sum( ( llik_nobs[ 1:S, i, k ] - mr$mean$llik_nobs[ i, k ] )^2 )
      } #if loop
    } #M loop
    # Calculate the summary statistic of individual CPO values by summing over sites that were sampled each year inside VNP:
    #adding a very small value to avoid numeric overload
    CPO_yobs_k[ k ] <- - sum( log( CPO_yobs_ik[ 1:M, k ] +  1e-10 ), na.rm = TRUE )  
    CPO_nobs_k[ k ] <- - sum( CPO_nobs_ik[ 1:M, k ], na.rm = TRUE ) 
    # sum log likelihood only across sites that were sampled each year inside VNP:
    ll_yobs_k[ k ] <- sum( ll_yobs_ik[ 1:M, k ], na.rm = TRUE )
    ll_nobs_k[ k ] <- sum( ll_nobs_ik[ 1:M, k ], na.rm = TRUE )
    # sum variances only across sites sampled each year inside VNP:
    var_yobs_k[ k ] <- sum( var_yobs_ik[ 1:M, k ] , na.rm = TRUE)
    var_nobs_k[ k ] <- sum( var_nobs_ik[ 1:M, k ], na.rm = TRUE )
  } #k loop

  # Calculate measures of predictive performance (CPO) for each model aspect:
  CPOyobs <- sum( CPO_yobs_k[ !is.na( CPO_yobs_k ) ] ) 
  CPOnobs <- sum( CPO_nobs_k[ !is.na( CPO_nobs_k ) ] )
  
  # Calculate elppd of WAIC:
  elppd <- sum( ll_yobs_k[ 1:K ], na.rm = TRUE ) + sum( ll_nobs_k[ 1:K ], na.rm = TRUE )
  # Calculate model complexity component of WAIC following Gelman et al. 2014 pWAIC2:
  pWAIC2 <- sum( var_yobs_k[ 1:K ], na.rm = TRUE ) + sum( var_nobs_k[ 1:K ], na.rm = TRUE )
  # Calculate WAIC
  WAIC <- -2 * elppd + 2 * pWAIC2 
  # Function output:
  return( data.frame( CPOyobs, CPOnobs, elppd, pWAIC2, WAIC ) )
}
# Run function and store:
EvalRes <- EvaluateModel( mr$sims.list$lik_yobs, mr$sims.list$llik_nobs )
EvalRes

# Define function to calculate Model deviances to be used in Bayesian p value calculations
CalcDevs <- function( ){
  # Define initial values:
  lik_yobs <- mr$sims.list$lik_yobs
  lik_yhat <- mr$sims.list$lik_yhat
  llik_nobs <- mr$sims.list$llik_nobs
  llik_nhat <- mr$sims.list$llik_nhat
  #assign temporary objects
  ll_yobs_k <- ll_yhat_k <- ll_nobs_k <- ll_nhat_k <- matrix( nrow = S, ncol = K )
  ll_yobs <- ll_yhat <- ll_nobs <- ll_nhat <- Dev_obs <- Dev_hat <- rep( NA, S )

  for( s in 1:S ){ #loop over MCMC iterations
    for( k in 1:K ){ #loop over years
      # sum log likelihoods only for sites sampled each year inside VNP:
      ll_yobs_k[ s, k ] <- sum( log( lik_yobs[ s, 1:M, k ] ) *
                               XM$inVNP[ 1:M ] * surv_ind[ 1:M, k ]  )
      ll_yhat_k[ s, k ] <- sum( log( lik_yhat[ s, 1:M, k ] ) *
                               XM$inVNP[ 1:M ] * surv_ind[ 1:M, k ] )
      ll_nobs_k[ s, k ] <- sum( llik_nobs[ s, 1:M, k ] *
                                 XM$inVNP[ 1:M ] * surv_ind[ 1:M, k ] )
      ll_nhat_k[ s, k ] <- sum( llik_nhat[ s, 1:M, k ] *
                                 XM$inVNP[ 1:M ] * surv_ind[ 1:M, k ] )
    } #close K loop
    #print( ll_yobs_k[ s,  ] )
    #sum log likelihoods across years
    ll_yobs[ s ] <- sum( ll_yobs_k[ s, 1:K ] ) #for observed 
    ll_yhat[ s ] <- sum( ll_yhat_k[ s, 1:K ] ) #for estimated 
    ll_nobs[ s ] <- sum( ll_nobs_k[ s, 1:K ] ) #from observed 
    ll_nhat[ s ] <- sum( ll_nhat_k[ s, 1:K ] ) #from estimated 
    # Calculate overall deviances
    Dev_obs[ s ] <- -2 * ( ll_yobs[ s ] + ll_nobs[ s ] )
    Dev_hat[ s ] <- -2 * ( ll_yhat[ s ] + ll_nhat[ s ] )
  } # close s loop
  
  return( data.frame( Dev_obs, Dev_hat, ll_yobs, ll_nobs, ll_yhat, ll_nhat ) )
  
} # close function

# Run function to calculate deviances
ModDevs <- CalcDevs(  )
#head( ModDevs ); dim( ModDevs )
# Bayesian p-value defined as mean number of times that Deviance of observed data was #
# greater than the deviance of predicted model observations:
Baypvalue <- mean( ModDevs$Dev_obs > ModDevs$Dev_hat )
plot( ModDevs$Dev_obs, ModDevs$Dev_hat, main = paste( "Bayesian P value = ", 
        round( Baypvalue, 3 ) ), xlab = "Deviance of observed data", 
      ylab = "Deviance of predicted data", tcl = 0.2, cex.lab = 1.5, bty = "l"  ) 
abline( 0, 1, col = 'black', lwd = 2 )

######################    END OF SCRIPT             ###########################################
##################################################################################################