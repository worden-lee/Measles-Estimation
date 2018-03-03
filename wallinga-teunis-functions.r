#!/usr/bin/R
## functions for wallinga-teunis estimation of infection numbers
library(tidyverse)
library(magrittr)

## define some "w" functions
## this is the probability density function of serial interval tau.

## an oversimple w
w.piecewise.constant <- function( tau, theta ) {
	return( case_when(
		tau < 4 ~ 0,
		tau < 8 ~ 1 / ( (22-8) + (26-4) ),
		tau < 22 ~ 2 / ( (22-8) + (26-4) ),
		tau < 26 ~ 1 / ( (22-8) + (26-4) ),
		T ~ 0
	) )
}

## the gamma distribution used in Paul's paper
## this is a little quick compared to the 8-21 days we usually use
w.gamma.pg <- function( tau, theta ) {
	m = 11.1  # mean = a*s
	sd = 2.47 # var = a*s^2
	s = sd * sd / m
	a = m / s
	return ( ifelse( tau <= 0,
		0,
		( pgamma( tau + 0.5, shape=a, scale=s ) -
		  pgamma( tau - 0.5, shape=a, scale=s )
		)
	) )
}

## the gamma distribution used in Paul's paper
## this is a little quick compared to the 8-21 days we usually use
w.gamma.sfa <- function( tau, theta ) {
	m = 14.5  # mean = a*s
	sd = 3.25 # var = a*s^2
	s = sd * sd / m
	a = m / s
	return ( ifelse( tau <= 0,
		0,
		( pgamma( tau + 0.5, shape=a, scale=s ) -
		  pgamma( tau - 0.5, shape=a, scale=s )
		)
	) )
}

## w matrix used in Wallinga-Teunis
## t: list of onset dates
## v: list of transmission sources (see below)
## w: serial interval density function
## returns: matrix of probability of transmission from j to i for each i,j
w.matrix <- function( t, v, w, penalty ) {
	## this is likelihood for infection of i given source j

	if ( is.null( penalty ) ) {
		penalty <- function(i,j) 1
	}

	## if source of i is not unknown, skip
	wij <- outer( seq_along(t), seq_along(t), FUN=function( i, j ) ifelse( is.na(v[i]), w( t[i] - t[j], c() ) * penalty(i,j), 0 ) )
	#cat( 'w_ij\n' )
	#print( wij )
	wij
}

## p matrix used in Wallinga-Teunis
## returns: matrix of relative likelihoods of each j as source for i
p.matrix <- function( t, v, wij ) {
	## p is 1 where j is known, 0 if i a primary case
	pij <- outer( seq_along(t), seq_along(t), FUN=function(i,j)
	    sapply( seq_along(i), FUN=function(k)
	        ifelse( is.na(v[i[k]]),
		    ifelse( wij[i[k],j[k]] == 0,
			0,
			(wij[i[k],j[k]] / sum(wij[i[k],]))
		    ),
		    ifelse( v[i[k]] == j[k], 1, 0 )
		)
	    )
	)
	#cat( 'p_ij\n' )
	#print( pij )
	pij
}

## Wallinga-Teunis formula (p.511 of W and T 2004, and appendix)
## t: list of onset dates, indexed by index "i" of individual cases
## v: list of transmission sources, indexed by i; NA if unknown; -1 if primary case
## w: density function for serial interval
## penalty: multiplier for likelihood depending on other details of the two cases
## returns: list of R values indexed by i
wallinga.teunis <- function( t, v, w, penalty ) {
	cat( 'here are ROD and epi links\n' )
	print( t )
	print( v )

	## record all w_ij = w( ti - tj )
	wij <- w.matrix( t, v, w, penalty )

	## construct all p_ij = w_ij / sum(k) w_ik
	## this is relative likelihood for j being source of i
	pij <- p.matrix( t, v, wij )

	## construct R_j = sum(ti>tj) p_ij
	## an expected number of cases caused by j
	Rj <- sapply( seq_along(t), FUN=function(j) sum( sapply( seq_along(t), FUN=function(i) pij[i,j] ) ) )
	#cat( 'R_j\n' )
	#print( Rj )
	Rj
}

## estimate likely generation membership of cases using W-T-like summation
estimate.generations <- function( t, v, w, penalty ) {
	## wij = prob( j -> i )
	wij <- w.matrix( t, v, w, penalty )

	## pij = L( j -> i | i )
	pij <- p.matrix( t, v, wij )

	g <- data.frame( i=c(), t=c(), generation=c(), p=c() )
	## p( i is in generation 0 ) = delta( i an index case )
	## p_0[i] = delta( v(i), -1 )
	p.i <- ifelse( !is.na(v) & v == -1, 1, 0 )
	gen <- 0
	## p( i is in generation n ) = sum_j( p( j in n ) L( j -> i ) )
	## p_n[i] = sum_j p_{n-1}[j] pij[i,j]
	## p_n = pij p_{n-1}
	while( sum( p.i ) > 0 ) {
	    g %<>% rbind( data.frame(
		i = seq_along(t),
		t = t,
		generation = gen,
		p = p.i
	    ) )
	    p.i <- sapply( seq_along(t), FUN=function(i) {
		sum( sapply( seq_along(t), FUN=function(j) pij[i,j]*p.i[j] ) )
	    } )
	    gen <- gen + 1
	}
	g
}

avg.by.week <- function( df ) {
    ( df
	%>% mutate( week_start = as.integer( ROD / 7 ) * 7 )
	%>% group_by( week_start )
	%>% summarize( R = mean(R) )
    )
}

