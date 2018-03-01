#!/usr/bin/R
library(tidyverse)

## replace ggplot's characteristic gray background with a more
## standard white one
theme_set(theme_bw(12)
   + theme(#legend.title=element_blank(), legend.position=c(0.88,0.22),
           legend.key=element_rect(linetype='blank',fill='white'),
           axis.line=element_line(linetype=1),
           #axis.line=element_blank(),
           #axis.ticks=element_blank(),
           #strip.background=element_rect(linetype='blank',fill=NULL),
           strip.background=element_blank(),
           panel.border=element_rect(linetype='blank',fill=NULL),
           #panel.grid.major=element_blank(),
           panel.grid.minor=element_blank()
   )
)

## given a w function
w.piecewise.constant <- function( tau, theta ) {
	return( case_when(
		tau < 4 ~ 0,
		tau < 8 ~ 1 / ( (22-8) + (26-4) ),
		tau < 22 ~ 2 / ( (22-8) + (26-4) ),
		tau < 26 ~ 1 / ( (22-8) + (26-4) ),
		T ~ 0
	) )
}

w.gamma <- function( tau, theta ) {
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

## and an epidemic curve
## i.e. list of times of symptom onset t
#t <- c( 0, 8, 11, 22, 29, 37 )

## reduced Wallinga-Teunis (p.511 of W and T 2004)
## t: list of onset dates, indexed by index "i" of individual cases
## v: list of transmission sources, indexed by i; NA if unknown; -1 if primary case
## w: density function for serial interval
## returns: list of R values indexed by i
reduced.wallinga.teunis <- function( t, v, w, penalty ) {
	cat( 'here are ROD and epi links\n' )
	print( t )
	print( v )

	if ( is.null( penalty ) ) {
		penalty <- function(i,j) 1
	}

	## record all w_ij = w( ti - tj )
	## this is likelihood for infection of i given source j
	## if source of i is not unknown, skip
	wij <- outer( seq_along(t), seq_along(t), FUN=function( i, j ) ifelse( is.na(v[i]), w( t[i] - t[j], c() ) * penalty(i,j), 0 ) )
	#cat( 'w_ij\n' )
	#print( wij )

	## construct all p_ij = w_ij / sum(k) w_ik
	## this is relative likelihood for j being source of i
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

	## construct R_j = sum(ti>tj) p_ij
	## an expected number of cases caused by j
	Rj = sapply( seq_along(t), FUN=function(j) sum( sapply( seq_along(t), FUN=function(i) pij[i,j] ) ) )
	#cat( 'R_j\n' )
	#print( Rj )
	Rj
}

## try the full Wallinga Teunis, with enumeration of all the networks,
## see how far we can get
big.wallinga.teunis <- function( t, v, w ) {
	## there are two formulations, I want to look at both

	## first an integration over all network structures
	## maybe do it monte carlo
	tf <- data.frame( t = t ) # df form of onset date sequence
	tq <- filter( tf, t > 0 ) # the non-primary cases
	for (i in 1:1000) {
	    ## consider all possible networks with equal weight
	    ## a network is: for each case i with unknown source,
	    ## assign a source j /= i.
	    ## networks that have likelihood zero might as well
	    ## be excluded, including any with cycles
	    ## or with sources whose onset date is later than target's.
	    ## so, for example, we can sample from networks by
	    ## choosing each j uniformly from all the earlier cases.

	    ## given a network, 
	}
}

#Rj = reduced.wallinga.teunis(t, w.gamma)

## get the disney epidemic curve
small.world <- read.csv( 'smallword_deid.csv', row.names=NULL )

## first run it with no smarts
## it overdoes initial R because trying to infect all the primary cases
Rj <- reduced.wallinga.teunis( small.world$ROD, NA*small.world$ROD, w.gamma, NULL )
p <- ( ggplot( data.frame( t=small.world$ROD, R=Rj ) )
	+ geom_point( aes( x=t, y=Rj ) )
	+ geom_line( aes( x=t, y=Rj ), width=0 )
	+ geom_line( aes( x=t ), y=1 )
	+ labs( x='Day', y='R' )
)
png( 'smallworld-naive.png' )
print(p)
dev.off()

## now with the primary cases left out of infection process
Rj <- reduced.wallinga.teunis( small.world$ROD, ifelse( small.world$Transmission == 'Disney primary case', -1, NA ), w.gamma, NULL )
p <- ( ggplot( data.frame( t=small.world$ROD, R=Rj ) )
	+ geom_point( aes( x=t, y=Rj ) )
	+ geom_line( aes( x=t, y=Rj ), width=0 )
	+ geom_line( aes( x=t ), y=1 )
	+ labs( x='Day', y='R' )
)
png( 'smallworld-simple.png' )
print(p)
dev.off()

## now with primary cases and epi links
## sort so that i == ID, so that epi links to ID work
ss <- ( small.world
    %>% arrange( IID )
    %>% mutate( v = ifelse(
        Transmission == 'Disney primary case',
        -1,
        EpilinkID
    ) )
)
print( select( ss, IID, ROD, Transmission, EpilinkID, v ) )
Rj <- reduced.wallinga.teunis( ss$ROD, ss$v, w.gamma, NULL )
p <- ( ggplot( data.frame( t=ss$ROD, R=Rj ) )
	+ geom_point( aes( x=t, y=Rj ) )
	+ geom_line( aes( x=t, y=Rj ), width=0 )
	+ geom_line( aes( x=t ), y=1 )
	+ labs( x='Day', y='R' )
)
png( 'smallworld-simple-epi.png' )
print(p)
dev.off()

## now try with penalty for cross-county transmission
penalty <- function(i,j) {
	ifelse( ss$LHJ[i] == ss$LHJ[j], 1, 0.01 )
}
print( select( ss, IID, ROD, Transmission, EpilinkID, v ) )
Rj <- reduced.wallinga.teunis( ss$ROD, ss$v, w.gamma, penalty )
p <- ( ggplot( data.frame( t=ss$ROD, R=Rj ) )
	+ geom_point( aes( x=t, y=Rj ) )
	+ geom_line( aes( x=t, y=Rj ), width=0 )
	+ geom_line( aes( x=t ), y=1 )
	+ labs( x='Day', y='R' )
)
png( 'smallworld-simple-penalty.png' )
print(p)
dev.off()

