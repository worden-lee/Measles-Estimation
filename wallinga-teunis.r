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
reduced.wallinga.teunis <- function( t, w ) {
	## record all w_ij = w( ti - tj )
	wij <- outer( t, t, FUN=function( ti, tj ) w( ti - tj, c() ) )
	cat( 'w_ij\n' )
	print( wij )

	## construct all p_ij = w_ij / sum(k) w_ik
	pij <- outer( seq_along(t), seq_along(t), FUN=function(i,j) sapply( seq_along(i), FUN=function(k) (wij[i[k],j[k]] / sum(wij[i[k],])) ) )
	cat( 'p_ij\n' )
	print( pij )

	## construct R_j = sum(ti>tj) p_ij
	## note, ti > tj may not be sufficient to avoid NaNs if multiple
	## importations
	Rj = sapply( seq_along(t), FUN=function(j) sum( sapply( seq_along(t), FUN=function(i) ifelse( t[i] > t[j], pij[i,j], 0 ) ) ) )
	cat( 'R_j\n' )
	print( Rj )
	Rj
}

## try the full Wallinga Teunis, with enumeration of all the networks,
## see how far we can get
big.wallinga.teunis <- function( t ) {

}
#Rj = reduced.wallinga.teunis(t, w.gamma)

## get the disney epidemic curve
small.world <- read.csv( 'smallword_deid.csv', row.names=NULL )
Rj <- reduced.wallinga.teunis( small.world$ROD, w.gamma )

p <- ( ggplot( data.frame( t=small.world$ROD, R=Rj ) )
	+ geom_point( aes( x=t, y=Rj ) )
	+ geom_line( aes( x=t ), y=1 )
	+ labs( x='Day', y='R' )
)
png( 'smallworld.png' )
print(p)
dev.off()
