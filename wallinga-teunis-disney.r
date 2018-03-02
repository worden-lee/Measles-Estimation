#!/usr/bin/R
library(tidyverse)
library(magrittr)

source( 'wallinga-teunis-functions.r' )

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

## get the disney epidemic curve
small.world <- read.csv( 'smallword_deid.csv', row.names=NULL )

## first run it with no smarts
## it overdoes initial R because trying to infect all the primary cases
Rj <- wallinga.teunis( small.world$ROD, NA*small.world$ROD, w.gamma, NULL )
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
Rj <- wallinga.teunis( small.world$ROD, ifelse( small.world$Transmission == 'Disney primary case', -1, NA ), w.gamma, NULL )
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
#print( select( ss, IID, ROD, Transmission, EpilinkID, v ) )
ss %<>% mutate( R = wallinga.teunis( ROD, v, w.gamma, NULL ) )
print( head(ss) )
p <- ( ggplot( ss, aes( x=ROD, y=R ) )
	+ geom_point()
	+ geom_line( width=0 )
	#+ geom_smooth()
	+ geom_segment( data=avg.by.week( ss ), aes( x = week_start, xend=week_start + 6, y=R, yend=R ), color='red', width=3 )
	+ geom_line( aes( x=ROD ), y=1 )
	+ labs( x='Day', y='R' )
)
png( 'smallworld-simple-epi.png' )
print(p)
dev.off()

p <- ( ggplot( ss, aes( x=Age, y=R ) )
	+ geom_point()
	+ scale_y_log10()
)
png( 'smallworld-by-age.png' )
print(p)
dev.off()

## now try with penalty for cross-county transmission
penalty <- function(i,j) {
	ifelse( ss$LHJ[i] == ss$LHJ[j], 1, 0.01 )
}
print( select( ss, IID, ROD, Transmission, EpilinkID, v ) )
ss %<>% mutate( R = wallinga.teunis( ROD, v, w.gamma, penalty ) )
p <- ( ggplot( ss, aes( x=ROD, y=R ) )
	+ geom_point()
	+ geom_line( width=0 )
	#+ geom_smooth()
	+ geom_segment( data=avg.by.week( ss ), aes( x = week_start, xend=week_start + 6, y=R, yend=R ), color='red', width=3 )
	+ geom_line( aes( x=ROD ), y=1 )
	+ labs( x='Day', y='R' )
)
png( 'smallworld-simple-penalty.png' )
print(p)
dev.off()

