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

small.world <- ( small.world
    ## sort so that i == ID, so that epi links to ID work
    %>% arrange( IID )
    ## add column of known sources for W-T algorithm.
    ## NA = not known
    ## -1 = primary case
    ## other number = known source
    %>% mutate( v = ifelse(
        Transmission == 'Disney primary case',
        -1,
	ifelse(
	    Transmission == 'Unknown' & ROD %in% c(0,4,7),
	    -1,
            EpilinkID
	)
    ) )
)
#print( small.world %>% filter( v == -1 ) %>% nrow )
## 45

plot.transmission.graph(
	small.world,
	transmute( small.world, s=EpilinkID, t=IID, w=1 ),
	'smallworld-transmission.pdf'
)

## which generation time function to use
#w.fn <- w.piecewise.constant
w.fn <- w.gamma.sfa #make.w.gamma( 9, 2.25 )

## first run it with no smarts
## it overdoes initial R because trying to infect all the primary cases
small.world %<>% mutate( R.naive = wallinga.teunis( ROD, NA*ROD, w.fn, NULL ) )
p <- ( ggplot( small.world, aes( x=ROD, y=R.naive ) )
	+ geom_point()
	+ geom_line( size=0 )
	+ geom_line( aes( x=ROD ), y=1 )
	+ labs( x='Day', y='R' )
)
png( 'smallworld-naive.png' )
print(p)
dev.off()

## now with the primary cases left out of infection process
small.world %<>% mutate( R.simple = wallinga.teunis( ROD, ifelse( v == -1, -1, NA ), w.fn, NULL ) )
p <- ( ggplot( small.world, aes( x=ROD, y=R.simple ) )
	+ geom_point()
	+ geom_line( size=0 )
	+ geom_line( aes( x=ROD ), y=1 )
	+ labs( x='Day', y='R' )
)
png( 'smallworld-simple.png' )
print(p)
dev.off()

## now with primary cases and epi links
small.world %<>% mutate( R.epi = wallinga.teunis( ROD, v, w.fn, NULL ) )
p <- ( ggplot( small.world, aes( x=ROD, y=R.epi ) )
	+ geom_point()
	+ geom_line( size=0 )
	#+ geom_smooth()
	+ geom_segment( data=avg.by.week( select( small.world, ROD, R=R.epi ) ), aes( x = week_start, xend=week_start + 6, y=R, yend=R ), color='red', size=3 )
	+ geom_line( aes( x=ROD ), y=1 )
	+ labs( x='Day', y='R' )
)
png( 'smallworld-simple-epi.png' )
print(p)
dev.off()

## transmission graph with these imputed link likelihoods
pij <- ( small.world
	%>% { p.matrix( .$ROD, .$v, w.matrix( .$ROD, .$v, w.fn, NULL ) ) }
	%>% tidy.matrix()
)
plot.transmission.graph( small.world,
	mutate(pij, s=j, t=i, w=p),
	'smallworld-transmission-epi.pdf'
)

## R.epi vs age
p <- ( ggplot( small.world, aes( x=Age, y=R.epi ) )
	+ geom_point()
	+ scale_y_log10()
)
png( 'smallworld-by-age.png' )
print(p)
dev.off()

## estimate generation membership of each case
gens <- ( small.world
	%>% { estimate.generations( .$ROD, .$v, w.fn, NULL ) }
	%>% mutate( generation = fct_other(factor(generation), keep=c(0,1,2,3,4), other_level='5+') )
)

## plot generations by time
g <- ( gens
	%>% group_by( t, generation )
	%>% summarize( p = sum(p) )
	%>% ungroup
)
p <- ( ggplot( g, aes( x=t, y=p ) )
	+ geom_col( aes( fill=generation ) )
)
png( 'smallworld-generations-epi.png' )
print(p)
dev.off()

## R by generation
R.gen <- ( inner_join( gens, small.world, by=c(i='IID') )
	#%T>% { print(head(.)) }
	%>% group_by( generation )
	%>% summarize( R = sum( R.epi * p ) / sum( p ) )
)
p <- ( ggplot( R.gen )
	+ geom_col( aes( x=generation, y=R, fill=generation ) )
	+ guides( fill='none' )
)
png( 'smallworld-R-by-generation-epi.png' )
print(p)
dev.off()

## transmission per vaccine status
small.world %<>% mutate( vaccinated = fct_explicit_na(factor(VaccineDoses > 0), na_level='-') )

v <- ( estimate.R.by.status( small.world$ROD, small.world$v, w.fn, NULL, small.world$vaccinated )
	%>% filter( type %in% c('FALSE -> FALSE','FALSE -> TRUE','TRUE -> FALSE','TRUE -> TRUE') )
)

(v
	%>% group_by( t.j, type )
	%>% summarize( R = sum(p) )
	%>% ggplot()
	#%>% add( geom_col( aes( x=t.j, y=R, fill=type ) ) )
	%>% add( geom_line( aes( x=t.j, y=R, color=type ) ) )
	%>% { png( 'smallworld-R-by-status.png' ); print(.); dev.off() }
)

print(head(v))
print(head(gens))

## transmission by generation and status
( inner_join( v, gens, by=c(j='i') )
	%>% mutate( p.ji=p.x, p.gen = p.y, vaccination=type )
	%T>% { print(head(.)) }
	%>% group_by( vaccination, generation )
	%>% summarize( R=sum(p.ji*p.gen) / sum(p.gen) )
	%T>% { print(head(.)) }
	%>% ggplot()
	%>% add( geom_col( aes( x=generation, y=R, fill=vaccination ), position='dodge' ) )
	%>% { png( 'smallworld-R-by-status-by-generation.png' ); print(.); dev.off() }
)

## now try with penalty for cross-county transmission
penalty <- function(i,j) {
	ifelse( small.world$LHJ[i] == small.world$LHJ[j], 1, 0.01 )
}
#print( select( ss, IID, ROD, Transmission, EpilinkID, v ) )
small.world %<>% mutate( R.penalty = wallinga.teunis( ROD, v, w.fn, penalty ) )
p <- ( ggplot( small.world, aes( x=ROD, y=R.penalty ) )
	+ geom_point()
	+ geom_line( size=0 )
	#+ geom_smooth()
	+ geom_segment( data=avg.by.week( select( small.world, ROD, R=R.penalty ) ), aes( x = week_start, xend=week_start + 6, y=R, yend=R ), color='red', size=3 )
	+ geom_line( aes( x=ROD ), y=1 )
	+ labs( x='Day', y='R' )
)
png( 'smallworld-simple-penalty.png' )
print(p)
dev.off()

write.csv( small.world, 'smallworld-estimates.csv', row.names=F, quote=F )
