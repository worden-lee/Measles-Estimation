export WW_THIS_PROJECT_SOURCE_FILES=Makefile wallinga-teunis.r
export RESOURCES=/usr/local/src/workingwiki/ProjectEngine/resources
export WW_OUTPUTFORMAT=tex
-include $(RESOURCES)/makefile-before
export DEFAULT_LATEX=pdflatex
-include $(RESOURCES)/makefile-after

%.r.out : %.Rout
	mv $< $@

%.R : %.r
	cp -f $< $@
