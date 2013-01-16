#
# Make MAPPINGS III
#
#
#
#    v1.0.0q
#-------------------------------
#--- Output (executable) name --
#-------------------------------
#
OUTNAME = map3r
#
#-------------------------------
#---------- Directories --------
#-------------------------------
#
EXEDIR = ../lab/
INCDIR = includes/
CODDIR = workcode/
BACKDIR = mastercode/
#
#-------------------------------
#---------- Compilers ----------
#-------------------------------
#
#  comment out and in as needed 
#   for different machines
#
#  MacOSX f2c
#
#FC     = fc
#CC     = mwccppc
#LDR    = mwldppc -exe
#OPTS   = -c -I${INCDIR} -f 
#LIB    = -L/usr/local/lib/ -lf2c 
#CPP    = /usr/bin/cpp -P -Iincludes/ -I.
#
#
# MacOSX Absoft FORTRAN
#        Note that Shock4 must be commented out
#
#FC     = f77
#LDR    =  -s -O3 -N109  -N16 -W -N34
#OPTS   = -c ${LDR} -I${INCDIR}
#CC     = cc
#LIB    = -lU77  -lf90math -lfio -lf77math
#CPP    = cpp -P -xassembler-with-cpp -Iincludes/ -I.
#
#
#  Linux IFC FORTRAN - replace ifc with f95i for later version
#
#FC     = ifc
#CC     = mwccppc
#LDR    = ifc -Vaxlib
#OPTS   = -I${INCDIR} -B18 -O3 -c -N109 -W -N34
#OPTS	= -i4 -r8 -132 -w95 -Vaxlib -I${INCDIR} -c -e -o $*.o -I${INCDIR}
#LIB    = 
#CPP    = cpp -P -Iincludes/ -I.
#
#
#     SunOS - xtarget can be dropped or changed depending on machine
#
#FC     = f77
#OPTS   = -c -e -o $*.o -I${INCDIR} -fast -xO5 -xtarget=ultra3
#OPTS   = -c -e -o $*.o -I${INCDIR} -g
#CC     = cc
#LDR    =  -fast -xO5 -xtarget=ultra3
#LDR    =  -g
#LIB    = 
#CPP    = /usr/bin/cpp -P -Iincludes/ -I.
#
#
#
#     g77 universal compilier
#	Note that optimization must be dropped for Mac OSX
#
FC     = g77
#FC	= gfortran -f77
#OPTS   = -c  -I${INCDIR} -ffixed-line-length-132
OPTS   = -c -ffixed-line-length-132  -I${INCDIR} -g -ffortran-bounds-check -W  -Wall -pedantic -Wno-globals -Wimplicit
CC     = gcc
#LDR    =   -ffixed-line-length-132
LDR    =  -g -ffortran-bounds-check -Wsurprising
LIB    = 
CPP    = /usr/bin/cpp -P -Iincludes/ -I.
#
#
#
.SUFFIXES: .f .o
#
#-------------------------------
#-- include/header files -------
#-------------------------------
#
INCS =  ${INCDIR}cblocks.inc \
	${INCDIR}const.inc
#
#-------------------------------
#-- object files ---------------
#-------------------------------
#
OBJ = ${CODDIR}mappings.o \
	${CODDIR}abecha.o \
	${CODDIR}absdis.o \
	${CODDIR}allrates.o \
	${CODDIR}averinto.o \
	${CODDIR}avrdata.o \
	${CODDIR}chacha.o \
	${CODDIR}charex.o \
	${CODDIR}cheat.o \
	${CODDIR}collion.o \
	${CODDIR}coloss.o \
	${CODDIR}compph3.o \
	${CODDIR}compph4.o \
	${CODDIR}compph5.o \
	${CODDIR}compsh2.o \
	${CODDIR}compsh3.o \
	${CODDIR}compton.o \
	${CODDIR}cool.o \
	${CODDIR}coolc.o \
	${CODDIR}copinto.o \
	${CODDIR}cosmic.o \
	${CODDIR}crosssections.o \
	${CODDIR}depcha.o \
	${CODDIR}difpop.o \
	${CODDIR}dumpp4.o \
	${CODDIR}dumps4.o \
	${CODDIR}dusttemp.o \
	${CODDIR}equion.o \
	${CODDIR}evoltem.o \
	${CODDIR}fine3.o \
	${CODDIR}fivelevel.o \
	${CODDIR}findtde.o \
	${CODDIR}fpol.o \
	${CODDIR}freebound.o \
	${CODDIR}freefree.o \
	${CODDIR}functions.o \
	${CODDIR}grainpar.o \
	${CODDIR}helioi.o \
	${CODDIR}helif.o \
	${CODDIR}hgrains.o \
	${CODDIR}hpahs.o \
	${CODDIR}hydrec.o \
	${CODDIR}hydro.o \
	${CODDIR}hydro2p.o \
	${CODDIR}inter.o \
	${CODDIR}inter2.o \
	${CODDIR}interpol.o \
	${CODDIR}intvec.o \
	${CODDIR}invrho.o \
	${CODDIR}iobal.o \
	${CODDIR}iohyd.o \
	${CODDIR}ionab.o \
	${CODDIR}ionsec.o \
	${CODDIR}ironii.o \
	${CODDIR}liftoff.o \
	${CODDIR}localem.o \
	${CODDIR}lookup.o \
	${CODDIR}mapinit.o \
	${CODDIR}mdiag.o \
	${CODDIR}ninelevel.o \
	${CODDIR}netgain.o \
	${CODDIR}newdif.o \
	${CODDIR}newfil.o \
	${CODDIR}oldmodels.o \
	${CODDIR}p3fheader.o \
	${CODDIR}p4fheader.o \
	${CODDIR}p5fheader.o \
	${CODDIR}pheat.o \
	${CODDIR}phion.o \
	${CODDIR}phocrv.o \
	${CODDIR}photo1.o \
	${CODDIR}photo2.o \
	${CODDIR}photo3.o \
	${CODDIR}photo4.o \
	${CODDIR}photo5.o \
	${CODDIR}photsou.o \
	${CODDIR}popcha.o \
	${CODDIR}preion.o \
	${CODDIR}qua.o \
	${CODDIR}quart.o \
	${CODDIR}rankhug.o \
	${CODDIR}ratec.o \
	${CODDIR}recom.o \
	${CODDIR}reson.o \
	${CODDIR}reson2.o \
	${CODDIR}roots.o \
	${CODDIR}rk43sph.o \
	${CODDIR}s3fheader.o \
	${CODDIR}s2fheader.o \
	${CODDIR}sdifeq.o \
	${CODDIR}shock1.o \
	${CODDIR}shock2.o \
	${CODDIR}shock3.o \
	${CODDIR}shock4.o \
	${CODDIR}sinsla.o \
	${CODDIR}sixlevel.o \
	${CODDIR}slab.o \
	${CODDIR}spectrum.o \
	${CODDIR}spec2.o \
	${CODDIR}speclocal.o \
	${CODDIR}sphderivs.o \
	${CODDIR}sphere1.o \
	${CODDIR}spotap.o \
	${CODDIR}sumdata.o \
	${CODDIR}taudist.o \
	${CODDIR}teequi.o \
	${CODDIR}timion.o \
	${CODDIR}timtqui.o \
	${CODDIR}totphot.o \
	${CODDIR}transferline.o \
	${CODDIR}twophoton.o \
	${CODDIR}velshock.o \
	${CODDIR}vphotsou.o \
	${CODDIR}wbal.o \
	${CODDIR}wionpop.o \
	${CODDIR}wionabal.o \
	${CODDIR}wpsou.o \
	${CODDIR}wpsoufile.o \
	${CODDIR}wmodel.o \
	${CODDIR}wrsppop.o \
	${CODDIR}zer.o \
	${CODDIR}zerbuf.o \
	${CODDIR}zetaeff.o

#-------------------------------
#	ADDITIONAL FILES
#	${CODDIR}p6fheader.o \
#	${CODDIR}photo6.o \
#	${CODDIR}compph6.o \
#	${CODDIR}testroutines.o \
#
#-------------------------------
#-------  implicit rules  -----
#-------------------------------
#
.f.o: ${INCS}
	${FC} ${OPTS} -o $*.o $<

.c.o:
	${CC} -c $<
#
#-------------------------------
#------  explicit rules  ------
#-------------------------------
#
#
#-------------------------------
#---------  Targets  ----------
#-------------------------------
#
help:
	@echo MAPPINGS III v1.0.0m make options:
	@echo ' '
	@echo Type "'make help'    to see this menu"
	@echo Type "'make build'   to create executable from scratch and clean"
	@echo ' '
	@echo Type "'make prepare' to copy over '*.f' and '*.inc' files"
	@echo Type "'make compile' to create executable"
	@echo Type "'make backup'  to backup '*.f' and '*.inc' files"
	@echo Type "'make clean'   to backup '*.f' and '*.inc', remove '*.o'"
	@echo Type "'make listing' to create listing of code"
#
#-----------------------------------------------------------
#
prepare:
	@echo ' Copying over .f files'
	@cp ${BACKDIR}*.f ${CODDIR}
	@echo ' Copying over include files'
	@cp ${BACKDIR}*.inc ${INCDIR}
	@echo ' Cleaning up src area'
	@rm -f *~
	@echo ' Done.'
#
#-----------------------------------------------------------
#
clean:
	@echo ' Removing object files'
	@rm -f ${OBJ}
	@echo ' Backing up and cleaning up .f files'
	@mv $(OBJ:.o=.f) ${BACKDIR}
	@echo ' Backing up and cleaning up include files'
	@mv $(INCS) ${BACKDIR}
	@echo ' Cleaning up src area'
	@rm -f *~
	@echo ' Done.'
#
#-----------------------------------------------------------
#
backup:
	@echo ' Backing up .f files'
	@cp $(OBJ:.o=.f) ${BACKDIR}
	@echo ' Backing up include files'
	@cp $(INCS) ${BACKDIR}
	@echo ' Cleaning up src area'
	@rm -f *~
	@echo ' Done.'
#
#-----------------------------------------------------------
#
listing:
	@cat Makefile map.prefs ${INCS} $(OBJ:.o=.f) ${CODDIR}mapinitmac.f  ${CODDIR}abechamac.f > map3r.lis
#
#------------------------------------------------------------
#
compile:${EXEDIR}${OUTNAME}

${EXEDIR}${OUTNAME}: ${INCS} ${OBJ}
	${FC} ${LDR} -o ${EXEDIR}${OUTNAME} ${OBJ} ${LIB}

#
#------------------------------------------------------------
#
build:
	@echo ' Copying over .f files'
	@cp ${BACKDIR}*.f ${CODDIR}
	@echo ' Copying over include files'
	@cp ${BACKDIR}*.inc ${INCDIR}
	@echo ' Cleaning up src area'
	@echo ' Compiling...'
	@make compile
	@echo ' Removing object files'
	@rm -f ${OBJ}
	@echo ' Backing up and cleaning up .f files'
	@mv $(OBJ:.o=.f) ${BACKDIR}
	@echo ' Backing up and cleaning up include files'
	@mv $(INCS) ${BACKDIR}
	@echo ' Cleaning up src area'
	@rm -f *~
	@echo ' Done.'

