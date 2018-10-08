# Makefile
#
# This file is part of time_clean, a CASA task for cleaning images in which the source varies significantly in flux over the time of the observation.
#
# See the file COPYRIGHT for authorship and intellectual property rights.

# Paths:
srcdir		= ${CURDIR}/src
datadir		= ${CURDIR}/data

ifneq (,$(wildcard ${CURDIR}/lib/.))
    LIBS += -L${CURDIR}/lib
endif
ifneq (,$(wildcard ${HOME}/lib/.))
    LIBS += -L${HOME}/lib
endif
ifneq (,$(wildcard /opt/local/lib/.))
    LIBS += -L/opt/local/lib
endif
ifneq (,$(wildcard /sw/lib/.))
    LIBS += -L/sw/lib
endif
ifneq (,$(wildcard /usr/local/lib/.))
    LIBS += -L/usr/local/lib
endif

#CCFLAGS += -I/usr/include/cfitsio

CPPFLAGS = -I${CURDIR}/include \
	   -I${CURDIR}/src \
	   -I${HOME}/include \
	   -I/opt/local/include \
	   -I/sw//include

CPPFLAGS += -I/usr/include/cfitsio

# Names of source files included:
include Makefile.srcs

##
## Do not change anything below unless you know what you are doing! 
##

TARGET  = sw_clean
CC	= gcc

CCFLAGS += -O3 -falign-loops=16 -fno-strict-aliasing
LDFLAGS += -lgsl -lgslcblas -lcfitsio -lfftw3 -lm 

# Comment out or remove the following 2 lines when debugging is no longer required:
CCFLAGS += -DTEST
CC += -g -Wall

OBJS = $(SRCS:.c=.o)

.PHONY: all clean objclean distclean

all:: ${TARGET} 

# Implicit rules:
%.o : %.c
	${CC} ${CCFLAGS} ${CPPFLAGS} -o $@ -c $<

${OBJS} : ${INCS}

${TARGET}: ${OBJS}
	${CC} -o $@ $^ ${LIBS} ${LDFLAGS}

objclean::
	rm -f ${srcdir}/*.o ${srcdir}/*/*.o

clean:: objclean
	rm -f ${TARGET}
	rm -f *~ ${srcdir}/*~ ${srcdir}/*/*~

distclean:: clean
	rm ${datadir}/test_CC.fits ${datadir}/test_resids.fits

