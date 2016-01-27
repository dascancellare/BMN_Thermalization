ACLOCAL_AMFLAGS=-I m4

bin_PROGRAMS=$(top_builddir)/bin/perturb $(top_builddir)/bin/jfix_finder
__top_builddir__bin_perturb_SOURCES= \
        src/conf.cpp \
        src/gauge_fix.cpp \
        src/matr.cpp \
        src/observables.cpp \
        src/random.cpp \
        src/theory.cpp \
        src/tools.cpp \
        src/perturb.cpp \
        src/update.cpp

__top_builddir__bin_jfix_finder_SOURCES= \
        src/jfix_finder.cpp

CXXFLAGS+=-I/usr/include/root -pthread -m64 -I/usr/include/root

LDADD=-lGui
LDADD+=-lCore
LDADD+=-lCint
LDADD+=-lRIO
LDADD+=-lNet
LDADD+=-lHist
LDADD+=-lGraf
LDADD+=-lGraf3d
LDADD+=-lGpad
LDADD+=-lTree
LDADD+=-lRint
LDADD+=-lPostscript
LDADD+=-lMatrix
LDADD+=-lPhysics
LDADD+=-lMathCore
LDADD+=-lThread
LDADD+=-lMinuit
LDFLAGS=-rdynamic -pthread
