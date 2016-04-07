ACLOCAL_AMFLAGS=-I m4

common_sources= \
        src/conf.cpp \
        src/gauge_fix.cpp \
        src/matr.cpp \
        src/observables.cpp \
        src/random.cpp \
        src/theory.cpp \
        src/tools.cpp  \
        src/update.cpp

bin_PROGRAMS= \
	$(top_builddir)/bin/perturb \
	$(top_builddir)/bin/pert_edge \
	$(top_builddir)/bin/jfix_finder \
	$(top_builddir)/bin/relax

__top_builddir__bin_perturb_SOURCES= \
	$(common_sources) \
	src/perturb.cpp

__top_builddir__bin_pert_edge_SOURCES= \
	$(common_sources) \
	src/pert_edge.cpp

__top_builddir__bin_relax_SOURCES= \
	$(common_sources) \
	src/relax.cpp

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
