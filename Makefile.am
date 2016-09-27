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
	$(top_builddir)/bin/pert_edge

__top_builddir__bin_pert_edge_SOURCES= \
	$(common_sources) \
        src/pert_edge.cpp
