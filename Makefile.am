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
	$(top_builddir)/bin/random_walk \
	$(top_builddir)/bin/relax

__top_builddir__bin_perturb_SOURCES= \
	$(common_sources) \
	src/perturb.cpp

__top_builddir__bin_pert_edge_SOURCES= \
	$(common_sources) \
	src/pert_edge.cpp

__top_builddir__bin_random_walk_SOURCES= \
	$(common_sources) \
	src/random_walk.cpp

__top_builddir__bin_relax_SOURCES= \
	$(common_sources) \
	src/relax.cpp
