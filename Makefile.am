ACLOCAL_AMFLAGS=-I m4

bin_PROGRAMS=$(top_builddir)/bin/main $(top_builddir)/bin/jfix_finder
__top_builddir__bin_main_SOURCES= \
        src/conf.cpp \
        src/main.cpp \
        src/matr.cpp \
        src/observables.cpp \
        src/random.cpp \
        src/theory.cpp \
        src/tools.cpp \
        src/update.cpp

__top_builddir__bin_jfix_finder_SOURCES= \
        src/jfix_finder.cpp

__top_builddir__bin_jfix_finder_CXXFLAGS= \
	-Wall -O3 -fPIC -I/usr/include/root -pthread -m64 -I/usr/include/root
__top_builddir__bin_jfix_finder_LDADD=-lGui
__top_builddir__bin_jfix_finder_LDADD+=-lCore
__top_builddir__bin_jfix_finder_LDADD+=-lCint
__top_builddir__bin_jfix_finder_LDADD+=-lRIO
__top_builddir__bin_jfix_finder_LDADD+=-lNet
__top_builddir__bin_jfix_finder_LDADD+=-lHist
__top_builddir__bin_jfix_finder_LDADD+=-lGraf
__top_builddir__bin_jfix_finder_LDADD+=-lGraf3d
__top_builddir__bin_jfix_finder_LDADD+=-lGpad
__top_builddir__bin_jfix_finder_LDADD+=-lTree
__top_builddir__bin_jfix_finder_LDADD+=-lRint
__top_builddir__bin_jfix_finder_LDADD+=-lPostscript
__top_builddir__bin_jfix_finder_LDADD+=-lMatrix
__top_builddir__bin_jfix_finder_LDADD+=-lPhysics
__top_builddir__bin_jfix_finder_LDADD+=-lMathCore
__top_builddir__bin_jfix_finder_LDADD+=-lThread
__top_builddir__bin_jfix_finder_LDADD+=-lMinuit
__top_builddir__bin_jfix_finder_LDFLAGS=-rdynamic -pthread
