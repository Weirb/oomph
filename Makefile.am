SUBDIRS = \
SourceFiles

#---------------------------------------------------------------------------

# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS = main

# Name of executable
check_PROGRAMS = main

#---------------------------------------------------------------------------

# Local sources that Joe's Poisson code depends on:
main_SOURCES = main.cc

# Required libraries
main_LDADD = -L@libdir@ -lpml_helmholtz -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
AM_CPPFLAGS += -I@includedir@  
