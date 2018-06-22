# Include commands for every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables: 
noinst_PROGRAMS= main

# Local sources
main_SOURCES = main.cc

# Requried libraries
main_LDADD = -L@libdir@ -lbilly -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Include path for library headers
AM_CPPFLAGS += -I@includedir@ -I$(abs_top_builddir)/user_drivers/billy/src
