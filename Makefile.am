LIBDIR=@LIBDIR@

ACLOCAL_AMFLAGS = -I m4

AUTOMAKE_OPTIONS = foreign
SUBDIRS = src examples doc

# DOXYGEN SUPPORT
include $(top_srcdir)/aminclude.am

# ensure the distribution of the doxygen configuration file
EXTRA_DIST = doxygen.cfg