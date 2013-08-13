ACLOCAL_AMFLAGS = -I m4 -I m4/gnulib

SUBDIRS = code EXAMPLES doc

docdir = $(pkgdatadir)/doc

dist_doc_DATA = NEWS README


cvstag:
	cd $(top_srcdir) && cvs tag $(CVSTAGFLAGS) version_`echo $(VERSION) | sed 's/[-.]/_/g'` 

