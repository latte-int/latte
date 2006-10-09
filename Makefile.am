ACLOCAL_AMFLAGS = -I m4

SUBDIRS = code EXAMPLES

docdir = $(pkgdatadir)/doc

dist_doc_DATA = manual_v1.1.pdf


cvstag:
	cd $(top_srcdir) && cvs tag $(CVSTAGFLAGS) version_`echo $(VERSION) | sed 's/[-.]/_/g'` 

