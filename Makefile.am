ACLOCAL_AMFLAGS = -I m4 -I m4/gnulib

SUBDIRS = code EXAMPLES doc

docdir = $(pkgdatadir)/doc

dist_doc_DATA = NEWS README


cvstag:
	cd $(top_srcdir) && cvs tag $(CVSTAGFLAGS) version_`echo $(VERSION) | sed 's/[-.]/_/g'` 

update-changelog: 
	(cd $(srcdir) && svn2cl --output=ChangeLog.svn2cl --authors=.svn2cl-authors -r HEAD:{2012-01-01} https://polyhedron.math.ucdavis.edu/usvn/svn/latte/trunk)

#svn2cl --output=ChangeLog.version_1_4_branch.svn2cl --include-rev --stop-on-copy --authors=.svn2cl-authors https://polyhedron.math.ucdavis.edu/usvn/svn/latte/branches/version_1_4_branch
