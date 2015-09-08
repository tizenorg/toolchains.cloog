cloog_srcdir = @CLOOG_SRCDIR@
cloog_builddir = @CLOOG_BUILDDIR@

if BUNDLED_CLOOG
    MAYBE_CLOOG = cloog-core
endif

SUBDIRS = . test
DIST_SUBDIRS = $(MAYBE_CLOOG) . test

ACLOCAL_AMFLAGS = -I m4

FORCE:
cloog-core/source/version.c: FORCE
	cd cloog-core; $(MAKE) $(AM_MAKEFLAGS) source/version.c
cloog-core/version.h: FORCE
	cd cloog-core; $(MAKE) $(AM_MAKEFLAGS) version.h
cloog-core/include/cloog/version.h: FORCE
	cd cloog-core; $(MAKE) $(AM_MAKEFLAGS) include/cloog/version.h

bin_PROGRAMS = cloog
lib_LTLIBRARIES = libcloog-ppl.la

if NEED_GET_MEMORY_FUNCTIONS
GET_MEMORY_FUNCTIONS=$(cloog_srcdir)/source/mp_get_memory_functions.c
endif

SOURCES_CORE = \
	$(GET_MEMORY_FUNCTIONS) \
	$(cloog_srcdir)/source/block.c \
	$(cloog_srcdir)/source/clast.c \
	$(cloog_srcdir)/source/matrix.c \
	$(cloog_srcdir)/source/state.c \
	$(cloog_srcdir)/source/input.c \
	$(cloog_srcdir)/source/int.c \
	$(cloog_srcdir)/source/loop.c \
	$(cloog_srcdir)/source/names.c \
	$(cloog_srcdir)/source/options.c \
	$(cloog_srcdir)/source/pprint.c \
	$(cloog_srcdir)/source/program.c \
	$(cloog_srcdir)/source/statement.c \
	$(cloog_srcdir)/source/stride.c \
	$(cloog_srcdir)/source/union_domain.c \
	$(cloog_srcdir)/source/matrix/constraintset.c \
	$(cloog_builddir)/source/version.c \
	$(cloog_builddir)/version.h

DEFAULT_INCLUDES = -I.
INCLUDES = -I$(top_builddir)/include -I$(top_srcdir)/include
AM_CFLAGS = $(CFLAGS_WARN)
libcloog_ppl_la_CPPFLAGS = @CLOOG_CPPFLAGS@ @PPL_CPPFLAGS@
libcloog_ppl_la_LDFLAGS = -version-info @versioninfo@ \
	-rpath $(libdir) @PPL_LDFLAGS@
libcloog_ppl_la_LIBADD = @PPL_LIBS@
nodist_libcloog_ppl_la_SOURCES = $(SOURCES_CORE)
libcloog_ppl_la_SOURCES = \
	source/domain.c \
	source/backend.c
LDADD = libcloog-ppl.la
cloog_CPPFLAGS = @CLOOG_CPPFLAGS@
cloog_DEPENDENCIES = libcloog-ppl.la
nodist_cloog_SOURCES = $(cloog_srcdir)/source/cloog.c

cloogincludedir = $(includedir)/cloog
nodist_clooginclude_HEADERS = \
	$(cloog_srcdir)/include/cloog/block.h \
	$(cloog_srcdir)/include/cloog/clast.h \
	$(cloog_srcdir)/include/cloog/cloog.h \
	$(cloog_srcdir)/include/cloog/input.h \
	$(cloog_srcdir)/include/cloog/int.h \
	$(cloog_srcdir)/include/cloog/matrix.h \
	$(cloog_srcdir)/include/cloog/state.h \
	$(cloog_srcdir)/include/cloog/domain.h \
	$(cloog_srcdir)/include/cloog/loop.h \
	$(cloog_srcdir)/include/cloog/constraints.h \
	$(cloog_srcdir)/include/cloog/names.h \
	$(cloog_srcdir)/include/cloog/options.h \
	$(cloog_srcdir)/include/cloog/pprint.h \
	$(cloog_srcdir)/include/cloog/program.h \
	$(cloog_srcdir)/include/cloog/statement.h \
	$(cloog_srcdir)/include/cloog/stride.h \
	$(cloog_srcdir)/include/cloog/union_domain.h \
	$(cloog_builddir)/include/cloog/version.h

cloogmatrixincludedir = $(cloogincludedir)/matrix
nodist_cloogmatrixinclude_HEADERS = \
	$(cloog_srcdir)/include/cloog/matrix/constraintset.h

cloogpplincludedir = $(cloogincludedir)/ppl
cloogpplinclude_HEADERS = \
	include/cloog/ppl/backend.h \
	include/cloog/ppl/cloog.h \
	include/cloog/ppl/domain.h

check:
	$(MAKE) test -C test

install-data-local:
	@test -z "$(pkgconfig_libdir)" || $(mkdir_p) "$(DESTDIR)$(pkgconfig_libdir)"
	$(INSTALL_DATA) $(pkgconfig_libfile) "$(DESTDIR)$(pkgconfig_libdir)/$(pkgconfig_libfile)"

uninstall-local:
	rm -f "$(DESTDIR)$(pkgconfig_libdir)/$(pkgconfig_libfile)"
