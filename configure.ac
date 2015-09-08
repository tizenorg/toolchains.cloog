AC_PREREQ(2.13)
AC_INIT(source/backend.c)
AC_CONFIG_AUX_DIR(autoconf)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE(cloog-parma, 0.16.1)
AC_SUBST(versioninfo)
versioninfo=1:1:0

AC_PROG_CC
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_CHECK_PROG(CD, cd)
AC_PROG_INSTALL
AC_PROG_LIBTOOL

AX_CC_MAXOPT
AC_SUBST(CFLAGS_WARN)
AX_CFLAGS_WARN_ALL(CFLAGS_WARN)

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

dnl Checks for header files.
AC_HEADER_STDC

dnl /**************************************************************************
dnl  *                            Where is GMP?                               *
dnl  **************************************************************************/

AX_SUBMODULE(gmp,system,system)

case "$with_gmp" in
system)
    if test "x$with_gmp_prefix" != "x"; then
	CPPFLAGS="-I$with_gmp_prefix/include $CPPFLAGS"
    fi
    
    if test "$with_gmp_exec_prefix" != "yes" ; then
	LDFLAGS="-L$with_gmp_exec_prefix/lib $LDFLAGS"
    fi
     
    AC_CHECK_HEADER(gmp.h,
                    [],
                    [AC_MSG_ERROR(Can't find gmp headers.)])
    AC_CHECK_LIB(gmp,
                 __gmpz_init,
                 [LIBS="$LIBS -lgmp"],
                 [AC_MSG_ERROR(Can't find gmp library.)])
    AC_CHECK_DECLS(mp_get_memory_functions,[],[
	    need_get_memory_functions=true
    ],[#include <gmp.h>])
    ;;
esac
AM_CONDITIONAL(NEED_GET_MEMORY_FUNCTIONS,
		test x$need_get_memory_functions = xtrue)

dnl /**************************************************************************
dnl  *                          Where is PPL?                                 *
dnl  **************************************************************************/

AX_SUBMODULE(ppl,system,system)

AC_SUBST(PPL_CPPFLAGS)
AC_SUBST(PPL_LDFLAGS)
AC_SUBST(PPL_LIBS)
case "$with_ppl" in
system)
	if test "x$with_ppl_prefix" != "x"; then
		PPL_CPPFLAGS="-I$with_ppl_prefix/include"
	fi
	if test "x$with_ppl_exec_prefix" != "x"; then
		PPL_LDFLAGS="-L$with_ppl_exec_prefix/lib"
	fi
	PPL_LIBS="-lppl_c -lppl -lgmpxx"
	SAVE_CPPFLAGS="$CPPFLAGS"
	CPPFLAGS="$PPL_CPPFLAGS $CPPFLAGS"
	AC_CHECK_HEADER(ppl_c.h,
			[],
			[AC_MSG_ERROR(Can't find PPL headers.)])
	CPPFLAGS="$SAVE_CPPFLAGS"
esac

dnl /**************************************************************************
dnl  *                          Where is CLooG?                               *
dnl  **************************************************************************/

AX_SUBMODULE(cloog,build|bundled,bundled)

AC_SUBST(CLOOG_SRCDIR)
AC_SUBST(CLOOG_BUILDDIR)
AC_SUBST(CLOOG_CPPFLAGS)
case "$with_cloog" in
bundled)
	CLOOG_SRCDIR="\$(top_srcdir)/cloog-core"
	CLOOG_BUILDDIR="\$(top_builddir)/cloog-core"

	;;
build)
	case "$cloog_srcdir" in
	/*)
		CLOOG_SRCDIR=$cloog_srcdir
		;;
	*)
		CLOOG_SRCDIR="\$(top_srcdir)/$cloog_srcdir"
	esac
	case "$with_cloog_builddir" in
	/*)
		CLOOG_BUILDDIR=$with_cloog_builddir
		;;
	*)
		CLOOG_BUILDDIR="\$(top_builddir)/$with_cloog_builddir"
	esac
	;;
esac
CLOOG_CPPFLAGS="-I$CLOOG_SRCDIR/include -I$CLOOG_BUILDDIR -I$CLOOG_BUILDDIR/include"
AM_CONDITIONAL(BUNDLED_CLOOG, test $with_cloog = bundled)


AC_DEFINE([CLOOG_INT_GMP], 1, [Use arbitrary precision integers])
BITS="gmp"


dnl /**************************************************************************
dnl  *                            Substitutions                               *
dnl  **************************************************************************/


dnl Substitutions to do in Makefile.in.
AC_SUBST(CC)
AC_SUBST(LN_S)

AC_SUBST(prefix)
AC_SUBST(exec_prefix)
AC_SUBST(INSTALL)

AC_SUBST(BITS)

PACKAGE_NAME="cloog-ppl"
PACKAGE_CFLAGS="-DCLOOG_INT_GMP=1"
AX_CREATE_PKGCONFIG_INFO

AC_CONFIG_FILES(Makefile test/Makefile)
if test $with_cloog = bundled; then
	AC_CONFIG_SUBDIRS(cloog-core)
fi
AC_CONFIG_COMMANDS_POST([
	old_args=""
	for arg in $ac_configure_args; do
		case arg in
		--with-ppl*)
			;;
		*)
			old_args="$old_args $arg"
			;;
		esac
	done
	ac_configure_args="$old_args --without-isl $cloog_ppl_configure_args"
])
AC_OUTPUT


echo "             /*-----------------------------------------------*"
echo "              *           CLooG configuration is OK           *"
echo "              *-----------------------------------------------*/"
echo "It appears that your system is OK to start CLooG compilation. You need"
echo "now to type \"make\". After compilation, you should check CLooG by typing"
echo "\"make check\". If no problem occur, you can type \"make uninstall\" if"
echo "you are upgrading an old version. Lastly type \"make install\" to install"
echo "CLooG on your system (log as root if necessary)."
