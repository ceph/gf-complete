# gf-complete autoconf template

# FIXME - add project url as the last argument
AC_INIT(gf-complete, 1.0)

AC_PREREQ([2.61])

AM_INIT_AUTOMAKE(no-dependencies)
LT_INIT # libtool

AC_CONFIG_HEADER(include/config.h)

dnl Needed when reconfiguring with 'autoreconf -i -s'
AC_CONFIG_MACRO_DIR([m4])

dnl Compiling with per-target flags requires AM_PROG_CC_C_O.
AC_PROG_CC

AX_EXT()

AC_CONFIG_FILES([Makefile src/Makefile tools/Makefile test/Makefile examples/Makefile])
AC_OUTPUT
