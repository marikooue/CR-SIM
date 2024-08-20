# Add --with-nccmp=... option to configure

# AX_CHECK_NCCMP([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ---------------------------------------------------------
#
# Upon success, sets the executable path for nccmp
AC_DEFUN([AX_CHECK_NCCMP], [
AC_ARG_WITH([nccmp],
           [AS_HELP_STRING([--with-nccmp=DIR],
                           [use nccmp executable from directory DIR])])
AS_IF([test "x$with_nccmp" != x],
      [AC_CHECK_PROG([ax_have_nccmp], [nccmp], [yes], [no], [$with_nccmp])
       AS_IF([test "x$ax_have_nccmp" == xyes],
              [NCCMP="$with_nccmp/nccmp" ax_have_nccmp=yes],
              [ax_have_nccmp=no])],
      [AC_CHECK_PROG([ax_have_nccmp], [nccmp], [yes], [no])
       AS_IF([test "x$ax_have_nccmp" == xyes],
              [NCCMP=nccmp ax_have_nccmp=yes],
              [ax_have_nccmp=no])] 
      )
AS_IF([test "x$ax_have_nccmp" = xyes],
     [AC_SUBST([NCCMP])
      $1],
     [$2])
])
dnl vim:set softtabstop=4 shiftwidth=4 expandtab:

