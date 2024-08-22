# ----------------------------------------------------------------------------

# AX_CHECK_NETCDF([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#
# Upon success, sets the variables NETCDF_LDFLAGS and NETCDF_CPPFLAGS.

dnl defines a custom macro
AC_DEFUN([AX_CHECK_NETCDF], [

      dnl provides a framework to handle the --with-{arg} values passed to configure on the command line      
      AC_ARG_WITH([netcdf],
            [AS_HELP_STRING([--with-netcdf=DIR], [use netCDF library from directory DIR])],
            netcdf_prefix="$with_netcdf"
            []
            )
      
      AS_IF([test x$netcdf_prefix != x],
            [AS_IF([test -d "$netcdf_prefix/lib64"],
                  [NETCDF_LDFLAGS="-L$netcdf_prefix/lib64 -Wl,-rpath,$netcdf_prefix/lib64"
                  NETCDF_CPPFLAGS="-I$netcdf_prefix/include"],
                  [test -d "$netcdf_prefix/lib"],
                  [NETCDF_LDFLAGS="-L$netcdf_prefix/lib -Wl,-rpath,$netcdf_prefix/lib"
                  NETCDF_CPPFLAGS="-I$netcdf_prefix/include"],
                  [AC_MSG_ERROR([
  -----------------------------------------------------------------------------
     --with-netcdf=$netcdf_prefix is not a valid directory
  -----------------------------------------------------------------------------])])],
      [AC_MSG_ERROR([
  -----------------------------------------------------------------------------
   Missing option `--with-netcdf=DIR`. Please pass the option to specify the
   location of the netCDF library
  -----------------------------------------------------------------------------])]
           )
      
      ax_saved_LDFLAGS=$LDFLAGS
      ax_saved_CPPFLAGS=$CPPFLAGS
      LDFLAGS="$NETCDF_LDFLAGS $LDFLAGS"
      CPPFLAGS="$NETCDF_CPPFLAGS $CPPFLAGS"
      ax_have_netcdf=yes
      dnl checks for netCDF library
      AC_SEARCH_LIBS([nc_open], [netcdf], [], [ax_have_netcdf=no])
      dnl checks for header files
      AC_CHECK_HEADERS([netcdf.h], [], [ax_have_netcdf=no])
      LDFLAGS=$ax_saved_LDFLAGS
      CPPFLAGS=$ax_saved_CPPFLAGS

      AS_IF([test "x$ax_have_netcdf" = xyes],
            dnl outputing netCDF library and header flags for C
            [AC_SUBST([NETCDF_LDFLAGS])
            AC_SUBST([NETCDF_CPPFLAGS])
            $1],
            [$2])
      ]
)

# AX_CHECK_NETCDF_FORTRAN([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#
# Upon success, sets the variables NETCDF_FORTRAN_LDFLAGS and
# NETCDF_FORTRAN_FCFLAGS.
AC_DEFUN([AX_CHECK_NETCDF_FORTRAN], [

      AC_ARG_WITH([netcdf-fortran],
            [AS_HELP_STRING([--with-netcdf-fortran=DIR], 
                            [use netCDF fortran library from directory DIR])],
            netcdf_fortran_prefix="$with_netcdf_fortran"
            [AS_IF([test "x$with_netcdf" != x], [with_netcdf_fortran=$with_netcdf])]
            )

      AS_IF([test x$netcdf_fortran_prefix != x],
            [AS_IF([test -d "$netcdf_fortran_prefix/lib64"],
                  [NETCDF_FORTRAN_LDFLAGS="-L$netcdf_fortran_prefix/lib64 -Wl,-rpath,$netcdf_fortran_prefix/lib64"
                  NETCDF_FORTRAN_FCFLAGS="-I$netcdf_fortran_prefix/include"],
                  [test -d "$netcdf_fortran_prefix/lib"],
                  [NETCDF_FORTRAN_LDFLAGS="-L$netcdf_fortran_prefix/lib -Wl,-rpath,$netcdf_fortran_prefix/lib"
                  NETCDF_FORTRAN_FCFLAGS="-I$netcdf_fortran_prefix/include"],
                  [AC_MSG_ERROR([
  -----------------------------------------------------------------------------
     --with-netcdf-fortran=$netcdf_fortran_prefix is not a valid directory
  -----------------------------------------------------------------------------])])],
      [AC_MSG_ERROR([
  -----------------------------------------------------------------------------
   Missing option `--with-netcdf-fortran=DIR`. Please pass the option to specify
   the location of the netCDF library
  -----------------------------------------------------------------------------])]
           )

      ax_saved_LDFLAGS=$LDFLAGS
      ax_saved_FCFLAGS=$FCFLAGS
      LDFLAGS="$NETCDF_FORTRAN_LDFLAGS $NETCDF_LDFLAGS $LDFLAGS"
      FCFLAGS="$NETCDF_FORTRAN_FCFLAGS $FCFLAGS"
      ax_have_netcdf_fortran=yes
      dnl defines language for checks
      AC_LANG_PUSH([Fortran])
      dnl checks for netCDF fortran library
      AC_SEARCH_LIBS([nf_open], [netcdff], [], [ax_have_netcdf_fortran=no])
      dnl checks for netCDF fortran modules
      AC_MSG_CHECKING([for Fortran module netcdf])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [      
            use netcdf
            character(len=*), parameter :: path = "dummy"
            integer :: mode, status, ncid
            status = nf90_open(path,mode,ncid)            
            ])],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no)
            ax_have_netcdf_fortran=no])
      AC_LANG_POP([Fortran])
      LDFLAGS=$ax_saved_LDFLAGS
      FCFLAGS=$ax_saved_FCFLAGS

      AS_IF([test "x$ax_have_netcdf_fortran" = xyes],
            dnl outputing netCDF library and module flags for fortran
            [AC_SUBST([NETCDF_FORTRAN_LDFLAGS])
            AC_SUBST([NETCDF_FORTRAN_FCFLAGS])
            $1],
            [$2])
      ]
)
