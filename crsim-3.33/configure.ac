# ---------------------------------------------------------------------------
# COPYRIGHT    : Aleksandra Tatarevic
# ----------------------------------------------------------------------------
# PROJECT      : CRSIM
# FILE         : configure.ac
# VERSION      : $Revision: $
# DATE         : $Date: $
# TYPE         : Automake file
# IDENTIFIER   : $Id: $
# ----------------------------------------------------------------------------
#
# This file contains the basic configuration for the generation of a configure
# script
# ----------------------------------------------------------------------------
#
# LICENSE
# ----------------------------------------------------------------------------
#
# You may obtain a copy of the License distributed with this s/w package, 
# software distributed under the License is distributed on an “AS IS” BASIS, 
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See
# the License for the specific language governing permissions and limitations 
# under the License.
# ----------------------------------------------------------------------------
#
#
# Running autoconf on this file will trigger a warning if
# autoconf is not at least the specified version.
AC_PREREQ([2.63])

# Initialize with name, version, and support email address.
AC_INIT([crsim],[3.32], [aleksandra.tatarevic@mcgill.ca])

# Check for the existance of a file in the source directory.
AC_CONFIG_SRCDIR([src/crsim.f90])
AC_CONFIG_HEADERS([config.h])

# Folder where some additional scripts are stored
AC_PREFIX_DEFAULT([`pwd`])
AC_CONFIG_AUX_DIR([build-aux])

# This call is required by automake.
AM_INIT_AUTOMAKE([foreign subdir-objects tar-ustar])

# For using custom m4 macros.
AC_CONFIG_MACRO_DIR([m4])

# Disable default CFLAGS='-g -O2' and FCFLAGS='-g -O2' set by AC_PROG_CC and AC_PROG_FC
: ${CFLAGS=""}
: ${FCFLAGS=""}

# Checks for programs.
AC_PROG_CC
AC_PROG_FC

# sets the variable OPENMP_FCFLAGS  to the fortran  compiler flags needed for supporting OpenMP
AC_OPENMP

# Checks for libraries.
AC_CHECK_LIB([m], [pow])
AC_CHECK_LIB([z], [uncompress])
# Checks for header files.
AC_CHECK_HEADERS([fcntl.h float.h inttypes.h limits.h malloc.h stddef.h stdlib.h string.h sys/statvfs.h sys/time.h unistd.h])
AC_CHECK_HEADERS([zlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_MODE_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([fdatasync gethostname memset mkdir pow putenv sqrt strchr strdup strerror strrchr strspn strstr strtol strtoul uname])

# debugging support
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [enable debugging support @<:@default = no@:>@])],
              [],
              [enable_debug=no])


# netCDF support
AX_CHECK_NETCDF([have_netcdf=yes], [have_netcdf=no])
AX_CHECK_NETCDF_FORTRAN([have_netcdf=yes], [have_netcdf=no])

if test "x$have_netcdf" = xno; then
    AC_MSG_ERROR([
  -----------------------------------------------------------------------------
   The netCDF C and Fortran libraries, header files, and modules are required
   to build $PACKAGE_NAME. They seem to be missing or incomplete on your system.
   Stopping the configuration now. Check 'config.log' for more information.
  -----------------------------------------------------------------------------])
fi

#AC_CHECK_LIB([gomp], [omp_get_wtime], [], [
#	echo "OpenMP library is mandatory for $PACKAGE_NAME. Stopping the configuration now."
#	exit -1
#])

# NCCMP tool
AX_CHECK_NCCMP([have_nccmp=yes], [have_nccmp=no])
if test "x$have_nccmp" = xno; then
    AC_MSG_NOTICE([
  -----------------------------------------------------------------------------
   The NCCMP is required to test $PACKAGE_NAME. It seems to be
   missing or incomplete on your system. Please be aware that one or more 
   tests will fail because the comparison between netcdf files will not be possible.
  -----------------------------------------------------------------------------])
fi
AC_SUBST([NCCMP])

# disabling  mpi
AC_ARG_ENABLE([openmp],
              [AS_HELP_STRING([--disable-openmp],
                              [disable openmp support @<:@default = no@:>@])],
              [disable_openmp=yes],
              [disable_openmp=no])

if test "x$disable_openmp" = xyes; then
    AC_MSG_NOTICE([
  -----------------------------------------------------------------------------
   The OpenMP library is mandatory for $PACKAGE_NAME. However, the option 
   '--disable-openmp' has been specified when invoking 'configure'. Note: 
   this wil reject OpenMP support.
  -----------------------------------------------------------------------------])
fi


# preload lut option
AC_ARG_WITH([preloadlut],
              [AS_HELP_STRING([--without-preloadlut],
                              [disable support for  preloading LUTs @<:@default = no@:>@])],
              [],
              [with_preloadlut=yes])




if test "x$with_preloadlut" != xyes; then
    AC_MSG_NOTICE([
  -----------------------------------------------------------------------------
   The option '--without-preloadlut' is specified when invoking 'configure'.
   Note: This option turns off preloading of all scattering lookup tables
   (LUTs) before starting running the main code.  
  -----------------------------------------------------------------------------])
if test "x$disable_openmp" != xyes; then

    AC_MSG_ERROR([
  -----------------------------------------------------------------------------
    You have specified the optiopn '--without-preloadlut'.  You have also to
    reject the OpenMP support by specifying the otion '--disable-openmp' when 
    invoking 'configure'.
  -----------------------------------------------------------------------------])
fi


fi


# tuning Fortran compiler flags (FC is set by AC_PROG_FC)
if test "x$FC" = xgfortran; then
#-std=f2003
     EXTRA_FCFLAGS=$OPENMP_CFLAGS" -O3 -cpp -Wunused -Wuninitialized -Waliasing "
     EXTRA_FCFLAGS1="-D__PRELOAD_LUT__"
#     EXTRA_FCFLAGS=$OPENMP_CFLAGS" -O3 -cpp -Wunused -Wuninitialized -Waliasing -Wconversion -Wsurprising -Wintrinsics-std -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,invalid"


#    EXTRA_FCFLAGS="-Wall -Wextra -cpp -ffree-line-length-none -fall-intrinsics -fno-range-check -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant -Wunused -Wno-maybe-uninitialized -Wuninitialized -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,invalid -O2"

       if test "x$with_preloadlut" != xyes; then
          EXTRA_FCFLAGS1=""
       fi

     EXTRA_FCFLAGS="$EXTRA_FCFLAGS $EXTRA_FCFLAGS1"
fi
AC_SUBST([EXTRA_FCFLAGS])

AC_CONFIG_FILES([Makefile
                 etc/Makefile
                 doc/Makefile
                 share/crsim/Makefile
                 share/crsim/test/Makefile
                 src/Makefile])

AC_OUTPUT

AC_MSG_RESULT([
--------------------------------------------------------------------------------
    Configure Information:

    $PACKAGE_NAME version $PACKAGE_VERSION

    Prefix (installation path)........................: $prefix

    C Compiler........................................: $CC
       C preprocessor flags (CPP).....................: $CPP
       C compiler flags (CFLAGS)......................: $CFLAGS
       Extra C flags..................................: $EXTRA_CFLAGS

    C/C++ preprocessor flags (CPPFLAGS)...............: $CPPFLAGS $NETCDF_CPPFLAGS

    Fortran Compiler..................................: $FC 
       Fortran compiler flags (FCFLAGS)...............: $FCFLAGS
       Extra Fortran flags............................: $EXTRA_FCFLAGS

    Fortran preprocessor flags........................: $NETCDF_FORTRAN_FCFLAGS

    Linker............................................: $LD
       Linker flags (LDFLAGS).........................: $LDFLAGS $NETCDF_LDFLAGS 
       Linker flags (fortran only)....................: $NETCDF_FORTRAN_LDFLAGS
       Linker Libraries (LIBS)........................: $LIBS

    preloading of luts................................: $with_preloadlut
    disable openmp....................................: $disable_openmp
    Fortan flags needed for supporting OpenMP.........: $OPENMP_CFLAGS
    NCCMP.............................................: $NCCMP
--------------------------------------------------------------------------------
])
