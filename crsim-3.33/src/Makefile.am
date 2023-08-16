# ---------------------------------------------------------------------------
# COPYRIGHT    : Aleksandra Tatarevic
# ----------------------------------------------------------------------------
# PROJECT      : CRSIM
# FILE         : Makefile.am
# VERSION      : $Revision: $
# DATE         : $Date: $
# TYPE         : Automake file
# IDENTIFIER   : $Id: $
# ----------------------------------------------------------------------------
#
# This file contains the basic configuration for the generation of a Makefile
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
AM_FCFLAGS = $(EXTRA_FCFLAGS) $(NETCDF_FORTRAN_FCFLAGS)
AM_CPPFLAGS = -DNOMAGICK -DNOHDF -DUNIX $(NETCDF_CPPFLAGS)
AM_CFLAGS = $(EXTRA_CFLAGS)
AM_LDFLAGS =  $(NETCDF_FORTRAN_LDFLAGS) $(NETCDF_LDFLAGS)
#bin_SCRIPTS = $(top_builddir)/bin/run_crsim.sh
bin_PROGRAMS = $(top_builddir)/bin/crsim

__top_builddir__bin_crsim_SOURCES = phys_param_mod.f90 crsim_luts_mod.f90 crsim_mod.f90 postprocess_mod.f90 wrf_var_mod.f90 wrf_rvar_mod.f90 crsim.f90 crsim_subrs.f90 crsim_subrs_postprocess.f90 ReadConfParameters.f90 ReadInpRAMSFile.f90 ReadInpWRFFile.f90 WriteOutNetcdf.f90

CLEANFILES = *.mod *.o 
