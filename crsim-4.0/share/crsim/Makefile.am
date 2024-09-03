# ---------------------------------------------------------------------------
# COPYRIGHT    : Mariko Oue and Aleksandra Tatarevic
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
#
nobase_dist_pkgdata_DATA = \
                aux/CHANGE_LOG \
                aux/LLUTS3/
                aux/P3/INFO \
                aux/P3/*dat\
                aux/SBM_20/INFO\
                aux/SBM_20/bulk*\
                aux/SBM_20/*asc\
                test/inp/wrf_MP10/PARAMETERS \
                test/inp/wrf_MP50/PARAMETERS \
                test/inp/rams_MP40/PARAMETERS \
                test/out_ref/


testdir=$(pkgdatadir)/test/
logs_refdir=$(pkgdatadir)/test/logs_ref
logsdir=$(pkgdatadir)/test/logs
outdir=$(pkgdatadir)/test/out

dist_logs_ref_DATA = test/logs_ref/*.log
dist_logs_DATA=test/logs/*ME
dist_out_DATA=test/out/*ME

