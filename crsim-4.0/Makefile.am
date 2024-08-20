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
# For using custom m4 macros
ACLOCAL_AMFLAGS = -I m4

print-%  : ; @echo $* = $($*)

dist_pkgdata_DATA = LICENSE INSTALL.txt README
# For checking Makefiles in folders
SUBDIRS = src/ etc/ doc/ share/crsim share/crsim/test

# For including extra files in distcheck
EXTRA_DIST = autogen.sh 

