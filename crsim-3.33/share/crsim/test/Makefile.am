# ---------------------------------------------------------------------------
# COPYRIGHT    : Aleksandra Tatarevic
# ----------------------------------------------------------------------------
# PROJECT      : DORSY
# FILE         : Makefile.am
# VERSION      : $Revision: $
# DATE         : $Date: $
# TYPE         : Automake file
# IDENTIFIER   : $Id: $
# ----------------------------------------------------------------------------
#
# This file contains the basic configuration for the generation of a Makefile
#
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
#
# ----------------------------------------------------------------------------

print-%  : ; @echo $* = $($*)

TESTS = run_test1.sh run_test2.sh run_test3.sh run_test4.sh

EXTRA_DIST = run_test1.sh.in \
             run_test2.sh.in \
             run_test3.sh.in \
             run_test4.sh.in 


auxdir = $(datadir)/crsim/test
ddir = .$(pwd)

edit = sed -e 's|@nccmp[@]|$(NCCMP)|g'

run_test1.sh: run_test1.sh.in Makefile
	$(edit) '$(srcdir)/$@.in' > $@
	chmod 755 run_test1.sh
run_test2.sh: run_test2.sh.in Makefile
	 $(edit) '$(srcdir)/$@.in' > $@
	chmod 755 run_test2.sh
run_test3.sh: run_test3.sh.in Makefile
	$(edit) '$(srcdir)/$@.in' > $@
	chmod 755 run_test3.sh
run_test4.sh: run_test4.sh.in Makefile
	$(edit) '$(srcdir)/$@.in' > $@
	chmod 755 run_test4.sh


nobase_dist_aux_DATA =	run_test1.sh \
                        run_test2.sh \
                        run_test3.sh \
                        run_test4.sh 


install-data-hook:
	chmod a+rx $(auxdir)/run_test1.sh
	chmod a+rx $(auxdir)/run_test2.sh
	chmod a+rx $(auxdir)/run_test3.sh
	chmod a+rx $(auxdir)/run_test4.sh

installcheck-local:
	cd $(auxdir) && sh run_test1.sh
	cd $(auxdir) && sh run_test2.sh
	cd $(auxdir) && sh run_test3.sh
	cd $(auxdir) && sh run_test4.sh

CLEANFILES = run_test1.sh run_test2.sh run_test3.sh run_test4.sh

