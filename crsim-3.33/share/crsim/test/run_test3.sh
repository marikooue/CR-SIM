#! /bin/bash
# ----------------------------------------------------------------------------
#
# ---------------------------------------------------------------------------
# COPYRIGHT    : Aleksandra Tatarevic
# ----------------------------------------------------------------------------
# PROJECT      : CRSIM
# FILE         : run_test3.sh.(in)
# VERSION      : $Revision: $
# DATE         : $Date: $
# TYPE         : Automake file
# IDENTIFIER   : $Id: $
# ----------------------------------------------------------------------------
#
# This file contains the test script running an example with RAMS microph.
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
# ----------------------------------------------------------------------------
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# given tolerance:
nccmp_tol=1.0e-4
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
echo "========================================================================================"
echo " TEST 3 "
echo "========================================================================================"
echo ""
#
testDir=`"pwd"`
binDir="${testDir}/../../../bin"

echo "binDir= " ${binDir}
echo "testDir= " ${testDir}

#  define test directories
inputDir=${testDir}/inp/rams_MP40/
logDir=${testDir}/logs
outDir=${testDir}/out

echo " "
echo "inputDir= "${inputDir}
echo "outDir=" ${outDir}
#
# Log File 
log_OutFile="$logDir/test3.log"

# Names of the input PARAMETERS file
InputParamFile="${inputDir}PARAMETERS"
InputFile1="${inputDir}a-A-2013-06-19-210000-g3.h5"
InputFile2="${inputDir}a-A-2013-06-19-210000-head.txt"
OutputFile="$outDir/rams_MP40_a-A-2013-06-19-210000-g3.nc"

echo "InputParamFile= "${InputParamFile}
echo " "

# remove any previously generated files
rm -fr ${outDir}/*.nc
#
# 
astatus=0 # status for all tests, =0 all tests passed, /=0 otherwise

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                 

nt_status=1  # nominal test case status,  =0 test passed, /=0 otherwise
echo '-->>     Running CRSIM ...'
echo ""
echo ${binDir}/crsim ${InputParamFile} ${InputFile1} ${OutputFile} &> ${log_OutFile}
${binDir}/crsim ${InputParamFile} ${InputFile1} ${OutputFile} &> ${log_OutFile} 
ExitCode=$?

 if [ "${ExitCode}" -eq 0 ] ; then 
     nt_status=0 
     echo '-->> CRSIM successfully executed'
 else
     echo '-->> CRSIM cannot be executed'
 fi

 echo "========================================================================================"
 echo ""
 echo ""
 echo '-->> dumping logs/test3.log'
 echo ""
 cat ${log_OutFile}
 echo ""
 echo ""
 echo ""
 echo ""
 echo ""


# save status if the test didn't pass
 if [ ${nt_status} -ne 0 ] ; then
    astatus=1
 fi
#  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#
 nccmp_status=0
# check if nccmp exists
if ! [ -x "$(command -v nccmp)" ]; then
  echo '-->> Error: nccmp is not installed. Tests 3-0 to 3-5 cannot be executed.' >&2
  nccmp_status=1
fi




if [ ${nccmp_status} -eq 0 ] ; then # if nccmp exists 


    if [ ${nt_status} -eq 0 ] ; then # if CRSIM succesefully executed

       status1=1 # =0 if the main output is succesefully compared with reference data set, /=0 otherwise
       echo "========================================================================================"
       echo "========================================================================================"
       echo " TEST 3-0  THE MAIN OUTPUT : CHECKING DIFFERENCES for tolerance '$nccmp_tol'(absolute units) "
       echo "========================================================================================"
       echo ""
       outputDiff="${logDir}/test3-0.log"
       file1=${testDir}/out/rams_MP40_a-A-2013-06-19-210000-g3.nc
       file2=${testDir}/out_ref/rams_MP40_a-A-2013-06-19-210000-g3.nc
       rm -f ${outputDiff}
    
       if (nccmp -F -d --tolerance=$nccmp_tol --statistics -f $file1  $file2  2>&1 ) > $outputDiff ; then
          ExitCode=$?
          status1=${ExitCode}
    
          numberOfDifferences=`grep 'DIFFER'  $outputDiff | wc -l`
          echo '-->> check of differences for' $file1 ' and ' $file2
          echo '-->> number of differences' $numberOfDifferences
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-0 Sucessful'
          echo '-->>'
       else
          echo '-->>'
          echo '-->> nccmp fails for test test3-0; dumping logs/test3-0.log'
          cat ${outputDiff}
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-0 Fails' 
          echo '-->>'
       fi

        # save status if the test didn't pass
        if [ ${status1} -ne 0 ] ; then
            astatus=1
        fi
        #  
    fi
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 
    echo ""
    echo ""
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


    if [ ${nt_status} -eq 0 ] ; then # if CRSIM succesefully executed
    
       status1=1 # =0 if CRSIM succesefully compared with reference data set, /=0 otherwise
       echo "========================================================================================"
       echo "========================================================================================"
       echo " TEST Test3-1  Hydrometeor=cloud : CHECKING DIFFERENCES for tolerance '$nccmp_tol' (absolute units) "
       echo "========================================================================================"
       echo ""
       outputDiff="${logDir}/test3-1.log"
       file1=${testDir}/out/rams_MP40_*cloud.nc
       file2=${testDir}/out_ref/rams_MP40_*cloud.nc
       rm -f ${outputDiff}

       if (nccmp -F -d --tolerance=$nccmp_tol --statistics -f $file1  $file2  2>&1 ) > $outputDiff ; then
          ExitCode=$?
          status1=$ExitCode 
    
          numberOfDifferences=`grep 'DIFFER'  $outputDiff | wc -l`
          echo '-->> check of differences for' $file1 ' and ' $file2
          echo '-->> number of differences' $numberOfDifferences
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-1 Sucessful'
          echo '-->>'
       else
          echo '-->>'
          echo '-->> nccmp fails for test test3-1; dumping logs/test3-1.log'
          cat ${outputDiff}
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-1 Fails' 
          echo '-->>'
       fi

       # save status if the test didn't pass
        if [ ${status1} -ne 0 ] ; then
            astatus=1
        fi
        #  
    fi
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    
    echo ""
    echo ""
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    
    if [ ${nt_status} -eq 0 ] ; then # if CRSIM succesefully executed

       status1=1 # =0 if CRSIM succesefully compared with reference data set, /=0 otherwise
       echo "========================================================================================"
       echo "========================================================================================"
       echo " TEST Test3-2  Hydrometeor=ice : CHECKING DIFFERENCES for tolerance '$nccmp_tol' (absolute units) "
       echo "========================================================================================"
       echo ""
       outputDiff="${logDir}/test3-2.log"
       file1=${testDir}/out/rams_MP40_*ice.nc
       file2=${testDir}/out_ref/rams_MP40_*ice.nc
       rm -f ${outputDiff}

       if (nccmp -F -d --tolerance=$nccmp_tol --statistics -f $file1  $file2  2>&1 ) > $outputDiff ; then
          ExitCode=$?
          status1=$ExitCode

          numberOfDifferences=`grep 'DIFFER'  $outputDiff | wc -l`
          echo '-->> check of differences for' $file1 ' and ' $file2
          echo '-->> number of differences' $numberOfDifferences
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-2 Sucessful'
          echo '-->>'
       else
          echo '-->>'
          echo '-->> nccmp fails for test test3-2; dumping logs/test3-2.log'
          cat ${outputDiff}
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-2 Fails' 
          echo '-->>'
       fi

       # save status if the test didn't pass
        if [ ${status1} -ne 0 ] ; then
            astatus=1
        fi
        #  
    fi

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    echo ""
    echo ""
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    if [ ${nt_status} -eq 0 ] ; then # if CRSIM succesefully executed

       status1=1 # =0 if CRSIM succesefully compared with reference data set, /=0 otherwise
       echo "========================================================================================"
       echo "========================================================================================"
       echo " TEST Test3-3  Hydrometeor=rain : CHECKING DIFFERENCES for tolerance '$nccmp_tol' (absolute units) "
       echo "========================================================================================"
       echo ""
       outputDiff="${logDir}/test3-3.log"
       file1=${testDir}/out/rams_MP40_*rain.nc
       file2=${testDir}/out_ref/rams_MP40_*rain.nc
       rm -f ${outputDiff}

       if (nccmp -F -d --tolerance=$nccmp_tol --statistics -f $file1  $file2  2>&1 ) > $outputDiff ; then
          ExitCode=$?
          status1=$ExitCode

          numberOfDifferences=`grep 'DIFFER'  $outputDiff | wc -l`
          echo '-->> check of differences for' $file1 ' and ' $file2
          echo '-->> number of differences' $numberOfDifferences
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-3 Sucessful'
          echo '-->>'
       else
          echo '-->>'
          echo '-->> nccmp fails for test test3-3; dumping logs/test3-3.log'
          cat ${outputDiff}
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-3 Fails' 
          echo '-->>'
       fi

       # save status if the test didn't pass
        if [ ${status1} -ne 0 ] ; then
            astatus=1
        fi
        #  
    fi

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    echo ""
    echo ""
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    if [ ${nt_status} -eq 0 ] ; then # if CRSIM succesefully executed

       status1=1 # =0 if CRSIM succesefully compared with reference data set, /=0 otherwise
       echo "========================================================================================"
       echo "========================================================================================"
       echo " TEST Test3-4  Hydrometeor=snow : CHECKING DIFFERENCES for tolerance '$nccmp_tol' (absolute units) "
       echo "========================================================================================"
       echo ""
       outputDiff="${logDir}/test3-4.log"
       file1=${testDir}/out/rams_MP40_*snow.nc
       file2=${testDir}/out_ref/rams_MP40_*snow.nc
       rm -f ${outputDiff}

       if (nccmp -F -d --tolerance=$nccmp_tol --statistics -f $file1  $file2  2>&1 ) > $outputDiff ; then
          ExitCode=$?
          status1=$ExitCode

          numberOfDifferences=`grep 'DIFFER'  $outputDiff | wc -l`
          echo '-->> check of differences for' $file1 ' and ' $file2
          echo '-->> number of differences' $numberOfDifferences
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-4 Sucessful'
          echo '-->>'
       else
          echo '-->>'
          echo '-->> nccmp fails for test test3-4; dumping logs/test3-4.log'
          cat ${outputDiff}
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-4 Fails' 
          echo '-->>'
       fi

       # save status if the test didn't pass
        if [ ${status1} -ne 0 ] ; then
            astatus=1
        fi
        #  
    fi

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    echo ""
    echo ""
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    if [ ${nt_status} -eq 0 ] ; then # if CRSIM succesefully executed

       status1=1 # =0 if CRSIM succesefully compared with reference data set, /=0 otherwise
       echo "========================================================================================"
       echo "========================================================================================"
       echo " TEST Test3-5  Hydrometeor=graupel : CHECKING DIFFERENCES for tolerance '$nccmp_tol' (absolute units) "
       echo "========================================================================================"
       echo ""
       outputDiff="${logDir}/test3-5.log"
       file1=${testDir}/out/rams_MP40_*graupel.nc
       file2=${testDir}/out_ref/rams_MP40_*graupel.nc
       rm -f ${outputDiff}

       if (nccmp -F -d --tolerance=$nccmp_tol --statistics -f $file1  $file2  2>&1 ) > $outputDiff ; then
          ExitCode=$?
          status1=$ExitCode

          numberOfDifferences=`grep 'DIFFER'  $outputDiff | wc -l`
          echo '-->> check of differences for' $file1 ' and ' $file2
          echo '-->> number of differences' $numberOfDifferences
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-5 Sucessful'
          echo '-->>'
       else
          echo '-->>'
          echo '-->> nccmp fails for test test3-5; dumping logs/test3-5.log'
          cat ${outputDiff}
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-5 Fails' 
          echo '-->>'
       fi

       # save status if the test didn't pass
        if [ ${status1} -ne 0 ] ; then
            astatus=1
        fi
        #  
    fi

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    echo ""
    echo ""
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    if [ ${nt_status} -eq 0 ] ; then # if CRSIM succesefully executed

       status1=1 # =0 if CRSIM succesefully compared with reference data set, /=0 otherwise
       echo "========================================================================================"
       echo "========================================================================================"
       echo " TEST Test3-6  Hydrometeor=hail : CHECKING DIFFERENCES for tolerance '$nccmp_tol' (absolute units) "
       echo "========================================================================================"
       echo ""
       outputDiff="${logDir}/test3-6.log"
       file1=${testDir}/out/rams_MP40_*hail.nc
       file2=${testDir}/out_ref/rams_MP40_*hail.nc
       rm -f ${outputDiff}

       if (nccmp -F -d --tolerance=$nccmp_tol --statistics -f $file1  $file2  2>&1 ) > $outputDiff ; then
          ExitCode=$?
          status1=$ExitCode

          numberOfDifferences=`grep 'DIFFER'  $outputDiff | wc -l`
          echo '-->> check of differences for' $file1 ' and ' $file2
          echo '-->> number of differences' $numberOfDifferences
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-6 Sucessful'
          echo '-->>'
       else
          echo '-->>'
          echo '-->> nccmp fails for test test3-6; dumping logs/test3-6.log'
          cat ${outputDiff}
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-6 Fails' 
          echo '-->>'
       fi

       # save status if the test didn't pass
        if [ ${status1} -ne 0 ] ; then
            astatus=1
        fi
        #  
    fi

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    echo ""
    echo ""
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    if [ ${nt_status} -eq 0 ] ; then # if CRSIM succesefully executed

       status1=1 # =0 if CRSIM succesefully compared with reference data set, /=0 otherwise
       echo "========================================================================================"
       echo "========================================================================================"
       echo " TEST Test3-7  Hydrometeor=drizzle : CHECKING DIFFERENCES for tolerance '$nccmp_tol' (absolute units) "
       echo "========================================================================================"
       echo ""
       outputDiff="${logDir}/test3-7.log"
       file1=${testDir}/out/rams_MP40_*drizzle.nc
       file2=${testDir}/out_ref/rams_MP40_*drizzle.nc
       rm -f ${outputDiff}

       if (nccmp -F -d --tolerance=$nccmp_tol --statistics -f $file1  $file2  2>&1 ) > $outputDiff ; then
          ExitCode=$?
          status1=$ExitCode

          numberOfDifferences=`grep 'DIFFER'  $outputDiff | wc -l`
          echo '-->> check of differences for' $file1 ' and ' $file2
          echo '-->> number of differences' $numberOfDifferences
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-7 Sucessful'
          echo '-->>'
       else
          echo '-->>'
          echo '-->> nccmp fails for test test3-7; dumping logs/test3-7.log'
          cat ${outputDiff}
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-7 Fails' 
          echo '-->>'
       fi

       # save status if the test didn't pass
        if [ ${status1} -ne 0 ] ; then
            astatus=1
        fi
        #  
    fi

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    echo ""
    echo ""
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    if [ ${nt_status} -eq 0 ] ; then # if CRSIM succesefully executed

       status1=1 # =0 if CRSIM succesefully compared with reference data set, /=0 otherwise
       echo "========================================================================================"
       echo "========================================================================================"
       echo " TEST Test3-8  Hydrometeor=aggregate : CHECKING DIFFERENCES for tolerance '$nccmp_tol' (absolute units) "
       echo "========================================================================================"
       echo ""
       outputDiff="${logDir}/test3-8.log"
       file1=${testDir}/out/rams_MP40_*aggregate.nc
       file2=${testDir}/out_ref/rams_MP40_*aggregate.nc
       rm -f ${outputDiff}

       if (nccmp -F -d --tolerance=$nccmp_tol --statistics -f $file1  $file2  2>&1 ) > $outputDiff ; then
          ExitCode=$?
          status1=$ExitCode

          numberOfDifferences=`grep 'DIFFER'  $outputDiff | wc -l`
          echo '-->> check of differences for' $file1 ' and ' $file2
          echo '-->> number of differences' $numberOfDifferences
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-8 Sucessful'
          echo '-->>'
       else
          echo '-->>'
          echo '-->> nccmp fails for test test3-8; dumping logs/test3-8.log'
          cat ${outputDiff}
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test3-8 Fails' 
          echo '-->>'
       fi

       # save status if the test didn't pass
        if [ ${status1} -ne 0 ] ; then
            astatus=1
        fi
        #  
    fi

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    echo ""
    echo ""
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


#    if [ ${nt_status} -eq 0 ] ; then # if CRSIM succesefully executed
#  
#     echo " Tests test3-0, test3-1, test3-2, test3-3, test3-4,"
#     echo " test3-5, test3-6, test3-7 and test3-8 passed"    #  
#    fi
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    
    echo ""
    echo ""
else

          echo '-->> Tests test3 test3-0, test3-1, test3-2, test3-3, test3-4, '
          echo '-->> test3-5,  test3-6, test3-7 and test3-8 cannot be executed'   
          astatus=1
fi
    
 # status for all tests
 if [ "${astatus}" -eq 0 ] ; then
    echo  '-->> TEST test3 PASSED'
  else
    echo  '-->> TEST test3 FAILED' 
    sh -c "false" #; echo $?
  fi
    
