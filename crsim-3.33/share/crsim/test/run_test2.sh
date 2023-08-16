#! /bin/bash
# ----------------------------------------------------------------------------
#
# ---------------------------------------------------------------------------
# COPYRIGHT    : Aleksandra Tatarevic
# ----------------------------------------------------------------------------
# PROJECT      : CRSIM
# FILE         : run_test2.sh.(in)
# VERSION      : $Revision: $
# DATE         : $Date: $
# TYPE         : Automake file
# IDENTIFIER   : $Id: $
# ----------------------------------------------------------------------------
#
# This file contains the test script running an example with P3 microph.
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
echo " TEST 2 "
echo "========================================================================================"
echo ""
#
testDir=`"pwd"`
binDir="${testDir}/../../../bin"

echo "binDir= " ${binDir}
echo "testDir= " ${testDir}

#  define test directories
inputDir=${testDir}/inp/wrf_MP50/
logDir=${testDir}/logs
outDir=${testDir}/out

echo " "
echo "inputDir= "${inputDir}
echo "outDir=" ${outDir}
#
# Log File 
log_OutFile="$logDir/test2.log"

# Names of the input PARAMETERS file
InputParamFile="${testDir}/inp/wrf_MP50/PARAMETERS"
InputFile1="${testDir}/inp/wrf_MP50/P3_wrfout_d01_2011-05-20_10_00_00"
InputFile2="${testDir}/inp/wrf_MP50/ABGDP3_wrfout_d01_2011-05-20_10_00_00"
OutputFile="$outDir/wrf_p3_wrfout_d01_2015-06-27_20:00:00.nc"

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
 echo '-->> dumping logs/test2.log'
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
  echo '-->> Error: nccmp is not installed. Tests 2-0 to 2-6 cannot be executed.' >&2
  nccmp_status=1
fi




if [ ${nccmp_status} -eq 0 ] ; then # if nccmp exists 


    if [ ${nt_status} -eq 0 ] ; then # if CRSIM succesefully executed

       status1=1 # =0 if the main output is succesefully compared with reference data set, /=0 otherwise
       echo "========================================================================================"
       echo "========================================================================================"
       echo " TEST 2-0  THE MAIN OUTPUT : CHECKING DIFFERENCES for tolerance '$nccmp_tol'(absolute units) "
       echo "========================================================================================"
       echo ""
       outputDiff="${logDir}/test2-0.log"
       file1=${testDir}/out/wrf_p3_wrfout_d01_2015-06-27_20:00:00.nc
       file2=${testDir}/out_ref/wrf_p3_wrfout_d01_2015-06-27_20:00:00.nc
       rm -f ${outputDiff}
    
       if (nccmp -F -d --tolerance=$nccmp_tol --statistics -f $file1  $file2  2>&1 ) > $outputDiff ; then
          ExitCode=$?
          status1=${ExitCode}
    
          numberOfDifferences=`grep 'DIFFER'  $outputDiff | wc -l`
          echo '-->> check of differences for' $file1 ' and ' $file2
          echo '-->> number of differences' $numberOfDifferences
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test2-0 Sucessful'
          echo '-->>'
       else
          echo '-->>'
          echo '-->> nccmp fails for test test2-0; dumping logs/test2-0.log'
          cat ${outputDiff}
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test2-0 Fails' 
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
       echo " TEST Test2-1  Hydrometeor=cloud : CHECKING DIFFERENCES for tolerance '$nccmp_tol' (absolute units) "
       echo "========================================================================================"
       echo ""
       outputDiff="${logDir}/test2-1.log"
       file1=${testDir}/out/wrf_p3_*cloud.nc
       file2=${testDir}/out_ref/wrf_p3_*cloud.nc
       rm -f ${outputDiff}

       if (nccmp -F -d --tolerance=$nccmp_tol --statistics -f $file1  $file2  2>&1 ) > $outputDiff ; then
          ExitCode=$?
          status1=$ExitCode 
    
          numberOfDifferences=`grep 'DIFFER'  $outputDiff | wc -l`
          echo '-->> check of differences for' $file1 ' and ' $file2
          echo '-->> number of differences' $numberOfDifferences
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test2-1 Sucessful'
          echo '-->>'
       else
          echo '-->>'
          echo '-->> nccmp fails for test test2-1; dumping logs/test2-1.log'
          cat ${outputDiff}
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test2-1 Fails' 
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
       echo " TEST Test2-2  Hydrometeor=parimedice : CHECKING DIFFERENCES for tolerance '$nccmp_tol' (absolute units) "
       echo "========================================================================================"
       echo ""
       outputDiff="${logDir}/test2-2.log"
       file1=${testDir}/out/wrf_p3_*parimedice.nc
       file2=${testDir}/out_ref/wrf_p3_*parimedice.nc
       rm -f ${outputDiff}

       if (nccmp -F -d --tolerance=$nccmp_tol --statistics -f $file1  $file2  2>&1 ) > $outputDiff ; then
          ExitCode=$?
          status1=$ExitCode

          numberOfDifferences=`grep 'DIFFER'  $outputDiff | wc -l`
          echo '-->> check of differences for' $file1 ' and ' $file2
          echo '-->> number of differences' $numberOfDifferences
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test2-2 Sucessful'
          echo '-->>'
       else
          echo '-->>'
          echo '-->> nccmp fails for test test2-2; dumping logs/test2-2.log'
          cat ${outputDiff}
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test2-2 Fails' 
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
       echo " TEST Test2-3  Hydrometeor=rain : CHECKING DIFFERENCES for tolerance '$nccmp_tol' (absolute units) "
       echo "========================================================================================"
       echo ""
       outputDiff="${logDir}/test2-3.log"
       file1=${testDir}/out/wrf_p3_*rain.nc
       file2=${testDir}/out_ref/wrf_p3_*rain.nc
       rm -f ${outputDiff}

       if (nccmp -F -d --tolerance=$nccmp_tol --statistics -f $file1  $file2  2>&1 ) > $outputDiff ; then
          ExitCode=$?
          status1=$ExitCode

          numberOfDifferences=`grep 'DIFFER'  $outputDiff | wc -l`
          echo '-->> check of differences for' $file1 ' and ' $file2
          echo '-->> number of differences' $numberOfDifferences
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test2-3 Sucessful'
          echo '-->>'
       else
          echo '-->>'
          echo '-->> nccmp fails for test test2-3; dumping logs/test2-3.log'
          cat ${outputDiff}
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test2-3 Fails' 
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
       echo " TEST Test2-4  Hydrometeor=smallice : CHECKING DIFFERENCES for tolerance '$nccmp_tol' (absolute units) "
       echo "========================================================================================"
       echo ""
       outputDiff="${logDir}/test2-4.log"
       file1=${testDir}/out/wrf_p3_*smallice.nc
       file2=${testDir}/out_ref/wrf_p3_*smallice.nc
       rm -f ${outputDiff}

       if (nccmp -F -d --tolerance=$nccmp_tol --statistics -f $file1  $file2  2>&1 ) > $outputDiff ; then
          ExitCode=$?
          status1=$ExitCode

          numberOfDifferences=`grep 'DIFFER'  $outputDiff | wc -l`
          echo '-->> check of differences for' $file1 ' and ' $file2
          echo '-->> number of differences' $numberOfDifferences
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test2-4 Sucessful'
          echo '-->>'
       else
          echo '-->>'
          echo '-->> nccmp fails for test test2-4; dumping logs/test2-4.log'
          cat ${outputDiff}
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test2-4 Fails' 
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
       echo " TEST Test2-5  Hydrometeor=unrimedice : CHECKING DIFFERENCES for tolerance '$nccmp_tol' (absolute units) "
       echo "========================================================================================"
       echo ""
       outputDiff="${logDir}/test2-5.log"
       file1=${testDir}/out/wrf_p3_*unrimedice.nc
       file2=${testDir}/out_ref/wrf_p3_*unrimedice.nc
       rm -f ${outputDiff}

       if (nccmp -F -d --tolerance=$nccmp_tol --statistics -f $file1  $file2  2>&1 ) > $outputDiff ; then
          ExitCode=$?
          status1=$ExitCode

          numberOfDifferences=`grep 'DIFFER'  $outputDiff | wc -l`
          echo '-->> check of differences for' $file1 ' and ' $file2
          echo '-->> number of differences' $numberOfDifferences
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test2-5 Sucessful'
          echo '-->>'
       else
          echo '-->>'
          echo '-->> nccmp fails for test test2-5; dumping logs/test2-5.log'
          cat ${outputDiff}
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test2-5 Fails' 
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
       echo " TEST Test2-6  Hydrometeor=graupel : CHECKING DIFFERENCES for tolerance '$nccmp_tol' (absolute units) "
       echo "========================================================================================"
       echo ""
       outputDiff="${logDir}/test2-6.log"
       file1=${testDir}/out/wrf_p3_*graupel.nc
       file2=${testDir}/out_ref/wrf_p3_*graupel.nc
       rm -f ${outputDiff}

       if (nccmp -F -d --tolerance=$nccmp_tol --statistics -f $file1  $file2  2>&1 ) > $outputDiff ; then
          ExitCode=$?
          status1=$ExitCode

          numberOfDifferences=`grep 'DIFFER'  $outputDiff | wc -l`
          echo '-->> check of differences for' $file1 ' and ' $file2
          echo '-->> number of differences' $numberOfDifferences
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test2-6 Sucessful'
          echo '-->>'
       else
          echo '-->>'
          echo '-->> nccmp fails for test test2-6; dumping logs/test2.log'
          cat ${outputDiff}
          echo '-->> log file stored: '$outputDiff
          echo '-->> Test test2-6 Fails' 
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
#    echo " Tests test2-0, test2-1, test2-2, test2-3, test2-4,test2-5 and test2-6 passed"   
#    fi
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    
    echo ""
    echo ""
else

          echo '-->> Tests test2-0, test2-1, test2-2, test2-3, test2-4,test2-5 and test2-6 cannot be executed'   
          astatus=1
fi
    
 # status for all tests
 if [ "${astatus}" -eq 0 ] ; then
    echo  '-->> TEST test2 PASSED'
  else
    echo  '-->> TEST test2 FAILED' 
    sh -c "false" #; echo $?
  fi
    
