#!/bin/bash

set -x

# Parse argument list
usage=0
while getopts ":s:i:o:d:y:U" opt; do
   case $opt in
     s    ) src_dir=$OPTARG ;;
     i    ) in_data_dir=$OPTARG ;;
     o    ) out_dir=$OPTARG ;;
     d    ) descriptor=$OPTARG ;;
     y    ) years=$OPTARG ;;
     \?   ) usage=1 ;;
   esac
done
shift $(($OPTIND - 1))

if [ "$src_dir" = "" ]; then
   echo "ERROR: no argument given for (python) source scripts directory."
   src_dir=$NBROOT/mom6/tools/analysis
else
  if [ ! -d "$src_dir" ]; then
     echo "ERROR: argument given for (python) source scripts directory is not a directory: "$src_dir
     usage=1
  fi
fi

if [ "$in_data_dir" = "" ]; then
   echo "ERROR: no argument given for input directory."
   usage=1
else
  if [ ! -d "$in_data_dir" ]; then
     echo "ERROR: argument given for input directory is not a directory: "$in_data_dir
     usage=1
  fi
fi

if [ "$out_dir" = "" ]; then
   echo "ERROR: no argument given for output directory."
   usage=1
fi

if [ "$descriptor" = "" ]; then
   echo "ERROR: no argument given for descriptor."
   usage=1
fi

if [ "$yr1" = "" ] || [ "$yr2" = "" ] ; then
   if [ "$years" = "" ]; then
      echo "ERROR: no argument given for years."
      usage=1
   else
      yrs=(`echo $years | sed -e "s/,/ /g"`)
      if [ ${#yrs[@]} == 2 ]; then
         yr1=${yrs[0]}
         yr2=${yrs[1]}
         databegyr=${yrs[0]}
         dataendyr=${yrs[1]}
         datachunk=`echo ${dataendyr}-${databegyr}+1 | bc`
      elif [ ${#yrs[@]} == 5 ]; then
         yr1=${yrs[0]}
         yr2=${yrs[1]}
         databegyr=${yrs[2]}
         dataendyr=${yrs[3]}
         datachunk=${yrs[4]}
      else
         echo "ERROR: invalid entry for years: "$years
         usage=1
      fi
   fi
fi

# Print usage message and exit in case of error
if [ $usage -ne 0 ]; then
   set cmdname = `basename $0`
   echo
   echo "USAGE:  $cmdname [-U] [-s sdir] [-i idir] [-o odir] [-d desc] [-y years]  files"
   echo
   echo "         sdir  = path of directory that contains the python scripts"
   echo "         idir  = path of directory containing input data"
   echo "                 (for example: -i /archive/user/MyExp/pp/atmos/ts/monthly/20yr)"
   echo "         odir  = path of directory for postscript figures"
   echo "         desc  = descriptor for input data added to figures"
   echo "         years = year range for observed data"
   echo "                 (for example: -y 1981,2000)"
   echo
   echo
   exit 1
fi


# Prepare workdir
workdir=$TMPDIR/nnz_ocean_ts_annual
rm -rf $workdir
mkdir -p $workdir
cd $workdir

# Component name
compName=`echo $in_data_dir | perl -e 'while(<>){print /pp.*\/(.*?)\/ts\//}'`

# List of variables
varlist='temp salt'

# Determine available input variables and generate list of input files
cd ${in_data_dir}
in_data_var=''
filelist='../../../ocean_annual.static.nc'
for var in ${varlist}
do

  file_present=0
  for file in `ls -1 ${compName}.????-????.${var}.nc 2> /dev/null`
  do

    # First year
    yra=${file##${compName}.}
    yra=${yra%%-????.${var}.nc}

    # Last year
    yrb=${file##${compName}.????-}
    yrb=${yrb%%.${var}.nc}

    # Does the time span of current file overlap with desired period yr1,yr2?
    # Note: use 10# to force decimal comparison, otherwise numbers with leading 
    # zeros will be interpreted as octal numbers.
    if (( 10#${yra} <= 10#${yr2} && 10#${yrb} >= 10#${yr1} ))
    then
      filelist=${filelist}' '${file}
      file_present=1
    fi

  done

  if (( file_present == 1 ))
  then
    in_data_var=${in_data_var}' '${var}
  fi

done

# dmget and copy files
cd ${in_data_dir}
dmget ${filelist}
gcp ${filelist} $workdir
cd $workdir

# Generate include file common to all scripts
cat > files.txt << EOF

   descriptor = "${descriptor}"
   years = "${yr1}-${yr2}"

EOF

# Generate plots using python
cp ${src_dir}/MOM6_ts_annual.py .
python MOM6_ts_annual.py

# Compress ps files
#gzip *.ps

# Copy output files to their destination
dest_dir=${out_dir}/${compName}_${yr1}_${yr2}/MOM6_test.ts/
gcp -cd *.ps* ${dest_dir}

exit 0
