# Get all the files in the current directory with "upload" in the name
all_files=`ls | grep 'upload'`

# Find all prefix of the files
for file in $all_files
do
  prefix=`echo $file | cut -d'-' -f1`
  # If the prefix is not in the list, add it
  if [[ ! " ${prefixes[@]} " =~ " ${prefix} " ]]
  then
    echo "Adding $prefix"
    prefixes+=($prefix)
  fi
done


# For each prefix, get the last file with that prefix and overwrite the original file (with only the prefix)
for prefix in ${prefixes[@]}
do
  lastFile=`(ls -Art | grep $prefix | tail -n 1)`
  if [ -e $lastFile ]
  then
    echo "Moving $lastFile to $prefix"
    mv $lastFile $prefix
  else
    echo "No file with prefix $prefix"
  fi
done

# Remove all files with "upload" in the name
echo "Removing all files with upload in the name"
# List all files with "upload" in the name that are not directories and remove them
rm *upload*

echo " "
echo "-------------------------------------"
echo "-------------------------------------"
echo "-------------------------------------"
echo " "

# Run the main.py file with the config.yaml file
python3 main.py config.yml
