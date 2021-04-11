pkg='drukbam'
array=( 3.7 3.8 3.9 )

#echo "Buildinconda package ..."
#conda skeleton pypi $pkg
for i in "${array[@]}"
do
	conda-build -c bioconda -c conda-forge --python $i $pkg
done
