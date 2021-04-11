platforms=( osx-64 linux-32 linux-64 )



find $HOME/miniconda3/conda-bld/linux-64/ -name *.tar.bz2 | while read file
do
    echo $file
    #conda convert --platform all $file  -o $HOME/conda-bld/
    for platform in "${platforms[@]}"
    do
       conda convert --platform $platform $file  -o $HOME/conda-bld/
    done

done

# upload packages to conda
find $HOME/miniconda3/conda-bld/ -name *.tar.bz2 | while read file
do
    echo $file
    anaconda upload $file
done

echo "Building conda package done!"
