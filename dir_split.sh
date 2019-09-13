#!/bin/bash

#dir_size=2000
dir_name="vm"
#n=$((`find . -maxdepth 1 -type f | wc -l`/$dir_size+1))
n=3
dir_size=$((`find . -maxdepth 1 -type f | wc -l`/3))

echo "creating "$n "directories of "$dir_size "elements"

for i in `seq 1 $n`;
do
    mkdir -p "$dir_name$i";
    find . -maxdepth 1 -type f | head -n $dir_size | xargs -i mv "{}" "$dir_name$i"
done

