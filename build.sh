#!/bin/bash

#No queremos copiar src
subdirs=$(find . -type d ! -name "src")

for dir in $subdirs
do
    if [ -f "$dir/main.f95" ]
    then
        #Almacenar los copiados
        copied=()
        for file in src/*
        do
            dest="$dir/${file##*/}"
            if [ ! -e "$dest" ]
            then
                cp -r "$file" "$dir"
                copied+=("$dest")
            fi
        done
        
        cd "$dir"
        makef
        
        #Quitar solo los copiados
        for file in "${copied[@]}"
        do
            rm -rf "$file"
        done

        cd ..
    fi
done














