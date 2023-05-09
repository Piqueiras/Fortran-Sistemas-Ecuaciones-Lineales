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
            cp -r "$file" "$dir"
            copied+=("${file##*/}")
        done
        
        cd "$dir"
        makef

        echo "${copied[@]}"
        
        #Quitar solo los copiados
        for file in "${copied[@]}"
        do
            rm -rf "$file"
        done

        cd ..
    fi
done














