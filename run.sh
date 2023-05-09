#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "Pasa el fichero por linea de comandos"
    exit 1
fi

if [ ! -f "$1" ]
  then
    echo "A ser posible un fichero que exista"
    exit 1
fi

# Vaciar fichero
> output

#Iterar ejecutables
for dir in */
do
  if [ "$dir" != "src" ] && [ -f "$dir/ejecutable.out" ]   #Comprobar que existe el ejecutable
  then
    "$dir/ejecutable.out" < "$1" >> output
  fi
done