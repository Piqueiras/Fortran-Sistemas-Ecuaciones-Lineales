#!/bin/make -f

CC:=f95
OUTPUT:=ejecutable
SRCS:=$(wildcard *.f95)

$(OUTPUT):
	$(CC) -o $(OUTPUT).out $(SRCS)
	
cleanall: clean
	rm -f $(OUTPUT).out

clean:
	rm -f *.o *~
