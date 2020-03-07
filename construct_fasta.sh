#!/bin/bash

file=$1
cancer=$2
path=$3
cd ${path}
for line in $(cat "${file}" | cut -f5); do 
	echo -n ">" >> ${file}.fa
	echo "${line}" >> ${file}.fa
	echo "${line}" >> ${file}.fa
done
	
