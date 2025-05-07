#!/usr/bin/bash

if [ $2 -lt $1 ]; then

	echo "Error, $2 is smaller than $1"
	exit 1
fi

num=$1

while [ $num -le $2 ]; do 

	echo -n "${num} "
	((num++))
done
