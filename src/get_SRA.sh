#!/bin/bash

while read i
do
	fasterq-dump $i --split-files --skip-technical --threads 8 --progress
done < $1
