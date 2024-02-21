#!/usr/bash -l

in_dir=$1

for gene in $in_dir/*
do
cp /home/chenshanshan/paml4.9j/bin/codeml $gene
done
	
