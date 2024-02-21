#!/bin/bash

in_dir=$1
t=$2

mkfifo bg1
exec 5<>bg1
rm -f bg1

{
for ((i=1;i<=$t;i++));do
	echo;
done
} >&5

for dic in $in_dir/*ed
do
read -u5
{
cd $dic
codeml *codeml.ctl
echo >&5
}&
done <&5

wait

exec 5>&-

