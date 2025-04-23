#!/bin/bash

if [ $# -eq 0 ]
then
	read -p "al catelea nr din sirul fibonacci vreti sa aflati: " x
fi

declare -i rez
rez=0

function fibonacci()
{
	local t1=$2
	local t2=$3

	if [ $1 -lt $x ]
	then
		let t1=t2
		let t2=rez
		rez=t1+t2
		fibonacci $(($1+1)) $t1 $t2
	fi
}

if [ $x -eq 1 ]
then
	rez=1
fi

if [ $x -eq 2 ]
then
	rez=1
fi

if [ $x -gt 2 ]
then
	let param=0
	let start1=1
	let start2=1
	fibonacci $param $start1 $start2
fi
echo $rez

# https://profs.info.uaic.ro/~vidrascu/SO/lectures/Linux/P4_shell2_web-ro.pdf
