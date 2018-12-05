#!/bin/bash

i=0

while read line
do
if [[ ${line:0:2} == "#p" ]]
then
	f[0]=`echo $line | awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "logL_autosomal") print i } }'`
	f[1]=`echo $line | awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "logL_haploid") print i } }'`
	f[2]=`echo $line | awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "logL_paralog") print i } }'`
	f[3]=`echo $line | awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "logL_xhemizygote") print i } }'`
	f[4]=`echo $line | awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "logL_xy") print i } }'`
	f[5]=`echo $line | awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "logL_zhemizygote") print i } }'`
	f[6]=`echo $line | awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "logL_zw") print i } }'`
#	echo -n "fields: "
#	for j in 0 1 2 3 4 5 6 7
#	do
#	echo -n "${f[$j]} "
#	done
#	echo
fi
if [[ ${line:0:1} == ">" ]] 
then
	if [[ $i -gt 0 ]]
	then
		echo -n "$contigname " 	
	sumL=0
	for j in 0 1 2 3 4 5 6 7

	do
	if [[ ${f[$j]} ]]
	then
		sumL=`echo ${L[$j]} | awk -v s="$sumL" '{ result = $1 + s ; print result}'`
	fi
	echo -n "${L[$j]} "
	done
	echo -n "$sumL "
	for j in 0 1 2 3 4 5 6 7

	do
		echo ${L[$j]} | awk -v s="$sumL" '{ result = $1 / s ; printf "%f ", result }'
	done
	echo
	fi
	contigname=`echo $line | awk '{print $1}'`
	i=$((i+1))
#	echo $contigname $i
	for j in 0 1 2 3 4 5 6 7
	do
	L[$j]=0
	done
elif [[ ${line:0:1} != "#" && $i -gt 0 ]]
then
#	echo -n "start: $L1"
#	echo $line | awk -v L1="$L1" -v f1="$f1" '{newl = L1 + exp($f1) ; print L1, f1, $f1, newl}'
	for j in 0 1 2 3 4 5 6 7
	do
#	echo $j ${f[$j]}
	if [[ ${f[$j]} ]]
	then	
#	echo "L: $j ${L[$j]}"
	L[$j]=`echo $line | awk -v L1="${L[$j]}" -v f1="${f[$j]}" '{newl = L1 + exp($f1) ; print newl}'`
#	echo "L: $j ${L[$j]}"
	fi
	done
#	echo -e "\tafter: $L1"
fi
done < $1

	if [[ $i -gt 0 ]]
	then
		echo -n "$contigname " 	
	sumL=0
	for j in 0 1 2 3 4 5 6 7

	do
	if [[ ${f[$j]} ]]
	then
		sumL=`echo ${L[$j]} | awk -v s="$sumL" '{ result = $1 + s ; print result}'`
	fi
	echo -n "${L[$j]} "
	done
	echo -n "$sumL "
	for j in 0 1 2 3 4 5 6 7

	do
		echo ${L[$j]} | awk -v s="$sumL" '{ result = $1 / s ; printf "%f ", result }'
	done
	echo
	fi

