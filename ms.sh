#!/bin/bash

if [ $# -ne 6 ]
then
	echo "$0 nsex ncontigs xyfrac pi l time"
	echo "nsex = number of individuals per sex"
	echo "ncontigs = total number of contigs"
	echo "xyfrac = fraction of sex-linked contigs"
	echo "pi = polymorphism"
	echo "l = sequence length"
	echo "time = time since the split of the sex chromosomes, expressed in units of 4*N_e"
	exit 1
fi

ncontigs=$2
nautocontigs=`bc <<< "$ncontigs-($2*$3/1)"`
time=$6
nsex=$1
theta=`bc <<< "scale=6; $4*$5/(1-$4)"`

nchrom=$(($nsex*4))
nX=$(($nsex*3))

echo "nsex=$nsex" "ncontigs=$ncontigs" "xyfrac=$3" "pi=$4" "length=$5" "time=$6" 
for ((i=1; i <= ncontigs ; i++))
do
if [[ $i -le $nautocontigs ]]
then
echo ">auto_contig${i}"
/home/jos/Software/msdir/ms $nchrom 1 -t $theta | tail -n +7
else
echo ">sl_contig${i}"
/home/jos/Software/msdir/ms $nchrom 1 -t $theta -I 2 $nsex $nX -n 1 0.25 -n 2 0.75 -ej $time 1 2 -eN $time 1 | tail -n +7
fi
done
