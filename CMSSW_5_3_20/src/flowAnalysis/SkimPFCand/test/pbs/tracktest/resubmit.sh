#!/bin/bash

###################################################Initialization################################################

trig=("300260" "260220" "220185" "185150" "150120")
data=("Pbp" "pPb" "pPbReReco")
nfilesperjob=10
array1g=();	array1j=();	array1r=()
array2g=();	array2j=();	array2r=()
array3g=();	array3j=();	array3r=()
array4g=();	array4j=();	array4r=()
array5g=();	array5j=();	array5r=()
for i in {0..2};do
nfiles=`cat filelist${data[$i]}.dat | wc -l`
njobs[$i]=`echo "($nfiles/$nfilesperjob)+1" | bc`
for j in {0..4};do
mr[$j+$i*5]=0
mj[$j+$i*5]=0
mg[$j+$i*5]=0
done
done

#############################################Looking at qstat of jobs############################################

rm -r temp.txt
qstat | grep xuq7 > temp.txt
n=0
jtrig=()
jdata=()
jnumber=()
while read line
do
        job=$line
if [[ `echo $job | awk '{print $5}'` == "R" || `echo $job | awk '{print $5}'` == "Q" ]];then
jobname=`echo $job | awk '{print $2}'`
for j in {0..4};do
for i in {0..2};do
if [[ $jobname == *${trig[$j]}${data[$i]}* ]];then
jobtrigname=${trig[$j]}
jobdataset=${data[$i]}
fi
done
done
if [[ $jobname == *pPbReReco* ]];then
jobtrigname="300260"
jobdataset="pPbReReco"
fi
jobnumber=${jobname##*$jobdataset}
jtrig+=($jobtrigname)
jdata+=($jobdataset)
jnumber+=($jobnumber)
n=`echo "$n+1" | bc`
fi
done < "temp.txt"

#############################################Select and handle job names##########################################

for file in `ls joboutput`; do
HLTTrig=${file%%_*}
HLTTrig=${HLTTrig#$"job"}
dataset=${file%_*}
dataset=${dataset#*_}
number=${file##*_}
number=${number%$".out"}
trigname=$HLTTrig
nfiles=`cat filelist${dataset}.dat | wc -l`
njob=`echo "$nfiles/$nfilesperjob" | bc`

######################################################resubmit jobs#################################################

flag=0
if grep -q "basket" "joboutput/$file" || grep -q "FileOpenError" "joboutput/$file" || ! grep -q "Done" "joboutput/$file"; then
for iter in $(seq 0 $(($n-1)));do
if [[ $dataset == ${jdata[$iter]} && $trigname == ${jtrig[$iter]} && $number == ${jnumber[$iter]} ]];then
flag=1
break
fi
done
i=$number
begin=`echo "$i*$nfilesperjob" | bc`
end=`echo "($i+1)*$nfilesperjob" | bc`
if [[ $i == $njob ]];then
end=$nfiles
fi
if [[ $flag == 0 ]];then
echo -n ""
#echo now resubmit ${HLTTrig} ${dataset} $i
#qsub -v I=$i,BEGIN=$begin,END=$end,DATASET=$dataset,TRIGNAME=$trigname -N ${HLTTrig}${dataset}$i -z jobsub.pbs
fi

#####################################################Count jobs#########################################################

if [[ $flag == 0 ]];then
for i in {0..2}; do
for j in {0..4}; do
	if [[ $HLTTrig == ${trig[$j]} ]];then
	if [[ $dataset == ${data[$i]} ]]; then
	mr[$j+$i*5]=`echo "${mr[$j+$i*5]}+1" | bc`
	if [[ $j == 0 ]];then array1r+=($number)
	elif [[ $j == 1 ]]; then array2r+=($number)
	elif [[ $j == 2 ]]; then array3r+=($number)
	elif [[ $j == 3 ]]; then array4r+=($number)
	elif [[ $j == 4 ]]; then array5r+=($number)
	fi
fi
fi
done
done
fi

else
for i in {0..2}; do
for j in {0..4}; do
        if [[ $HLTTrig == ${trig[$j]} ]];then
        if [[ $dataset == ${data[$i]} ]]; then
        mg[$j+$i*5]=`echo "${mg[$j+$i*5]}+1" | bc`
        if [[ $j == 0 ]];then array1g+=($number)
        elif [[ $j == 1 ]]; then array2g+=($number)
        elif [[ $j == 2 ]]; then array3g+=($number)
        elif [[ $j == 3 ]]; then array4g+=($number)
        elif [[ $j == 4 ]]; then array5g+=($number)
        fi
fi
fi
done
done

fi
done

for i in {0..2}; do
for j in {0..4}; do
for iter in $(seq 0 $(($n-1)));do
        if [[ ${jtrig[$iter]} == ${trig[$j]} ]];then
        if [[ ${jdata[$iter]} == ${data[$i]} ]]; then
        mj[$j+$i*5]=`echo "${mj[$j+$i*5]}+1" | bc`
        if [[ $j == 0 ]];then array1j+=(${jnumber[$iter]})
        elif [[ $j == 1 ]]; then array2j+=(${jnumber[$iter]})
        elif [[ $j == 2 ]]; then array3j+=(${jnumber[$iter]})
        elif [[ $j == 3 ]]; then array4j+=(${jnumber[$iter]})
        elif [[ $j == 4 ]]; then array5j+=(${jnumber[$iter]})
        fi
fi
fi
done
done
done

##################################################Print out results#######################################################

echo
echo ${data[*]}
for j in {0..0}; do
	echo -e "\t\t\t\\e[4m${trig[$j]}\e[24m"
	echo -ne "\e[34mall: ${njobs[*]}"
	echo -e "\e[39m"
	echo -ne "\e[33mNow running or queueing:"
	for i in {0..2};do
		echo -n " ${mj[$j+$i*5]}"
	done
	echo
	k=-1
	for i in {0..2};do
		if [[ $i ==  1 ]] ;then echo -ne "\e[33m${data[$i]}: ["
		else	echo -ne "\e[39m${data[$i]}: [";	fi
		for k in $(seq $(($k+1)) $(($k+${mj[$j+$i*5]})));do
			echo -n ${array1j[$k]}", "
		done
	if [[ ${mj[$j+$i*5]} != 0 ]];then echo -ne "\b\b";fi
	echo "]"
	done

        echo -ne "\e[31mNeed to resubmit:"
        for i in {0..2};do
                echo -n " ${mr[$j+$i*5]}"
        done
        echo
        k=-1
        for i in {0..2};do
                if [[ $i ==  1 ]] ;then echo -ne "\e[31m${data[$i]}: ["
                else    echo -ne "\e[39m${data[$i]}: [";        fi
                for k in $(seq $(($k+1)) $(($k+${mr[$j+$i*5]})));do
                        echo -n ${array1r[$k]}", "
                done
        if [[ ${mr[$j+$i*5]} != 0 ]];then echo -ne "\b\b";fi
        echo "]"
        done

	echo -ne "\e[32mGood: "
	for i in {0..2};do
		echo -n " ${mg[$j+$i*5]}"
	done
	echo -e "\e[39m"
	k=-1
	for i in {0..2};do
		if [[ $i ==  1 ]] ;then echo -ne "\e[32m${data[$i]}: ["
		else	echo -ne "\e[39m${data[$i]}: [";	fi
		for k in $(seq $(($k+1)) $(($k+${mg[$j+$i*5]})));do
			echo -n ${array1g[$k]}", "
		done
	if [[ ${mg[$j+$i*5]} != 0 ]];then echo -ne "\b\b"; fi
	echo "]"
	done
done
