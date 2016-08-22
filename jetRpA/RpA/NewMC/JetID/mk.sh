
arr=("chMax"  "chSum"  "neuMax"  "neuSum"  "phoMax"  "phoSum"  "chMaxpt"  "chSumpt"  "neuMaxpt"  "neuSumpt"  "phoMaxpt"  "phoSumpt" "SumSumpt" "SumSumrawpt" "neuMaxr" "chN" "neuN" "phoN")

cd pic
for var in "${arr[@]}"
do
  mkdir "${var}"
  # do something on $var
done
