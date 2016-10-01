#/bin/bash
#i=1
#j=1
typearr=(ptIn ptfiner centIn per per_ptIn)
for i in {0..4};do
  #j=1
  type=${typearr[$i]}
  for j in {0..2};do
      ii=$((j))
      echo $ii, $type
root -l << EOF
.x calccorr.C($ii,"$type")
EOF
#j=$((j+1))
done
#root -l << EOF
#.x plotcnbbccsouth.C("$type")
#EOF
#i=$((i+1))
done


