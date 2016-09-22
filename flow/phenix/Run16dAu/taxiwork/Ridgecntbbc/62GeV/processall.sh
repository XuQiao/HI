#/bin/bash
i=1
j=1
typearr=(ptIn25_4 ptIn ptcoarser ptfiner centIn ptccentc)
for i in {0..5};do
  j=1
  type=${typearr[$i]}
  for j in {0..1};do
      ii=$((j-1))
root -l << EOF
.x calccorr.C($ii,"$type")
EOF
j=$((j+1))
done
root -l << EOF
.x plotcnbbccsouth.C("$type")
EOF
i=$((i+1))
done


