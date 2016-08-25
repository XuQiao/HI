#/bin/csh
@ i = 1
set typearr = (ptIn25_4 ptIn ptcoarser ptfiner centIn ptccentc)
while( $i <= 6 ) 
  set type = $typearr[$i]
echo " Do "$type
root -l << EOF
.x calccorr_sandv.C("$type")
EOF
echo "calculation done!"
root -l << EOF
.x plotcnsandv.C("$type")
EOF
@ i++
end
