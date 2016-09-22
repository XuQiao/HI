#/bin/csh
@ i = 1
@ j = 1
set typearr = (ptIn25_4 ptIn ptcoarser ptfiner centIn)
while( $i <= 4 ) 
  @ j = 1
  set type = $typearr[$i]
  while( $j <= 2 ) 
    @ ii = $j - 1
root -l << EOF
.x calccorr.C($ii,"$type")
EOF
    @ j++
  end
root -l << EOF
.x plotcnbbccsouth.C("$type")
EOF
@ i++
end


