#!/bin/csh
root -b << EOF
#include "RpPar.h"
.L Getvn.C
int code = $1
int iangle1 = 0
int iangle2 = 0
bool usingCNTEP = 0
Getvn(code, iangle1, iangle2, usingCNTEP)
EOF
echo "Done!"
#.L Getdphi.C
#.L Getvn.C+
#Getvn(code, iangle1, iangle2, UseCNTEP)
