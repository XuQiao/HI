#!/bin/csh
root -b << EOF
#include "RpPar.h"
.L Getdphi.C
int code = $1
int iangle1 = 0
int iangle2 = 0
bool usingCNTEP = 0
Getdphi(code, iangle1, iangle2, usingCNTEP)
EOF
echo "Done!"
#.L Getvn.C+
#Getvn(code, iangle1, iangle2, UseCNTEP)
