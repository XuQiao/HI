#!/bin/csh
root -b << EOF
#include "RpPar.h"
.L Getdphi.C+
int code = $1
int iangle1 = 0
int iangle2 = 0
int ihar = code % nhar
int icent = code / nhar
bool usingCNTEP = 0
Getdphi(icent, ihar, iangle1, iangle2, usingCNTEP)
EOF
echo "Done!"
#.L Getdphi.C+
#.L Getvn.C+
#Getdphi(icent, ihar, iangle1, iangle2, usingCNTEP)
#Getvn(icent, ihar, iangle1, iangle2, usingCNTEP)
