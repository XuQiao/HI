#!/bin/csh
root -b << EOF
.L Getvn.C+
#include "RpPar.h"
int code = $1
int iangle1 = 0
int iangle2 = 0
bool UseCNTEP = 0
Getvn(code, iangle1, iangle2, UseCNTEP)
EOF
echo "Done!"
