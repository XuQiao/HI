#!/bin/sh
filename=(
# "/store/caf/user/tuos/pPb/Glauber/Default/Phob_Glau_pPb_sNN70mb_v15_1M_dmin04.root"
# "/store/caf/user/tuos/pPb/Glauber/SigmaNN/sigma65/Phob_Glau_pPb_sNN65mb_v15_1M_dmin04.root"
# "/store/caf/user/tuos/pPb/Glauber/SigmaNN/sigma75/Phob_Glau_pPb_sNN75mb_v15_1M_dmin04.root"
# "/store/caf/user/tuos/pPb/Glauber/Radius/R675/Phob_Glau_pPb_sNN70mb_v15_1M_R675.root"
# "/store/caf/user/tuos/pPb/Glauber/Radius/R649/Phob_Glau_pPb_sNN70mb_v15_1M_R649.root"
# "/store/caf/user/tuos/pPb/Glauber/dmin/dmin00/Phob_Glau_pPb_sNN70mb_v15_1M_dmin00.root"
# "/store/caf/user/tuos/pPb/Glauber/dmin/dmin08/Phob_Glau_pPb_sNN70mb_v15_1M_dmin08.root"
# "/store/caf/user/tuos/pPb/Glauber/skinDepth/D04914/Phob_Glau_pPb_sNN70mb_v15_1M_D04914.root"
# "/store/caf/user/tuos/pPb/Glauber/skinDepth/D06006/Phob_Glau_pPb_sNN70mb_v15_1M_D06006.root"

"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega055_D04914.root"
"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega055_D06006.root"
"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega055_R649.root"
"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega055_R675.root"
"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega055_default.root"
"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega055_dmin00.root"
"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega055_dmin08.root"
"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega055_sigma75.root"
"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega055_sigmaNN65.root"

"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega101_D04914.root"
"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega101_D06006.root"
"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega101_R649.root"
"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega101_R675.root"
"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega101_default.root"
"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega101_dmin00.root"
"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega101_dmin08.root"
"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega101_sigma75.root"
"/store/caf/user/tuos/pPb/glauberGribov/glaubergribovpPbv1omega101_sigmaNN65.root"
)

	cd ~/CMSSW_5_3_8_HI/src
	eval `scramv1 runtime -sh`
	cd ~/CMSSW_5_3_8_HI/src/UserCode/Centrality/Correction
for file in ${filename[*]}
do
echo $file
xrdcp root://xrootd.unl.edu/$file /lio/lfs/cms/store/user/qixu/Glauber
done
