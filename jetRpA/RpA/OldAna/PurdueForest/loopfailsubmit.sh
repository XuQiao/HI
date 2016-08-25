cd /home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/PurdueForest/
grid-proxy-init -valid 24:00
grid_result=$?
echo $grid_result
if [ $grid_result -eq 0 ]
        then
        if [ -f /tmp/x509up_u171345 ]
                then
                cp /tmp/x509up_u171345 /home/xuq7/tmp/
        else
                echo "No grid file exists!"
                exit
        fi

trigname="MB"
cd $trigname/script
for pbsname in `ls`
do
if [[ -f $pbsname && $pbsname =~ "fail" ]]; then
	qsub $pbsname
fi
done
fi
