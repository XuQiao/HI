#!/bin/sh
read -r -p "Are you sure to submit all with ${1}? [y/N] " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
then
        cd 200GeV
        ./submit.sh ${1}
        cd ../62GeV
        ./submit.sh ${1}
        cd ../39GeV
        ./submit.sh ${1}
        cd ../20GeV
        ./submit.sh ${1}
  else
        echo "exit"
fi
