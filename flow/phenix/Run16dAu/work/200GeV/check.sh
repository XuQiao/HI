#!/bin/sh
for i in `cat tree1.lst`; do
    if [[ ! -e $i ]]; then
        echo -n  $i " "
    fi
done
