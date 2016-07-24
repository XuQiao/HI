#!/bin/sh
	for j in {0..10}
	do
		echo -n `sed -n '5809p' out_job0$j`
		echo -n ',' 
	done
		echo " "
	for j in {0..10}
	do
                echo -n `sed -n '923p' out_job1$j`
                echo -n ',' 
	done
		echo " "
	for j in {0..10}
	do
                echo -n `sed -n '923p' out_job2$j`
                echo -n ',' 
	done
		echo " "
