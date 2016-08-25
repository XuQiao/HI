#!/bin/sh
for i in {0..10};do
 	echo -n `grep "{" "out0$i"`,
done
echo
for i in {0..10};do
 	echo -n `grep "{" "out1$i"`,
done
echo
for i in {0..10};do
 	echo -n `grep "{" "out2$i"`,
done
echo
