#! /bin/bash

make -j;
#for n in 10 20 40 80;
for n in 400;
#for n in 4;
do
./main $n 0 0.02;
#./main $n 0 0.055;
done

