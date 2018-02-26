#! /bin/bash

make -j;
for n in 10 20 40 80 160;
do
./main $n 0 1;
done

