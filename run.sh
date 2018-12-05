#!/bin/sh   
date
cd T0859-D1/outputFolder/
for i in T0859-D1\
        T0862-D1\
        T0863-D1\
        T0863-D2\
        T0864-D1\
        T0866-D1\
        T0868-D1\
        T0869-D1\
        T0870-D1\
        T0886-D1
do
{
    echo $i
    cd ../../$i/outputFolder/
    ../../../rosetta_src_2016.32.58837_bundle/main/source/bin/AbinitioRelax.linuxgccrelease @../inputFolder/flags
    ../../../spicker/spicker
    ../../../pulchra/pulchra.pl combo1.pdb
    mv pul_center.pdb combo1_pul.pdb
    ../../../pulchra/pulchra.pl combo2.pdb
    mv pul_center.pdb combo2_pul.pdb
    ../../../pulchra/pulchra.pl combo3.pdb
    mv pul_center.pdb combo3_pul.pdb
    ../../../pulchra/pulchra.pl combo4.pdb
    mv pul_center.pdb combo4_pul.pdb
    ../../../pulchra/pulchra.pl combo5.pdb
    mv pul_center.pdb combo5_pul.pdb
} &
done
wait
date
