#!/bin/bash
# Rotate a tensor using each of the six Euler rotation conventions.

which CIJ_rot_euler 2>/dev/null 1>&2 ||
	{ echo "Cannot find CIJ_rot_euler in PATH"; exit 1; }

modes="xzx xyx yxy yzy zyz zxz"

a=10
b=$a
c=$a

ECs="5.48771e+07  1.99771e+07  2.06714e+07  6.00000e+04 -1.14286e+04  2.85714e+04  1.99771e+07  6.77371e+07  2.09914e+07 -8.85714e+04  7.14286e+04  1.74286e+05  2.06714e+07  2.09914e+07  5.96429e+07 -8.57143e+04  6.57143e+04 -8.00000e+04  6.00000e+04 -8.85714e+04 -8.57143e+04  2.07286e+07  2.85714e+03  1.08571e+05 -1.14286e+04  7.14286e+04  6.57143e+04  2.85714e+03  1.80800e+07  2.57143e+04  2.85714e+04  1.74286e+05 -8.00000e+04  1.08571e+05  2.57143e+04  1.95714e+07"

for mode in $modes; do
	echo $ECs | CIJ_rot_euler -t $mode $a $b $c | CIJ_plot -t $mode -r 1
done
