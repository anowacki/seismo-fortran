#!/bin/bash
#
# Test CIJ_flip? visually

ecs_orig=`mktemp /tmp/ecs_orig.XXXXXXXXXX`
ecs_flipped=`mktemp /tmp/ecs_flipped.XXXXXXXXXX`

# Show us the original constants
./test_flip x $ecs_orig $ecs_flipped

echo "Plotting constants before transformation"
CIJ_plot $ecs_orig 2>&1 | grep -v Density &
sleep 1

for m in x y z; do
	./test_flip $m $ecs_orig $ecs_flipped
	echo "Constants flipped in plane normal to $m"
	CIJ_plot $ecs_flipped 2>&1 | grep -v Density
done

wait
