#!/bin/sh
# Produce a map of the misfit for a range of values of phi and dt compared
# to a constant set

phi_fixed=0
dt_fixed=1

phi_low=-90
phi_high=90
dphi=5

dt_low=0
dt_high=5
ddt=0.2

spol=5  # Keep it to whole numbers for shell arithmetic

R=$phi_low/$phi_high/$dt_low/$dt_high
I=$dphi/$ddt

GRD=$(mktemp /tmp/test_sw_misfit.grdXXXXXX || { echo "Can't make grid file"; exit 1; })
CPT=$(mktemp /tmp/test_sw_misfit.cptXXXXXX || { echo "Can't make cpt file"; exit 1; })
FIG=$(mktemp /tmp/test_sw_misfit.psXXXXXX || { echo "Can't make PS file"; exit 1; })

makecpt -I -Cwysiwyg -T0/1/0.01 -Z > $CPT && 
awk 'BEGIN{
	for (p='$phi_low'; p<='$phi_high'; p+='$dphi') {
		for(d='$dt_low'; d<='$dt_high'; d+='$ddt') print '$phi_fixed','$dt_fixed',p,d
	}}' |
	../progs/bin/split_misfit -spol $spol - 2>/dev/null |
	awk '{print $3,$4,$6}' 2>/dev/null |
	xyz2grd -R$R -I$I -G$GRD &&
{
	grdimage $GRD -JX10c/10c -R-90/90/0/5 -C$CPT \
		-Ba30f10:"@~f":/a1f0.2:"@~d@~t":nSeW -P -K
	echo "0 1" | psxy -N -J -R -O -K -Sx0.3c -W2p
	cat <<-END | psxy -J -R -O -W1p -m
	>
	$spol $dt_low
	$spol $dt_high
	>
	$((spol+90)) $dt_low
	$((spol+90)) $dt_high
	>
	$((spol-90)) $dt_low
	$((spol-90)) $dt_high
	END
} > $FIG && gv $FIG 2>/dev/null

rm -f $GRD $CPT $FIG
