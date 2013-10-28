#!/bin/bash
# Run test and show Gnuplot graph to display results

./test_FFFTW &&
{
for a in real complex; do

	# Plot of input and output plus the residual*1e5
	cat << END | gnuplot --persist
	set xlabel "Time / s"
	set ylabel "Amplitude"
	set xrange [0:5]
	set title "Time-domain trace (${a})"
	p   "input_${a}.xy" w l lt 2 title "Input trace", \
	    "output_${a}.xy" w l lt 0 title "Output trace", \
	    "< paste input_${a}.xy output_${a}.xy | awk '{print \$1,1e5*(\$4-\$2)}'" \
			w l lt 3 title "10**5 * resid"
END

	# Plot of the amplitdue spectrum
	cat << END | gnuplot --persist
	set xlabel "Frequency / Hz"
	set ylabel "Amplitude"
	set title "Amplitude spectrum (${a} input)"
	p "amplitude_${a}.xy" w l lt 1 title "Amplitude"
END
done
}