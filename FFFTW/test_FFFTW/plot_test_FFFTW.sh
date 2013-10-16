#!/bin/bash
# Run test and show Gnuplot graph to display results

./test_FFFTW &&
{
# Plot of input and output plus the residual*1e5
cat << END | gnuplot --persist
set xlabel "Time / s"
set ylabel " Amplitude"
set title "Time-domain trace"
p   "input.xy" w l lt 2 title "Input trace", \
    "output.xy" w l lt 0 title "Output trace", \
    "< paste input.xy output.xy | awk '{print \$1,1e5*(\$4-\$2)}'" \
		w l lt 3 title "10**5 * resid"
END

# Plot of the amplitdue spectrum
cat << END | gnuplot --persist
set xlabel "Frequency / Hz"
set ylabel "Amplitude"
set title "Amplitude spectrum"
p "amplitude.xy" w l lt 1 title "Amplitude"
END
}