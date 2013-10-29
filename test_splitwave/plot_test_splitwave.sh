#!/bin/bash
# Using SAC and SHEBA, plot the results of the test_splitwave program.

macrodir=~/Applications/sacmacros

./test_splitwave && {

	cat <<-END | sac
	setmacro $macrodir
	window 1 x 0 0.5 y 0.3 1
	bd x
	m sheba file wave_fdsplit nwind one pick no
	m sheba file wave_tdsplit nwind one pick no
	q
	END
	
	cleansheba

} || { echo "Problem running test_splitwave"; exit 1; }
