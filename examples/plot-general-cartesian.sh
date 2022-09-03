#!/bin/bash
#
#  Plot a 2D slice for the Cartesian interpolant of the general_tunnel.x solution
#
tab="${1:-"missing.table"}"
gpl="gpl/$(basename -s ".table" "$tab").gpl"
pdf="pdf/$(basename -s ".table" "$tab").pdf"
#
cat > "${gpl}" <<-eoi
	set terminal x11 0 enhanced
	# set terminal pdf enhanced color
	# set output "$pdf"
	set view equal xy
	set grid
	set surface
	set contour both
	set title "$tab"
	set xlabel "X, a_0"
	set ylabel "Z, a_0"
	set zlabel "Wavefunction"
	splot \
	  "< awk '!/^#/&&(\$2==0)&&(\$3!=oz){print \"\";}!/^#/&&(\$2==0){print;oz=\$3}' '${tab}'" u 4:6:8 ti "Re" with lines ls 1, \
	  1/0 noti
	pause mouse keypress
	splot \
	  "< awk '!/^#/&&(\$2==0)&&(\$3!=oz){print \"\";}!/^#/&&(\$2==0){print;oz=\$3}' '${tab}'" u 4:6:9 ti "Im" with lines ls 1, \
	  1/0 noti
	pause mouse keypress
	set terminal x11 1 enhanced
	set logscale z
	splot \
	  "< awk '!/^#/&&(\$2==0)&&(\$3!=oz){print \"\";}!/^#/&&(\$2==0){print;oz=\$3}' '${tab}'" u 4:6:(column(8)**2+column(9)**2) ti "ρ" with lines ls 1, \
	  1/0 noti
	pause mouse keypress
eoi
gnuplot -persist "$gpl"
