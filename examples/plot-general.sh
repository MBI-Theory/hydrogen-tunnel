#!/bin/bash
#
#  Plot 2D section of the total solution of general_tunnel.x, parabolic coordinates
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
	  "${tab}" u (+column(6)):7:8 ti "Re" with lines ls 1, \
	  "${tab}" u (-column(6)):7:8 noti    with lines ls 1, \
	  1/0 noti
	pause mouse keypress
	splot \
	  "${tab}" u (+column(6)):7:9 ti "Im" with lines ls 2, \
	  "${tab}" u (-column(6)):7:9 noti    with lines ls 2, \
	  1/0 noti
	pause mouse keypress
	set terminal x11 1 enhanced
	set logscale z
	splot \
	  "${tab}" u (+column(6)):7:(column(8)**2+column(9)**2) ti "{/Symbol r}" with lines ls 3, \
	  "${tab}" u (-column(6)):7:(column(8)**2+column(9)**2) noti             with lines ls 3, \
	  1/0 noti
	pause mouse keypress
eoi
gnuplot -persist "$gpl"
