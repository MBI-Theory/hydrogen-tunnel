#!/bin/bash
#
#  Plot single-channel continuum solution from general_tunnel.x
#
tab="${1:-"missing.table"}"
gpl="gpl/$(basename -s ".table" "$tab").gpl"
pdf="pdf/$(basename -s ".table" "$tab").pdf"
#
cat > "${gpl}" <<-eoi
	set terminal x11 0 enhanced
	# set terminal pdf enhanced color
	# set output "$pdf"
	set grid
	set title "$tab"
	set xlabel "Coordinate η, a_0"
	set ylabel "Wavefunction"
	plot \
	  "${tab}" u 2:(column(3)) ti "Re" with lines, \
	  "${tab}" u 2:(column(4)) ti "Im" with lines, \
	  1/0 noti
	set terminal x11 1 enhanced
	set yrange [*:*]
	set logscale y
	plot \
	  "${tab}" u 2:(column(3)**2+column(4)**2) ti "{/Symbol r}" with lines, \
	  1/0 noti
eoi
gnuplot -persist "$gpl"
