#!/bin/bash
#
#  Plot single-channel bound solutions for general_tunnel.x
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
	set surface
	set contour both
	set title "$tab"
	set xlabel "Coordinate η, a_0"
	set xlabel "Coordinate ξ, a_0"
	set zlabel "Wavefunction"
	splot \
	  "${tab}" u 2:4:5 ti "Re" with lines, \
	  "${tab}" u 2:4:6 ti "Im" with lines, \
	  "< awk '/#  *ETA = /{e=\$4}/#  *ipt_stop =/{ir=\$4}!/^#/&&(\$3==ir){print e, \$4, \$5**2+\$6**2}' '${tab}'" u 1:2:(sqrt(\$3)) ti "stop" with lines lw 3.0 lc "red", \
	  1/0 noti
	pause mouse keypress
	set terminal x11 1 enhanced
	set logscale z
	splot \
	  "${tab}" u 2:4:(column(5)**2+column(6)**2) ti "{/Symbol r}" with lines, \
	  "< awk '/#  *ETA = /{e=\$4}/#  *ipt_stop =/{ir=\$4}!/^#/&&(\$3==ir){print e, \$4, \$5**2+\$6**2}' '${tab}'" u 1:2:3 ti "stop" with lines lw 3.0 lc "red", \
	  1/0 noti
	pause mouse keypress
eoi
gnuplot -persist "$gpl"
