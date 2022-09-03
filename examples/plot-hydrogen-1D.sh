#!/bin/bash
#
#  Plot 1D wavefunction
#
tab="${1:-"missing.table"}"
gpl="gpl/$(basename -s ".table" "$tab").gpl"
pdf="pdf/$(basename -s ".table" "$tab").pdf"
#
rstop="$(awk '/^#  *R.ipt_stop. = /{print $4}' "$tab")"
fpbas="$(awk '/^#  *FP radix = /{print $5}' "$tab")"
ymin="$(awk 'function min(x,y){return x<y?x:y}!/^#/&&($2<rstop){x=min(x,min($4,$5)*fpbas**$3)}END{print x}' rstop="$rstop" fpbas="$fpbas" "$tab")"
ymax="$(awk 'function max(x,y){return x>y?x:y}!/^#/&&($2<rstop){x=max(x,max($4,$5)*fpbas**$3)}END{print x}' rstop="$rstop" fpbas="$fpbas" "$tab")"
#
cat > "${gpl}" <<-eoi
	set terminal x11 0 enhanced
	# set terminal pdf enhanced color
	# set output "$pdf"
	set grid
	set title "$tab"
	set xlabel "Coordinate, a_0"
	set ylabel "Wavefunction"
	set arrow 100 from ${rstop}, graph 0.0 to ${rstop}, graph 1.0 nohead back lw 2.0
	ymin=$ymin
	ymax=$ymax
	set yrange [ymin-0.02*(ymax-ymin):ymax+0.02*(ymax-ymin)]
	plot \
	  "${tab}" u 2:(column(4)*${fpbas}**column(3)) ti "Re" with lines, \
	  "${tab}" u 2:(column(5)*${fpbas}**column(3)) ti "Im" with lines, \
	  1/0 noti
	set terminal x11 1 enhanced
	set yrange [*:*]
	set logscale y
	plot \
	  "${tab}" u 2:((column(4)**2+column(5)**2)*${fpbas}**(2.0*column(3))) ti "{/Symbol r}" with lines, \
	  1/0 noti
eoi
gnuplot -persist "$gpl"
