#!/bin/awk -f
/ Real kind = / {
  digits = substr($5,2) + 1 ;
  if (digits<=16) {
    en_fmt="%.11f" ;
    amp_fmt="%.8f" ;
    printf "Real\n" ;
    }
  else {
    en_fmt="%.23f" ;
    amp_fmt="%.18f" ;
    printf "Quad\n" ;
    }
  }
# reporting
function magic(fmt,val, str) {
  str = sprintf(fmt,val+0.0) ;
  if ( str ~ /^[0.+-]*$/ ) str = "0.0" ;
  return str ;
  }
function energy(re,im) {
  printf "Energy = %s %s\n", magic(en_fmt,re), magic(en_fmt,im) ;
  }
function amp(tag,re,im) {
  printf "%s = %s %s\n", tag, magic(amp_fmt,re), magic(amp_fmt,im) ;
  }
# common
/Outgoing solution is at the energy = / { energy($8,$9) ; }
# hydrogen_tunnel.x
/Outgoing amplitude = / { amp("Out",$4,$5) ; }
/Incoming amplitude = / { amp("In", $4,$5) ; }
# general_tunnel.x
/Located solution at energy = / { energy($6,$7) ; }
/Optimized solution [^ ]* is at the energy = / { energy($9,$10) ; }
/Outgoing amplitude in the main channel = / { amp("Out",$8,$9) ; }
/Incoming amplitude in the main channel = / { amp("In",$8,$9) ; }
