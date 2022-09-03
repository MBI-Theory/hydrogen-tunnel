#!/bin/awk -f
/ Real kind = / {
  digits = substr($5,2) + 1 ;
  if (digits<=16) {
    en_fmt="%.12f" ;
    amp_fmt="%.10f" ;
    }
  else {
    en_fmt="%.24f" ;
    amp_fmt="%.20f" ;
    }
  }
function magic(fmt,val, str) {
  str = sprintf(fmt,val+0.0) ;
  if ( str ~ /^[0.+-]*$/ ) str = "0.0" ;
  return str ;
  }
/Outgoing solution is at the energy = / {
  printf "Energy = %s %s\n", magic(en_fmt,$8), magic(en_fmt,$9) ;
  }
/Outgoing amplitude in the main channel = / {
  printf "Out = %s %s\n", magic(amp_fmt,$8), magic(amp_fmt,$9) ;
  }
/Incoming amplitude in the main channel = / {
  printf "In = %s %s\n", magic(amp_fmt,$8), magic(amp_fmt,$9) ;
  }
