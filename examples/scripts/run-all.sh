#!/bin/bash
type="${1:-dble}"
place="$(dirname ${0})"
#
hyd="../hydrogen_tunnel_v2.x"
gen="../general_tunnel.x"
chk="${place}/check.sh"
#
export OMP_THREAD_LIMIT=$(nproc)
#
function run_test () {
  local exe=$1 ; shift ;
  local inp=$1 ; shift ;
  local out ;

  out="$(basename -s ".inp" "$inp").out"
  if [ -r ${out} ] ; then
    echo "${out} already exists, reusing"
  else
    echo "Running ${exe} < ${inp} > ${out}"
    ${exe} < ${inp} > ${out}
  fi
  echo "Verifying ${out}"
  ${chk} ${out}
  if [ $? -eq 0 ] ; then
    pass="$pass $inp"
  else
    fail="$fail $inp"
  fi
  }
#
pass=""
fail=""
for inp in [^g]*-${type}.inp ; do
  run_test ${hyd} ${inp}
done
for inp in g*-${type}.inp ; do
  run_test ${gen} ${inp}
done
#
[ ! -z "$pass" ] && echo "Passed: $pass"
[ ! -z "$fail" ] && echo "Failed: $fail"
#
if [ -z "$fail" ] ; then
  exit 0
else
  exit 1
fi
