#!/bin/bash
out="$1"
#
place="$(dirname "$0")"
#
[ -z "$TMPDIR" ] && TMPDIR="/tmp"
trap "rm -f $TMPDIR/$$.run $TMPDIR/$$.ref" 0
#
awk -f ${place}/summarize.awk     "$out" > "${TMPDIR}/$$.run"
awk -f ${place}/summarize.awk "ref/$out" > "${TMPDIR}/$$.ref"
#
diff -c "${TMPDIR}/$$.run" "${TMPDIR}/$$.ref"
pass=$?
if [ $pass -eq 0 ] ; then
  echo "OK"
else
  echo "FAIL"
fi
exit $pass
