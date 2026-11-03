#!/bin/bash
#
# short ConMan run as a build test
#
# generates a small thermal convection input deck with the CIG awk
# script, runs it in a temporary directory, and checks that the
# temperature field output was written
#
# usage: install/test.conman.sh [nel] [nsteps] [conman_dir]
#
nel=${1-16}			# elements in z
nsteps=${2-20}			# time steps
conman_dir=${3-modules/mc/ConMan}	# ConMan checkout relative to this directory

here=$(pwd)
conman_dir=$(cd "$conman_dir" && pwd)
exe=$conman_dir/conman
awkscript=$conman_dir/input_gen/make_conman_thermal_in.awk

if [ ! -x "$exe" ]; then
    echo "$0: $exe not found, run make conman first"
    exit 1
fi
if ! command -v gawk > /dev/null 2>&1; then
    echo "$0: gawk is required"
    exit 1
fi

tmpdir=$(mktemp -d /tmp/seatree.conman.test.XXXXXX)
cd "$tmpdir" || exit 1

name=test
opt="-v nel=$nel -v aspect=1 -v rayleigh=1e5 -v heating=0 -v viscosityE=0"
opt="$opt -v ntimestep=$nsteps -v nprintstep=$nsteps"

gawk $opt -v print_geom=1 -f "$awkscript" > geom.$name
gawk $opt -v print_geom=0 -f "$awkscript" > in.$name

cat <<EOF > run.$name
in.$name
geom.$name
out.$name
r.in
r.$name
tser.$name
temp.$name
stress.$name
comp.$name
geoid.$name
EOF

"$exe" < run.$name > conman.log 2>&1
status=$?

if [ $status -ne 0 ]; then
    echo "$0: conman exited with status $status, see $tmpdir/conman.log"
    exit 1
fi
if [ ! -s temp.$name ]; then
    echo "$0: no temperature output written, see $tmpdir/conman.log"
    exit 1
fi

nlines=$(wc -l < temp.$name)
echo "$0: ok, conman ran $nsteps steps on a $nel x $nel grid, $nlines lines of field output"
cd "$here" || exit 1
rm -rf "$tmpdir"
exit 0
