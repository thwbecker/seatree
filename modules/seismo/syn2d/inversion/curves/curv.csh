#!/bin/csh
#given f(x), calculate f'(x), f"(x), and the curvature of f(x)
#argument 1: xy file with x, f(x)
./bin/dumbder<<EOF
"$1"
2
2
1
EOF
./bin/dumbder<<EOF
"$1.der.02.wrt.01"
3
3
1
EOF
./bin/dumbcurv<<EOF
"$1.der.02.wrt.01.der.03.wrt.01"
EOF
