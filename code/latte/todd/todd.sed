## Sed script to change Maple expressions to C++ expressions using GMP rationals.
# Change powers 
s/\(x\[[0-9]*\]\)^\([0-9]*\)/pow(\1,\2)/g
# Rational literals
s|\([0-9][0-9]*\)/\([0-9][0-9]*\)|mpq_class(\1,\2)|g;
