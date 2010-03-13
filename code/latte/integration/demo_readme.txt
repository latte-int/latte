Usage: integrate_demo [-m] [inFile] [outFile]
Arguments:
The -m option specifies that the input will consist of monomial sums and simplices. Otherwise, an input of sums of powers of linear forms and simplices will be expected.
If inFile is specified, the program will attempt to read the input from the file. Otherwise, the input is assumed to come from the console.
Also, if inFile is specified, an outFile may be specified into which the resulting output will be written. Otherwise, all output is displayed to the console.
Input formatting:
The input is expected to alternate between a string representation of a polynomial (as a sum of either monomials or powers of linear forms) and a string representation of a simplex.
For an arbitrary sum of monomials a_{0}x_{0}^e_{0, 0}x_{1}^e_{0, 1}...x_{n}^e_{0, n} + ... + a_{m}x_{m}^e_{m, 0}x_{1}^e_{m, 1}...x_{n}^e_{m, n},
    the corresponding string representation is: [ [ a_{0}, [e_{0, 0}, e_{0, 1}, ..., e_{0, n}] ], ..., [ a_{m}, [e_{m, 0}, e_{m, 1}, ..., e_{m, n}] ] ]
For an arbitrary sum of powers of linear forms a_{0}(c_{0, 0}x_{0} + c_{0, 1}x_{1}... + c_{0, n}x_{n})^e_{0] + ... + a_{m}(c_{m, 0}x_{0} + c_{m, 1}x_{1}... + c_{m, n}x_{n})^e_{m},
    the corresponding string representation is: [ [ a_{0}, [ e_{0}, [c_{0, 0}, c_{0, 1}, ..., c_{0, n}] ] ], ..., [ a_{m}, [ e_{m}, [c_{m, 0}, c_{m, 1}, ..., c_{m, n}] ] ]
For a simplex consisting of n-dimensional points (a_{0, 1}, ..., a_{0, n}), ..., (a_{n, 1}, ..., a_{n, n}),
    the corresponding string representation is: [ [ a_{0, 1}, ..., a_{0, n} ], ..., [ a_{n, 1}, ..., a_{n, n} ] ]