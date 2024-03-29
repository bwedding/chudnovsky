# Chudnovsky Algorithm for Calculating Pi

This program calculates the value of Pi using the [Chudnovsky algorithm](https://en.wikipedia.org/wiki/Chudnovsky_algorithm),
an extremely fast method to calculate the digits of π. It was published by the Chudnovsky brothers in 1987.

## Prerequisites

The program requires two libraries:

1. **Boost:** A set of peer-reviewed portable C++ source libraries that extend the functionality of C++. You can download it from [here](https://www.boost.org/users/download/). 
2. **MPFR:** A library for multiple-precision floating-point computations with correct rounding. You can download it from [here](http://www.mpfr.org/mpfr-current/#download).

Make sure to install these libraries before compiling and running the program. Both packages are also available on [brew.sh](https://brew.sh).

## How To Compile The Program?
You need a compiler that supports C++20 standard because this code uses features introduced in this version.
Here is an example compile line: `c++ -std=c++20 -O3 -o chudnovsky chudnovsky.cpp -L/usr/local/opt/mpfr/lib/ -lmpfr`
This was used on my i9 MacBook Pro.

On an M2 Mac mini, I used the following line: `c++ -std=c++20 -O3 -I /opt/homebrew/include -L /opt/homebrew/lib -o chudnovsky chudnovsky.cpp -lmpfr`

This is because HomeBrew installs packages on Intel based Macs differently than M* based Macs. Or perhaps it is because of
Big Sur. Whatever. Who cares? It's different and it's annoying!

## Usage
To run the program, you should pass one argument which is the number of terms to be computed in the series.
Usage example: `./chudnovsky 10` 

It will compute Pi up to 10 terms.

The output will display Pi computed up to `n` terms (where `n` was your input) and then print out a calculated approximation of Pi.

## Algorithm Explanation
The Chudnovsky algorithm calculates 1/π as:

![12 Σ ((-1)^k * (6k)! * (545140134k + 13591409)) / ((3k)!(k!)^3 * (640320)^(3k+3/2))](https://github.com/Fudmottin/Chudnovsky/blob/main/chudnovsky.svg)

Where the summation Σ is from k = 0 to ∞.

The program uses Boost's multiprecision and special function libraries to calculate factorials and powers. It also uses MPFR for multiple-precision floating-point computations.

## Note

Remember that the more terms you use, the more accurate your calculation will be, but it will also take longer to compute. For instance, using 10 terms provides a very accurate approximation of Pi.

This program is a great way to understand how numerical methods work and how C++ can be used for scientific computing tasks. Feel free to experiment with it and optimize it as needed!

## Addendum

```
Pi-1000000.txt
pi-diff.sh
```

These files are for testing the output for correctness. The Pi-1000000.txt file is
the first million digits of π in a single line of text. The shell script is intended
to compare the output of chudnovsky to the canonical π text to see where disagreement
begins.

```
./chudnovsky 50 > pi.txt  0.00s user 0.00s system 75% cpu 0.007 total
$ sh pi-diff.sh                                                                    
cmp: EOF on pi.txt
Difference found at position: 711
Displaying ten characters from each file starting five chars before the difference...
           9   9   5   6   1   1   2   1   2   9                        

           9   9   5   6   2   4   3   8   6   5                        

$ wc -c pi.txt 
     752 pi.txt
$ time ./chudnovsky 100 > pi.txt
./chudnovsky 100 > pi.txt  0.01s user 0.00s system 84% cpu 0.014 total
$ sh pi-diff.sh                 
cmp: EOF on pi.txt
Difference found at position: 1421
Displaying ten characters from each file starting five chars before the difference...
           0   9   3   4   1   7   2   1   6   4                        

           0   9   3   4   9   4   4   5   8   6                        

$ wc -c pi.txt                  
    1502 pi.txt
$ time ./chudnovsky 200 > pi.txt
./chudnovsky 200 > pi.txt  0.04s user 0.00s system 95% cpu 0.045 total
$ sh pi-diff.sh                 
cmp: EOF on pi.txt
Difference found at position: 2838
Displaying ten characters from each file starting five chars before the difference...
           8   4   7   8   4   8   9   6   8   3                        

           8   4   7   8   5   2   7   1   0   4                        

$ wc -c pi.txt                  
    3002 pi.txt
```
As you can see, there is a point where the digits are start to be just plain wrong.
that's just the nature of this implementation. Still, it was fun to write. Debugging
was a bit of a challenge. Integer overflows do terrible things! If you follow the
link to the wikipedia page, you will see that there is more room for optimization.
Not only can the number of operations be reduced, the algorithm can be parallelized
to take advantage of multiple cores. This implementation is parallelized using the
C++17 standard threading mechanism. Even so, further optimization is possible.
