# Counting spinc structures on real Bott manifolds

The applications counts the numbers of real Bott manifolds with spinc and spin structures in dimensions from **3 to 10**.

## Requirements

- OpenMP
- gcc compiler
- GNU make

## Compiling

Simply type `make`.

## Invoking

After compiling, for both `count` and `backtrack` commands you may use the following switches:

`-d`: with argument being the dimension you want to check  
`-j`: with argument being number of threads to use  
`-v`: increase communicates that are printed after job is done  
`-h`: show short help  

For `backtrack` you have also:

`-s`: dimension you start the calculations with (default one less than the target dimension, but less than 12)
`-a`: wether you want to calculate number of manifolds with spin structures also
