To compile: 
```
g++ -fopenmp -I../../../src/ -I../../../utils/ -I/usr/include/eigen3/ test.cc -o test -fopenmp
```

To run: 
```
./test
```
* `RungeKutta4.hh` has to be replaced in `utils` folder for only single integrator dynamics for 2x faster computation.
