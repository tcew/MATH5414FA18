
### MATH 5414 Fall 2018 (Virginia Tech)

This repository is under active development as part of the Finite Elements and GPU computing graduate topics course offered by the Math Department at Virginia Tech by [Tim Warburton](http://www.paranumal.com). Each function is live coded during class to illustrate aspects of developing components of a multi-GPU accelerated finite element code using MPI and OCCA.

A complete multi-GPU accelerated finite element solver called libParanumal is available here [here](https://github.com/paranumal/libparanumal).

 ```To build:  
  
make  ```
 
```To run on 4 processes:  
  
mpiexec -n 4 ./meshMain meshes/Lshape2H01.msh ```



