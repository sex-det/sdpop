SDpop
=====

Author: Jos Kafer, jos.kafer@cnrs.fr

SDpop infers sex-linkage from genotyping data of several individuals of both sexes, collected in panmictic populations. 

To install, make sure you have a sufficiently recent C++ compiler (understanding C++11) and type:
make all
If you're on a computer cluster, you might want to compile statically, which you could try to do with:
make CXXFLAGS="--static" all

The software is provided "as is" and is still under development. 
It is developed by me (Jos Käfer); I'm a research scientist interested in plant evolutionary biology, not a software 
developer interested in optimizing code. Thus, the programs provided here are aimed to work on my data and with my 
pipelines. I've tried to describe the data formats and caveats, but certainly not in an exhaustive manner. Error 
handling in the code is rudimentary, so running the programs on you data will very likely cause segmentation faults (or 
worse).

You are free to use the code, as long as you cite the paper in your work. If you run into any issue, and you're not 
familiar enough with C/C++ to solve it yourself (or not enduring enough to scroll through my thousands of lines of 
ill-commented code, which I understand), feel free to contact me. I'm happy to make the code more robust, and learn 
something about your research. 

If you're a software developer interested in optimizing the code, I welcome your help and suggestions. The code can be 
written more concisely, and several algorithmic techniques could greatly diminish the running time and memory usage. One
of the first things to do would be to compress the data and to do calculations for sites with the same observed 
genotypes only once.
