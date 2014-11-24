JETJet
=================================================================
An implementation of JET algorithm described in arXiv:1411.3705.

Compiling (Requiring CMake and compiler with C++11 support):
--------------------------------------------
```sh
$ mkdir build
$ cd build
$ cmake ../ -DCMAKE_BUILD_TYPE=Release
$ make
```
A binary file "jetjet" will be generated in the build directory.

Usage:
------
To test the algorithm with default setting, simply enter the build directory and run:
```
$ ./jetjet < [input_file]
```
input_file is the file contains the four momentum of the particles in one event. There are a few sample files in the sample directory for various scenarios. The output of jetjet uses the same format as print_jets_for_root function in FastJet.

Algorithms:
----------
Three different algorithms have been implemented so far in the code, by default jetjet will use the jet function Et - beta*m^2/Et, with beta=6. To use a different beta, find the line below in src/example.cpp
```C++
JETJet::JetDefinition * jetDefinition = new JETJet::EtConeDefinition(6.0);
```
and change the number 6.0 to some other value.

You can replace the line above with
```C++
JETJet::JetDefinition * jetDefinition = new JETJet::EConeDefinition(beta);
```
This will invoke the jet algorithm defined in arXiv:1408.1161.

The third option is to use
```C++
JETJet::JetDefinition * jetDefinition = new JETJet::EtAlphaConeDefinition(alpha, beta);
```
This is the generalized algorithm described in Appendix A of arXiv:1411.3705.
