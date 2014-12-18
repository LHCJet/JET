JETJet
=================================================================
An implementation of JET algorithm described in arXiv:1411.3705.

Compiling (Requiring CMake and compiler with C++11 support):
--------------------------------------------
```sh
mkdir build
cd build
# Replace "/usr" with any path you prefer to install the libraries and headers
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr
make
# Optional. Install the libraries and headers
make install
```
A binary file "jetjet_example" will be generated in the build directory.

Usage:
------
To test the algorithm with default setting, simply enter the build directory and run:
```sh
./jetjet_example < ../sample/ttbar.txt
```
ttbar.txt is a sample file contains the four momentum of the particles in one ttbar event. There are a few more sample files in the sample directory for various scenarios. The output of jetjet_example uses the same format as print_jets_for_root function in FastJet.

Algorithms:
----------
Three different algorithms have been implemented so far in the code, by default jetjet_example will use the jet function Et - beta*m^2/Et, with beta=6. To use a different beta, find the line below in src/example.cpp
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
