JETJet
=================================================================
An implementation of JET algorithm described in arXiv:1411.3705.

Requirements:
--------------------------------------------
```
CMake > 2.8
GCC > 4.8 (or other compilers with C++11 support)
FastJet > 3.0 (optional)
Cython > 0.2 (optional)
```

Quick Start:
--------------------------------------------
After Unpack the source codes, enter the newly created directory, and execute the following command:
```sh
mkdir build
cd build
# Replace "/usr" with any path you prefer to install the libraries and headers
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr
make
# Optional. Install the libraries and headers
make install
#The output of jetjet_example uses the same format as print_jets_for_root function in FastJet
examples/jetjet_example < ../sample/dijet.txt
# If FastJet plugin is compiled
examples/fastjet_plugin_example < ../sample/ttbar.txt
# If Python Module is installed, the script simply print the four momentum of the clustered jets
python ../examples/python_example.py ../sample/large_eta.txt
```

Algorithms:
----------
Three different algorithms have been implemented so far in the code, by default jetjet_example will use the jet function Et - beta*m^2/Et, with beta=6. To use a different beta, find the line below in examples/jetjet_example.cpp
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
