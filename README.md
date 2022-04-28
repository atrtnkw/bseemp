BSEEMP is a binary population synthesis code to follow binary
evolution. Please find the manual here:
https://github.com/atrtnkw/bseemp/tree/main/doc/ms.pdf

The following is a tentative manual of the initial condition generator
in the "init" directory. The manual will be combined with the BSEEMP
manual in https://github.com/atrtnkw/bseemp/tree/main/doc/ms.pdf soon.

After entering into the "init" directory, please do "make". If
successful, you will find the executable file "run". After you do
"./run input.dat", you will find files "ibinary.dat" and
"ibinary.txt". The "ibinary.dat" is a file containing binary initial
conditions, and the "ibinary.txt" is a file describing a breif
explation of the file "ibinary.dat". By default, the binary initial
conditions include 10^6 binary stars with Z=2e-04. If you use them as
initial conditions of popbin2.bonn/popbin2.geneva, you will calculate
100Gyr evolution of the 10^6 binary stars. The binary initial
conditions follow the Kroupa (2001)'s IMF for the primary stellar
masses, and Sana et al. (2012)'s binary models for the distributions
of mass ratios, periods, and eccentricities. If you want to change the
binary initial conditions, please edit the file "input.dat". You can
change the number of binary stars, metallicities, IMFs, binary
parameters, and so on. The detail explanation will be described in
https://github.com/atrtnkw/bseemp/tree/main/doc/ms.pdf later.
