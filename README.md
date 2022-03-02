BSEEMP is a binary population synthesis code to follow binary
evolution. This code is based on the BSE code (Hurley et al. 2000;
Hurley et al. 2002), and supports single star evolution in a wider
parameter range than the original BSE: stellar metallicity down to
Z=2x10^-10, and stellar mass up to 10^5Msun. There are two kinds of
single star evolution models (the M and L models), whose detail
modeling can be seen in Yoshida et al. (2019, ApJ, 881, 16), and whose
features can be seen in Appendix of Tanikawa et al. (2022, ApJ, 926,
83,
https://ui.adsabs.harvard.edu/abs/2022ApJ...926...83T/abstract). The
difference between the M and L models is more pronouced in more
metal-poor stars (please see Tanikawa et al. 2021, MNRAS, 505, 217,
https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.2170T/abstract). The
detail implementation of this code is described in Tanikawa et
al. (2020, MNRAS, 495, 4170,
https://ui.adsabs.harvard.edu/abs/2020MNRAS.495.4170T/abstract).

This readme provides the usage of this code. The usage is described in
the following order: bse.geneva (a code to follow a binary star in the
M model), popbin2.geneva (a code to follow a bunch of binary stars in
the M model), bse.bonn (a code to follow a binary star in the L
model), and popbin2.bonn (a code to follow a bunch of binary stars in
the L model).

-- bse.geneva

You can compile "bse.geneva", and use it as follows.

1. Move to the "src" directory.

2. Excute the command "make clean" if you already compile bse.bonn or popbin2.bonn.

3. Compile "bse.geneva" by executing the command "make bse". Note that
gfortran and g++ compilers are required.

4. Move to the "example/bse.geneva" directory.

5. Execute the command "./bse.geneva".

6. You can get the following output:

      TIME      M1       M2   K1 K2        SEP    ECC  R1/ROL1 R2/ROL2  TYPE  
     0.0000   143.924   134.604  1  1   0.1188E+04     0.189024   0.012   0.012  INITIAL   0.6093E+00  0.6523E+00    1.000    1.000 NOSN NOSN     0.000     0.000  
     1.9122   143.852   134.539  4  1   0.1188E+04     0.189024   0.042   0.049  KW_CHNGE  0.6048E+00  0.6479E+00    1.000    1.000 NOSN NOSN     4.756     4.717  
     1.9581   143.848   134.534  4  4   0.1188E+04     0.189024   0.043   0.042  KW_CHNGE  0.6041E+00  0.6474E+00    1.000    1.000 NOSN NOSN     4.755     4.758  
     2.1677   137.665   134.654  5  4   0.1212E+04     0.188824   0.228   0.095  KW_CHNGE  0.2782E-01  0.6531E+00    1.000    1.000 NOSN NOSN     4.394     4.577  
     2.1748   143.852   134.695 15  4   0.0000E+00     0.000000   0.000  -1.000  NO_REMNT    Infinity  0.6528E+00    1.000    1.000 PISN NOSN     4.397     4.559  
     2.2177     0.000   129.815 15  5   0.0000E+00    -1.000000  -1.000   0.000  KW_CHNGE    Infinity  0.5020E-01    1.000    1.000 PISN NOSN     4.397     4.413  
     2.2250     0.000   134.534 15 15   0.0000E+00     0.000000   0.000  -2.000  NO_REMNT    Infinity    Infinity    1.000    1.000 PISN PISN     4.397     4.409  

time [Myr]  
star 1's mass [Msun]  
star 2's mass [Msun]  
star 1's type (described later)  
star 2's type  
binary semi-major axis  
binary eccentricity  
star 1's radius scaled by its Roche-lobe radius  
star 2's radius scaled by its Roche-lobe radius  
binary status (described later)  
star 1's dimensionless spin  
star 2's dimensionless spin  
star 1's spin-orbit inclination  
star 2's spin-orbit inclination  
star 1's supernova type (discribed later)  
star 2's supernova type  
star 1's effective temparature  
star 2's effective temparature  

7. You can set the initial condition of a binary star, editting the
file "binary.in". Its format is as follows.

mass0(1),mass0(2),tphysf,tb,kstar(1),kstar(2),z,ecc  
neta,bwind,hewind,alpha1,lambda,betaacc  
ceflag,tflag,ifflag,wdflag,bhflag,nsflag,psflag,mxns,idum  
NewStarModel,WindEnhanced,RadiusShrinkage,NewDynTide,NewMassTransfer 
pts1,pts2,pts3  
sigma,beta,xi,acc2,epsnov,eddfac,gamma  

mass0(1): star 1's mass [Msun]  
mass0(2): star 2's mass [Msun]  
tphysf: ending time [Myr]  
tb: binary period [day]  
kstar(1): star type  
kstar(2): star type  
z: metallicity (not normalized by the solar metallicity)  
ecc: binary eccentricity  
neta, bwind, hewind: stellar wind parameters  
alpha1, lambda: commone envelope parameters  
betaacc: efficiency of mass accretion  
bhflag: 0 (no BH natal kick) 1 (fallback BH natal kick)  
nsflag: 3 (Fryer's rapid model) 4 (Fryer's delayed model)  
psflag: 0 (no pair instability) 1 (standard pair instability) 2 (3-sigma pair instability)  
mxns: the maximum neutron star's mass  
NewStarModel: .true. (using the M model), .false. (using the original single star model)  
WindEnhanced: .false. (recommended)  
RadiusShrinkage: .false. (recommended)  
NewDynTide: .true. (Kinugawa et al. 2020's dynamical tide), .false. (the original dynamical tide)  
NewMassTransfer: .true. (Kinugawa et al. 2020's mass transfer), .false. (the original dynamical tide)  
pts1, pts2, pts3: parameters for timesteps  
sigma: NS kick velocity  
beta: efficiency of wind accretion  
xi:  
acc2:  
epsnov:  
eddfac: Eddington factor  
gamma:  

-- popbin2.geneva

You can compile "popbin2.geneva", and use it as follows.

1. Move to the "src" directory.

2. Excute the command "make clean" if you already compile bse.bonn or popbin2.bonn.

3. Compile "popbin2.geneva" by executing the command "make
popbin2". Note that gfortran and g++ compilers are required.

4. Move to the "example/popbin2.geneva" directory.

5. Execute the command "./popbin2.geneva".

6. You can get two files "binary0000001.txt" and "binary0000002.txt"
in the directory "output". The file "binary0000001.txt" includes 10000
binary star's evolution. The file format is the same as the file
format of "bse.geneva", except that "### 1" indicates binary ID (here,
binary 1).

7. You can set parameter sets and initial conditions of binary stars,
editting the files "header.in" and binaries.in", respectively.

The format of "header.in" is the same as "binary.in" except for the
first line.

z  
neta,bwind,hewind,alpha1,lambda,betaacc  
ceflag,tflag,ifflag,wdflag,bhflag,nsflag,psflag,mxns,idum  
NewStarModel,WindEnhanced,RadiusShrinkage,NewDynTide,NewMassTransfer 
pts1,pts2,pts3  
sigma,beta,xi,acc2,epsnov,eddfac,gamma  

The format of "binaries.in" is as follows.

nbinary  
mass0(1) mass0(2) tb ecc tphysf tphys  
mass0(1) mass0(2) tb ecc tphysf tphys  
mass0(1) mass0(2) tb ecc tphysf tphys  
...  

nbinary: the number of binary stars
mass0(1): star 1's mass [Msun]  
mass0(2): star 2's mass [Msun]  
tb: binary period [day]  
ecc: binary eccentricity  
tphysf: ending time [Myr]  
tphys: beginning time [Myr]  

-- bse.bonn

You can compile "bse.bonn", and use it as follows.

1. Move to the "src" directory.

2. Excute the command "make clean" if you already compile bse.geneva or popbin2.geneva.

3. Edit the file "Makefile", changing "MODEL  = GENEVA" to "#MODEL  = GENEVA".

4. Compile "bse.bonn" by executing the command "make bse". Note that
gfortran and g++ compilers are required.

5. Move to the "example/bse.bonn" directory.

6. Execute the command "./bse.bonn".

7. Analyze the output in a similar way to "bse.geneva"

-- popbin2.bonn

You can compile "popbin2.bonn", and use it as follows.

1. Move to the "src" directory.

2. Excute the command "make clean" if you already compile bse.bonn or popbin2.bonn.

3. Edit the file "Makefile", changing "MODEL  = GENEVA" to "#MODEL  = GENEVA".

4. Compile "popbin2.bonn" by executing the command "make
popbin2". Note that gfortran and g++ compilers are required.

5. Move to the "example/popbin2.bonn" directory.

6. Execute the command "./popbin2.bonn".

7. Analyze the output in a similar way to "popbin2.geneva"

Appendix

-- Star type

The star types are the same as those in the original BSE code.

0: Low-mass main-sequence star  
1: high-mass main-sequence star  
2: Hertzsprung gap star  
3: First Giant Branch  
4: Core helium burning star  
5: Early asymptotic giant branch  
6: Thermally pulsing asymptotic giant branch  
7: Naked helium main-sequence star  
8: Naked helium Hertzsprung gap star  
9: Naked helium giant branch  
10: Helium white dwarf  
11: Carbon-oxygen white dwarf  
12: Oxygen-neon white dwarf  
13: Neutron star  
14: Black hole  
15: Massless remnant  

-- Binary status is the same as those in the original BSE code.

The binary status is the 

INITIAL: Initial time  
KW_CHNGE: Changing star type  
BEG_RCHE: Beginning Roche-lobe overflow  
END_RCHE: Ending Roche-lobe overflow  
CONTACT: Contact of two stars  
COELESCE: Coalescence of two stars  
COMENV: Common envelope evolution  
GNTAGE:  
NO_REMNT: No remnant  
MAX_TIME: Ending time  
DISRUPT: Binary disruption  
BEG_SYMB: Beginning symbiotic phase  
END_SYMB: Ending symbiotic phase  
BEG_BSS: Beginning blue straggler phase  

-- Supernova type

NOSN: No supernova  
CCSN: Core-collapse supernova  
DC: Direct collapse to black hole  
PPI: Pulsational pair instability  
PISN: Pair instability supernova  
AIC: Accretion-induced collapse  
