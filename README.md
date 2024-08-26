# CNEO QM/MM Example
Example of how to run CNEO-QM/MM with GROMACS and PySCF

The constrained nuclear-electronic orbital (CNEO) method treats nuclei quantum mechanically by applying constraints to the expectation of nuclear position operator. CNEO QM/MM is a QM/MM method that uses the CNEO to treat the QM region of the QM/MM system, so that the nuclei in the QM region can be treated quantum mechanically. CNEO QM/MM is implemented in a modified version of GROMACS (https://github.com/theorychemyang/gromacs/tree/pyscf_testii) that uses modified PySCF package with CNEO capability (https://github.com/theorychemyang/pyscf) to treat the QM region.

## Installing modified PySCF 

Prerequisites for manual install are 
* C compiler 
* C++ compiler(optional, but required for XCFun and some extensions)
* CMake >= 3.10
* Python >= 3.7
* Numpy >= 1.13
* Scipy >= 0.19
* h5py >= 2.7 

The installation is the same as the build from the source method for the official PySCF package (https://pyscf.org/install.html), except downloading the locally modified version from our group GitHub:  

 	git clone https://github.com/theorychemyang/pyscf.git 

 	cd pyscf 
  
Next, build the C extensions in pyscf/lib: 

	cd pyscf/lib

	mkdir build 
	
	cd build 
 
	cmake .. 

	make 
Finally, to allow Python to find the PySCF package, add the top-level PySCF directory (not the pyscf/pyscf subdirectory) to `PYTHONPATH.`: 

	export PYTHONPATH=/path/to/pyscf:$PYTHONPATH 

## Install modified GROMACS
The prerequisites to install the modified version of GROMACS are the same as the official GROMACS package (https://manual.gromacs.org/documentation/current/install-guide/), except that the source code needs to be downloaded from our GitHub repository:
	git clone https://github.com/theorychemyang/gromacs.git

 	cd gromacs
	
	git checkout pyscf_testii
The next step is to compile the modified GROMACS. The `-DGMX_PYSCF` flag needs to be turned on along with other flags required to perform GROMACS QM/MM calculations:
	
 	mkdir -p build/
	
 	cd build/
 	
  	cmake .. -DBUILD_SHARED_LIBS=OFF -DGMXAPI=OFF -DGMX_INSTALL_NBLIB_API=OFF -DGMX_DOUBLE=ON -DGMX_BUILD_OWN_FFTW=ON -DGMX_PYSCF=ON
	
 	make
## Run a CNEO QM/MM calculation with GROMACS and PySCF:
The procedure of CNEO QM/MM is similar to that of a generic QM/MM calculation using GROMACS-CP2K packages (https://github.com/bioexcel/gromacs-2022-cp2k-tutorial). The required user inputs in a working directory are: (1) initial configuration (`.gro` file); (2) force field parameters (`.top`, `.itp`, `.rtp` files, etc.); (3) calculation parameters (`.mdp` file); (4) atom group indices (used for selecting QM atoms and other user-defined groups, `.ndx` file). With these files prepared in the working directory, a portable run input file (`.tpr`) needs to be generated with grompp command.
	
 	gmx_d grompp -f md_params.mdp -c conf.gro -n index.ndx -p topol.top -o md.tpr

For running QM/MM calculation with GROMACS-PySCF packages, a PySCF driver `pyscfdriverii.py` needs to be present in the working directory, which is called by the modified version of GROMACS during the calculation. The parameters can be set with the global variables at the beginning of this driver. With this driver set up, the next step is to run the 
	
 	gmx_d mdrun -deffnm md

