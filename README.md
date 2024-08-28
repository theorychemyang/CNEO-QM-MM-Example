# How to run CNEO-QM/MM with GROMACS-PySCF

## Download and install modified GROMACS and PySCF:

The source codes of [modified GROMACS](https://github.com/theorychemyang/gromacs) and [modified PySCF](https://github.com/theorychemyang/pyscf) can be downloaded from our GitHub repository, respectively. The procedure for installing PySCF is the same as "Build from source" method on PySCF's [official installation guide](https://pyscf.org/install.html#build-from-source). For GROMACS, the procedure is the same as the [official guide](https://manual.gromacs.org/current/install-guide/), except that when configuring GROMACS with `cmake`, the flag `-DGMX_PYSCF=ON` need to be set:

## Run a CNEO-QM/MM calculation with GROMACS and PySCF:

The procedure of CNEO-QM/MM is similar to that of a generic QM/MM calculation using GROMACS-CP2K packages (https://github.com/bioexcel/gromacs-2022-cp2k-tutorial). The required user inputs in a working directory are: 
1. initial configuration (`.gro` file)
2. force field parameters (`.top`, `.itp`, `.rtp` files, etc.)
3. calculation parameters (`.mdp` file)
4. atom group indices (used for selecting QM atoms and other user-defined groups, `.ndx` file)
5. a PySCF driver `pyscfdriver.py`

### Prepare the initial configuration and force field parameters
For a certain system, the initial configuration and the force field parameters for QM/MM calculation are similar to those for pure MM calculation, but there are a few aspects that need attention. 
* First, it is better to correctly label the atoms with `gmx editconf` command when preparing the initial configuration, because the atom indices in the initial configuration will be used to select QM atoms in the '.ndx' file later. 
* Second, the non-bonding parameters for hydrogen atoms on oxygen and nitrogen are often zero in common force fields, but finite values are usually necessary for QM/MM calculations.[[1]](#1)[[2]](#2) It is recommended to modify these parameters [[3]](#3) to prevent the polar hydrogen atoms from being exccessively attracted by the MM point charges near the QM/MM boundary.
* Third, it is also better to avoid imposing too many constraints on the QM system in the force field, therefore, the constraints on C-H, N-H, or O-H bonds in the QM system should be disabled in the force field.

### Enable QM/MM and select QM atoms in `.mpd` and `.ndx` files
With the initial configuration and force field parameters prepared, the next step is to enable the QM/MM calculation and select QM atoms with the options in `.mdp` file. Here is an example:

    qmmm-pyscf-active  = True
    qmmm-pyscf-qmgroup = QMAtoms

The QM/MM calculation is enabled by setting the `qmmm-pyscf-active` option to `Ture`. QM atoms are selected by first creating a group for QM atom indices in the `.ndx` file, and then setting the option `qmmm-pyscf-qmgroup` to refer to this group. In this example, the `[ QMAtoms ]` group are defined in the `.ndx` file and then used here to select QM atoms.

### Set QM calculation parameters in `pyscfdriver.py`
The PySCF driver can be copied from `gromacs_dir/src/gromacs/applied_forces/qmmm/pyscfdriver.py` and must be present in the working directory. The parameters for the QM calculation are set with global variables in `pyscfdriver.py` and are explained with the following example:

    QM_METHOD = "cneo"
    QM_CHARGE = -1
    QM_MULT = 1
    QM_E_BASIS = "aug-cc-pvdz"
    DFT_ELE_XC = "B3LYP"
    DFT_DF_FIT = True
    QM_E_BASIS_AUX = "aug-cc-pvdz-ri"
    QM_NUC_BASIS = "pb4d"
    MM_CHARGE_MODEL = "point"
    QMMM_CUT = 10

    ..* The methods for treating QM system is selected with `QM_METHOD`, and `CNEO` and `DFT` are available as for now. Additional methods can be added and customized by users. 
    
* The charge and multiplicity of the QM system are set with the global variables `QM_CHARGE` and `QM_MULT`, respectively. 
* `QM_E_BASIS` selects the electronic basis set, and `DFT_ELE_XC` selects the electronic exchange-correlation functional. 
* `DFT_DF_FIT` controls whether density fitting is to be used for electrons, and `QM_E_BASIS_AUX` sets the auxillary basis if density fitting is enabled.
* For CNEO calculation, `QM_NUC_BASIS` selects protonic basis set. `MM_CHARGE_MODEL` determines whether point-charge model (`point`) or a Gaussian-smeared charge model (`Gauss`) will be used for MM charges. 
* `QMMM_CUT` determines the range within which from QM region MM charges are to be included in the QM/MM calculation, the unit of this range is Å. In this example, the MM charges within 10 Å from any QM atom will be included in the QM/MM calculation.
## References

<a id="1">[1]</a> Riccardi, D.; Li, G.; Cui, Q. Importance of van Der Waals Interactions in QM/MM Simulations. J. Phys. Chem. B 2004, 108 (20), 6467–6478. https://doi.org/10.1021/jp037992q.

<a id="2">[2]</a> Dohn, A. O. Multiscale Electrostatic Embedding Simulations for Modeling Structure and Dynamics of Molecules in Solution: A Tutorial Review. International Journal of Quantum Chemistry 2020, 120 (21), e26343. https://doi.org/10.1002/qua.26343.

<a id="3">[3]</a> 
Freindorf, M.; Shao, Y.; Furlani, T. R.; Kong, J. Lennard–Jones Parameters for the Combined QM/MM Method Using the B3LYP/6-31G*/AMBER Potential. Journal of Computational Chemistry 2005, 26 (12), 1270–1278. https://doi.org/10.1002/jcc.20264.
