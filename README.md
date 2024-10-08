# How to run CNEO-QM/MM with GROMACS-PySCF

## Download and install modified GROMACS and PySCF:

The source codes of [modified GROMACS](https://github.com/theorychemyang/gromacs) and [modified PySCF](https://github.com/theorychemyang/pyscf) can be downloaded from our GitHub repository, respectively. The procedure for installing PySCF is the same as "Build from source" method on PySCF's [official installation guide](https://pyscf.org/install.html#build-from-source). For GROMACS, the procedure is the same as the [official guide](https://manual.gromacs.org/current/install-guide/), except that when configuring GROMACS with `cmake`, the flag `-DGMX_PYSCF=ON` needs to be set.

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
* Second, the non-bonding parameters for hydrogen atoms on oxygen and nitrogen are often zero in common force fields, but finite values are usually necessary for QM/MM calculations.[^1][^2] It is recommended to modify these parameters [^3] to prevent the polar hydrogen atoms from being exccessively attracted by the MM point charges near the QM/MM boundary.
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
    DFT_E_XC = "B3LYP"
    DFT_DF = True
    QM_E_BASIS_AUX = "aug-cc-pvdz-ri"
    QM_NUC_BASIS = "pb4d"
    MM_CHARGE_MODEL = "point"
    QMMM_CUT = 10

* The methods for treating QM system is selected with `QM_METHOD`, and `CNEO` and `DFT` are available as for now. Additional methods can be added and customized by users. 
* The charge and multiplicity of the QM system are set with the global variables `QM_CHARGE` and `QM_MULT`, respectively. 
* `QM_E_BASIS` selects the electronic basis set, and `DFT_E_XC` selects the electronic exchange-correlation functional. 
* `DFT_DF` controls whether density fitting is to be used for electrons, and `QM_E_BASIS_AUX` sets the auxillary basis if density fitting is enabled.
* For CNEO calculation, `QM_NUC_BASIS` selects protonic basis set.
* `MM_CHARGE_MODEL` determines whether point-charge model (`point`) or a Gaussian-smeared charge model (`gauss`) will be used for MM charges. 
* `QMMM_CUT` determines the range within which from QM region MM charges are to be included in the QM/MM calculation, the unit of this range is Å. In this example, the MM charges within 10 Å from any QM atom will be included in the QM/MM calculation.

#### Select quantum nuclei in CNEO calculations
Two parameters are used to select quantum nuclei in CNEO calculations

    QM_NUC_SELECT = 'all'
    QM_NUC_INDEX = []

`QM_NUC_SELECT` can be eiter `all` or `custom`. In the `all` mode, all hydrogen nuclei except for link atoms (see below) are quantized. In the `custom` mode, the user selects quantum nuclei through option `QM_NUC_INDEX`. It needs attention to create the `QM_NUC_INDEX` as PySCF creates its own index for atoms. `QM_NUC_INDEX` uses this index, not the index in GROMACS.
The procedure to create the list for quantized nuclei in `pyscfdriver.py` is briefly explain here. Before the data production run, the user should first run a 1-step test run to let the PySCF driver produce a list of QM atoms and print their coordinates in `.xyz` format. The user can then use these coordinates to select the hydrogen atoms to be protonated.

### QM parameters related to link atoms in `pyscfdriver.py`
Link atoms are used in our implementation to saturate the QM system when there are covalent bonds crossing the QM/MM boundary. When selecting QM atoms, it is preferable to break non-polar C-C bonds (e.g., between the α-C and β-C of an amino acide residue) to minimize the effect of altering the system with the introduction of link atoms. A link atom is a hydrogon atom that is placed along the broken C-C bond, and its coordinate is solely determined by the coordinates of two atoms forming the broken bond, so no new degrees of freedom are introduced. Here are the explanations for the parameters related to link atoms in `pyscfdriver.py`:

#### Parameters related to determining the coordinates for link atoms

    LINK_COORD_CORR_METHOD = "scale"
    LINK_COORD_SCALE = 0.7246
    LINK_COORD_RFLAT = 1.1

* Option `LINK_COORD_CORR_METHOD` selects how link atoms' coordinates are determined. It can be either `scale` or `flat`. Both methods place link atoms along the broken covalent bonds, but the C-H bond length between the QM carbon atoms and the link atoms are determined differently.
* In `scale` method, the C-H link bond length is `LINK_CORR_SCALE` multiplies the broken bond's length. [^4]
* In `flat` method, the C-H link bond length is always `LINK_CORR_RFLAT` (the unit is Å). [^5]

#### Parameters related to corrections made to MM charges 

    LINK_CHARGE_CORR_METHOD = "global"
    SYSTEM_CHARGE = 0
    LINK_MMHOST_NEIGHBOR_RANGE = 1.7

* Option `LINK_CHARGE_CORR_METHOD` determines how MM charges are modified before QM/MM calculation. The user can choose from `global`, `local`, and `delete`. 
* In the `global` method, the residue charge is spread over the MM atoms that do not form a crossing covalent bond with QM atoms. The residue charge is calculated by adding the classical charge of the QM atoms and MM atoms froming crossing covalent bonds, and then subtracting the `QM_CHARGE`. The charge of QM atoms is itself calculated by subtracting the sum of MM charges from the `SYSTEM_CHARGE`, so this parameter needs to be set correctly. [^5]
* In the `local` method, the charges of MM atoms forming crossing bonds will be spread to the neighbor MM atoms for each of them. The distance cut-off to search for neighbors for a MM atom forming crossing bond is set by `LINK_MMHOST_NEIGHBOR_RANGE` (unit is in Å). [^6]
* The `delete` method simply deletes the charges of the MM atoms forming crossing bonds.

 #### Additional notes about link atoms
 * In GROMACS, all bonds consisting of 2 QM atoms, angles and settles containing 2 or 3 QM atoms, and dihedrals containing 3 or 4 QM atoms, are all excluded from the forcefield evaluation.
 * The force on a link atoms is partitioned to the two atoms of the crossing bond. [^4][^6]


## References

[^1]: [Riccardi, D.; Li, G.; Cui, Q. Importance of van Der Waals Interactions in QM/MM Simulations. J. Phys. Chem. B 2004, 108 (20), 6467–6478.](https://doi.org/10.1021/jp037992q)

[^2]: [Dohn, A. O. Multiscale Electrostatic Embedding Simulations for Modeling Structure and Dynamics of Molecules in Solution: A Tutorial Review. International Journal of Quantum Chemistry 2020, 120 (21), e26343.](https://doi.org/10.1002/qua.26343)

[^3]: [Freindorf, M.; Shao, Y.; Furlani, T. R.; Kong, J. Lennard–Jones Parameters for the Combined QM/MM Method Using the B3LYP/6-31G*/AMBER Potential. Journal of Computational Chemistry 2005, 26 (12), 1270–1278.](https://doi.org/10.1002/jcc.20264)

[^4]: [Eichinger, M.; Tavan, P.; Hutter, J.; Parrinello, M. A Hybrid Method for Solutes in Complex Solvents: Density Functional Theory Combined with Empirical Force Fields. The Journal of Chemical Physics 1999, 110 (21), 10452–10467.](https://doi.org/10.1063/1.479049)

[^5]: [Amber Reference Manual (Chapter 10).](https://ambermd.org/doc12/Amber24.pdf)

[^6]: Sherwood, P. Hybrid Quantum Mechanics/Molecular Mechanics Approaches. In Modern Methods and Algorithms of Quantum Chemistry Proceedings; NIC series; John von Neumann Institute for Computing: Jülich, 2000; pp 285–305

