import numpy as np
import time, copy, os, warnings
import resource
import json
from functools import reduce
from operator import index as integer_index
try:
    import cupy as cp
    HAVE_GPU = True
except ImportError:
    HAVE_GPU = False
from pyscf import lib, gto, scf, grad, dft, neo, lo
from pyscf.qmmm import itrf
from pyscf.data import elements, radii, nist
from pyscf.qmmm.mm_mole import create_mm_mol
from collections import defaultdict
# Enable SCF-guess prediction only for fixed-timestep MD:
# True selects projected G-Ext with four snapshots; 'aspc' selects ASPC.
# Projected G-Ext fits historical AO overlaps and applies the fitted weights
# to occupied-space projectors. Method: Wilhelm et al., Sec. X.2.1,
# Eq. 100, doi:10.48550/arXiv.2607.22916. Implementation follows CP2K's
# wfi_gext_proj_nr in src/qs_wf_history_methods.F at revision
# dfe8787b4f0a22da8e73e9abe2b384fa021d2f53:
# https://github.com/cp2k/cp2k/blob/dfe8787b4f0a22da8e73e9abe2b384fa021d2f53/src/qs_wf_history_methods.F
# Four snapshots are a benchmark-selected default for this implementation.
# A dict exposes predictor options, e.g. dict(scheme='aspc', history_length=6).
# Leave this False for geometry optimization. GROMACS keeps history in this
# imported module; restarting GROMACS starts a fresh history in memory.
MD_EXTRAPOLATION = False
USE_GPU = False
SCF_CONV_TOL = 1e-9
SCF_CONV_TOL_GRAD = None
SCF_MAX_CYCLE = 50
# Print one timing record per force evaluation. Set a JSONL path to also
# append records to disk; this file contains diagnostics only.
PRINT_TIMING = True
TIMING_FILE = None
_extrapolator = None
_extrapolation_options = None


def reset_extrapolation():
    """Start a new MD prediction history after restarting or changing runs."""
    global _extrapolator, _extrapolation_options
    _extrapolator = None
    _extrapolation_options = None


def run_scf_extrapolation(mf, step):
    """Apply driver settings and manage the history retained between GROMACS calls."""
    global _extrapolator, _extrapolation_options
    # GROMACS retains this imported module across force calls. Each call
    # supplies a fresh SCF object with the current QM and MM coordinates.
    if isinstance(step, (bool, np.bool_)):
        raise ValueError('gmxstep must be an integer MD step number')
    step = integer_index(step)
    option = MD_EXTRAPOLATION
    if option is None or option is False:
        options = None
    elif option is True:
        options = {}
    elif isinstance(option, str):
        options = dict(scheme=option.lower())
    elif isinstance(option, dict):
        options = option.copy()
    else:
        raise ValueError('MD_EXTRAPOLATION must be False, True, a method name, or a dict')
    mf.conv_tol = SCF_CONV_TOL
    mf.conv_tol_grad = SCF_CONV_TOL_GRAD
    mf.max_cycle = SCF_MAX_CYCLE
    if DFT_DF_NN and not (DFT_DF and DFT_DF_NE):
        raise ValueError('DFT_DF_NN requires DFT_DF=True and DFT_DF_NE=True')
    try:
        scf_start = wall_time()
        if options is None or step < 0:
            reset_extrapolation()
            energy = mf.kernel()
        else:
            from pyscf.neo.extrapolation import Extrapolator
            if _extrapolator is None or options != _extrapolation_options:
                _extrapolator = Extrapolator(**options)
                _extrapolation_options = options
            last = _extrapolator.last_step
            if last is not None and (step < last or step > last+1):
                _extrapolator.clear()
            if (step == _extrapolator.last_step and
                not np.array_equal(mf.mol.atom_coords(), _extrapolator.last_coords)):
                warnings.warn('QM coordinates changed within one MD step; restarting SCF history',
                              RuntimeWarning)
                _extrapolator.clear()
            dm0 = _extrapolator.predict(mf, mf.mol, step)
            with _extrapolator.record_fock(mf):
                energy = mf.kernel(dm0=dm0)
        scf_seconds = wall_time() - scf_start
        energy = float(energy)
        if not mf.converged or not np.isfinite(energy):
            raise RuntimeError(f'SCF failed at GROMACS step {step}; forces were not evaluated')
        if _extrapolator is not None:
            _extrapolator.update(mf, step)
    except Exception:
        reset_extrapolation()
        raise
    return energy, scf_seconds


def wall_time():
    # GPU kernels finish asynchronously; synchronize the current stream at
    # section boundaries so wall times include the completed calculation.
    if USE_GPU:
        cp.cuda.get_current_stream().synchronize()
    return time.perf_counter()


def record_timing(mf, step, scf_seconds, grad_seconds):
    """Record one completed SCF and its QM/MM force evaluation, in seconds."""
    record = dict(step=integer_index(step), method=QM_METHOD,
                  backend='GPU' if USE_GPU else 'CPU', energy=float(mf.e_tot),
                  converged=bool(mf.converged), scf_cycles=int(mf.cycles),
                  scf_seconds=scf_seconds, grad_seconds=grad_seconds)
    line = json.dumps(record, allow_nan=False)
    if PRINT_TIMING:
        print('SCF/gradient: ' + line, flush=True)
    if TIMING_FILE is not None:
        with open(TIMING_FILE, 'a') as output:
            output.write(line + '\n')


def prepare_qmatoms(kinds, coords, indices=None):
    """Build an independent atom list and select zero-based QM nuclear indices."""
    if len(kinds) != len(coords):
        raise ValueError('QM element and coordinate counts differ')
    coordinates = np.asarray(coords, dtype=float)
    if coordinates.shape != (len(kinds), 3) or not np.isfinite(coordinates).all():
        raise ValueError('QM coordinates must be a finite (number of QM atoms, 3) array in Angstrom')
    qmatoms = [[kind.strip(), *coord.tolist()] for kind, coord in zip(kinds, coordinates)]
    if QM_METHOD.lower() != 'cneo':
        return qmatoms, []
    selection = QM_NUC_SELECT.lower()
    if selection == 'all':
        selected = [i for i, kind in enumerate(kinds)
                    if kind.strip().upper() == 'H' and
                    (indices is None or not str(indices[i]).startswith('L'))]
    elif selection == 'custom':
        if any(isinstance(i, (bool, np.bool_)) for i in QM_NUC_INDEX):
            raise ValueError('QM_NUC_INDEX must contain integer QM-local indices')
        selected = [integer_index(i) for i in QM_NUC_INDEX]
        if len(set(selected)) != len(selected):
            raise ValueError('QM_NUC_INDEX contains duplicate indices')
        for i in selected:
            if i < 0 or i >= len(kinds) or kinds[i].strip().upper() != 'H':
                raise ValueError('QM_NUC_INDEX must contain zero-based QM-local hydrogen indices')
            if indices is not None and str(indices[i]).startswith('L'):
                raise ValueError('Link atoms must remain classical nuclei')
    else:
        raise ValueError("QM_NUC_SELECT must be 'all' or 'custom'")
    for i in selected:
        qmatoms[i][0] = 'H~'
    return qmatoms, selected


def force_array(value, count):
    """Return the contiguous float64 CPU force buffer read by GROMACS."""
    if hasattr(value, 'get'):
        value = value.get()
    value = np.asarray(value, dtype=np.float64, order='C')
    if value.shape != (count, 3) or not np.isfinite(value).all():
        raise RuntimeError(f'Expected finite forces with shape ({count}, 3), got {value.shape}')
    return value

# DEBUG = True
SYSTEM_CHARGE = 0
QM_CHARGE = -1
QM_MULT = 1
QM_METHOD = 'cneo'
# select from {'cneo', 'dft', 'export': for export qm and mm xyz files }
# export options: uses EXPORT_QM_FNAME and EXPORT_MM_FNAME to specify the output file names, optional
QM_E_BASIS = 'def2-svp'
DFT_E_XC = 'B3LYP'
QM_NUC_BASIS = 'pb4d'
QM_NUC_SELECT = 'all' # select from {'all', 'custom'}
QM_NUC_INDEX = []
DFT_DF = True
DFT_DF_NE = True
DFT_DF_NN = True
QM_E_BASIS_AUX = 'def2-svp-jkfit'
MM_CHARGE_MODEL = 'point' # select from {'point', 'gaussian'}
QMMM_CUT = None # None includes all supplied MM charges; finite distances are in Angstrom
QMMM_SWITCH = 3 # Start tapering MM charges here when QMMM_CUT is finite

LINK_CHARGE_CORR_METHOD = 'global' # select from {'global', 'local', 'delete'}
LINK_MMHOST_NEIGHBOR_RANGE = 1.7
LINK_COORD_CORR_METHOD = 'scale' # select from {'scale', 'flat'}
LINK_COORD_SCALE = 0.7246
LINK_COORD_RFLAT = 1.1
LINK_PRINT_COORD_CORR = False
LINK_PRINT_FORCE_CORR = False
LINK_PRINT_CHARGE_CORR = False
QMMM_PRINT = False
QM_XYZ = False
EXPORT_QM_FNAME = 'qm.xyz'
EXPORT_MM_FNAME = 'mm.xyz'

def print_memory_usage(step, point):
    # CPU memory
    cpu_memory_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024

    memory_str = f"Step {step}, {point}:\n  Peak CPU Memory: {cpu_memory_mb:.1f} MB"

    # GPU memory if available
    if HAVE_GPU:
        try:
            gpu_memory_mb = cp.get_default_memory_pool().used_bytes() / 1024 / 1024
            gpu_cached_mb = cp.get_default_memory_pool().total_bytes() / 1024 / 1024
            memory_str += f"\n  GPU Memory: {gpu_memory_mb:.1f} MB (Cached: {gpu_cached_mb:.1f} MB)"

        except:
            memory_str += "\n  GPU memory query failed"

    print(memory_str)


def qmmmCalc(
    gmxstep,
    qmbasis,
    qmmult,
    qmcharge,
    qmindex,
    qmkinds,
    qmcoords,
    mmindex,
    mmkinds,
    mmcharges,
    mmcoords,
    links,
):
    if QMMM_CUT is not None:
        if not np.isfinite(QMMM_CUT) or not 0 <= QMMM_SWITCH < QMMM_CUT:
            raise ValueError('Require 0 <= QMMM_SWITCH < finite QMMM_CUT, or QMMM_CUT=None')
        warnings.warn('Finite QMMM_CUT scales MM charges with distance; forces omit '
                      'derivatives of this scaling. Use QMMM_CUT=None for MD.', RuntimeWarning)
    if MM_CHARGE_MODEL.lower() != 'point' and not MM_CHARGE_MODEL.lower().startswith('gauss'):
        raise ValueError("MM_CHARGE_MODEL must be 'point' or 'gaussian'")
    if len(qmindex) != len(qmkinds) or len(qmkinds) != len(qmcoords):
        raise ValueError('QM indices, elements and coordinate counts differ')
    if not len(mmindex) == len(mmkinds) == len(mmcharges) == len(mmcoords):
        raise ValueError('MM indices, elements, charges and coordinate counts differ')
    print(f"gmxstep is {gmxstep}")
    print_memory_usage(gmxstep, "start of qmmmCalc")
    t0 = time.time()
    qmindex_pyscf = [x for x in range(len(qmkinds))]
    if gmxstep == 0 or gmxstep == -1:
        prop_print_xzy("QM atoms", qmindex_pyscf, qmkinds, qmcoords)
    # Generate link atoms coordinates, and
    # extend kinds, coords, and index lists of qm atoms. Default mode:

    [qmkinds_link, qmcoords_link, qmindex_link, link_scale] = link_coord_corr(
        links,
        qmindex,
        qmkinds,
        qmcoords,
        mmindex,
        mmkinds,
        mmcoords,
        printflag=LINK_PRINT_COORD_CORR,
        method=LINK_COORD_CORR_METHOD,
        rflat0=LINK_COORD_RFLAT,
        scale0=LINK_COORD_SCALE,
    )

    qmatoms, qmnucindex = prepare_qmatoms(qmkinds_link, qmcoords_link, qmindex_link)
    # Generate the list of quantum nuclei for CNEO or NEO calculation.
    # The last (# of links H atoms) are link atoms, therefore should be
    # removed from the list of quantum nuclei for CNEO calculation

    # Redistribute the classical QM and MM link host charges.
    # Default mode: spread the residue classical charge over the MM
    # atoms less MM link host atoms, i.e., QM atoms and MM link host
    # atoms are 'zeroed-out'. This is the only mode currently implemented
    mmcharges = link_charge_corr(
            mmcharges=mmcharges,
            mmindex=mmindex,
            links=links,
            qm_charge=QM_CHARGE,
            system_charge=SYSTEM_CHARGE,
            printflag=LINK_PRINT_CHARGE_CORR,
            charge_corr_mode=LINK_CHARGE_CORR_METHOD,
            mmcoords=mmcoords,
            mmkinds=mmkinds,
            mmneighbor_thrsh=LINK_MMHOST_NEIGHBOR_RANGE)

    if QMMM_CUT is None:
        mmdist = {i: [None, kind, None, 1., charge]
                  for i, (kind, charge) in enumerate(zip(mmkinds, mmcharges))}
    else:
        mmdist = defaultdict(list)
        for qmindex_, qmkind, qmcoord in zip(qmindex, qmkinds, qmcoords):
            # print(f"{qmindex_=} {qmkind=} {qmcoord=}")
            for mmindex_, mmcoord in enumerate(mmcoords):
                # print(f"{mmindex_=} {mmcoord=}")
                dist = lib.norm(qmcoord - np.array(mmcoord))
                if dist <= QMMM_CUT:
                    switch_factor = switch(dist, QMMM_SWITCH, QMMM_CUT)
                    charge_switched = mmcharges[mmindex_] * switch_factor
                    # charge_old = mmdist[mmindex_][2]
                    mmkind = mmkinds[mmindex_]
                    if mmindex_ not in mmdist:
                        mmdist[mmindex_] = [qmindex_, mmkind, dist, switch_factor, charge_switched]
                    elif mmdist[mmindex_][2] > dist:
                        mmdist[mmindex_] = [qmindex_, mmkind, dist, switch_factor, charge_switched]
    # mmindex_incut = np.array(sorted(list(set(mmdist.keys()))))
    mmkinds_incut, mmcharges_incut, mmcoords_incut, mmindex_incut = [], [], [], []
    for i, mmindex_ in enumerate(sorted(mmdist.keys())):
        # print(f"{i=} {mmindex_=} {mmdist[mmindex_]=}")
        qmindex_, mmkind, dist, switch_factor, charge_switched = mmdist[mmindex_]
        # print(f'{qmindex_=} {mmindex_=} {mmkind=} {dist:.3f} {switch_factor:.3f} old charge {mmcharges[mmindex_]:.3f} new charge {charge_switched:.3f}')
        mmkinds_incut.append(mmkind)
        mmcharges_incut.append(charge_switched)
        mmcoords_incut.append(mmcoords[mmindex_])
        mmindex_incut.append(mmindex_)
    mmkinds_incut = np.array(mmkinds_incut)
    mmcharges_incut = np.array(mmcharges_incut)
    mmcoords_incut = np.asarray(mmcoords_incut, dtype=float).reshape(-1, 3)
    mmindex_incut = np.asarray(mmindex_incut, dtype=int)

    # mmkinds_incut = np.array(mmkinds)[mmindex_incut]
    # mmcharges_incut = np.array(mmcharges)[mmindex_incut]
    # mmcoords_incut = np.array(mmcoords)[mmindex_incut]

    # PySCF selects point charges when radii=None; supplying radii selects
    # normalized Gaussian charge distributions.
    mmradii = None
    if MM_CHARGE_MODEL.lower().startswith('gauss'):
        mmradii = [radii.COVALENT[elements.charge(kind)] * nist.BOHR
                   for kind in mmkinds_incut]

    if QMMM_PRINT:
        print("qm nuc index:", QM_NUC_INDEX)
        print("mmcharges\n", mmcharges)
        print("mmcoords\n", mmcoords)
        print("qmatomds\n", qmatoms)
        print(f"method we used in this calc is {QM_METHOD}")
        print(mmradii)

    if QM_XYZ:
            with open('qm.xyz', 'w') as f:
                f.write(str(len(qmatoms))+f"\n QM xyz at step {gmxstep}\n")
                for i in range(len(qmatoms)):
                    f.write(f"{qmatoms[i][0]} {qmatoms[i][1]} {qmatoms[i][2]} {qmatoms[i][3]}\n")

    if QM_METHOD.upper() == 'CNEO':
        [energy, qmforces, mmforces_incut] = qmmmCalc_cneo(
            qmatoms, mmcoords_incut, mmcharges_incut, mmradii, qmnucindex=qmnucindex,
            gmxstep=gmxstep
        )
    elif QM_METHOD.upper() == 'DFT':
        [energy, qmforces, mmforces_incut] = qmmmCalc_dft(
            qmatoms, mmcoords_incut, mmcharges_incut, mmradii, gmxstep=gmxstep
        )
    elif QM_METHOD.upper() == 'MULLIKEN':
        [energy, qmforces, mmforces_incut] = qmmmCalc_mulliken(
            qmatoms, mmcoords_incut, mmcharges_incut, mmradii
        )

    elif QM_METHOD.upper() == 'DFTB':
        if os.path.isfile('qm.xyz'):
            pass
        else:
            with open('qm.xyz', 'w') as f:
                f.write(str(len(qmatoms))+"\nfake QM xyz\n")
                for i in range(len(qmatoms)):
                    f.write(f"{qmatoms[i][0]} {0.0} {0.0} {0.0}\n")
        [energy, qmforces, mmforces_incut] = qmmmCalc_dftb(
            qmcoords_link, mmcoords_incut, mmcharges_incut)
    elif QM_METHOD.upper() == 'EXPORT':
        try:
            qmCalc_export(qmatoms, mmkinds_incut, mmcoords_incut, mmcharges_incut,
                          export_qm_fname=EXPORT_QM_FNAME, export_mm_fname=EXPORT_MM_FNAME)
        except:
            qmCalc_export(qmatoms, mmkinds_incut, mmcoords_incut, mmcharges_incut)
        [energy, qmforces, mmforces_incut] = qmmmCalc_dftb(
            qmcoords_link, mmcoords_incut, mmcharges_incut)
    else:
        raise ValueError(f'Unknown QM_METHOD: {QM_METHOD}')

    qmforces = force_array(qmforces, len(qmatoms))
    mmforces_incut = force_array(mmforces_incut, len(mmindex_incut))
    mmforces = np.zeros((len(mmindex), 3))
    mmforces[mmindex_incut] = mmforces_incut

    # force correction, partition the force on link atoms to
    # its host atoms with chain rule
    [qmindex_force_corr, qmforces_link, mmforces_link] = link_force_corr(
        links,
        qmindex_link,
        qmkinds_link,
        qmforces,
        mmindex,
        mmkinds,
        mmforces,
        link_scale,
        printflag=LINK_PRINT_FORCE_CORR,
    )

    print(f"time for this step = {time.time() - t0} seconds")
    print_memory_usage(gmxstep, "end of qmmmCalc")

    return energy, qmforces_link, mmforces_link


def qmCalc(gmxstep, qmbasis, qmmult, qmcharge, qmkinds, qmcoords):
    print_memory_usage(gmxstep, "start of qmCalc")
    t0 = time.time()
    qmindex_pyscf = [x for x in range(len(qmkinds))]
    if gmxstep == 0 or gmxstep == -1:
        prop_print_xzy("QM atoms", qmindex_pyscf, qmkinds, qmcoords)

    qmatoms, qmnucindex = prepare_qmatoms(qmkinds, qmcoords)

    if QM_METHOD.upper() == 'CNEO':
        [energy, qmforces] = qmCalc_cneo(qmatoms, qmnucindex=qmnucindex, gmxstep=gmxstep)
    elif QM_METHOD.upper() == 'DFT':
        [energy, qmforces] = qmCalc_dft(qmatoms, gmxstep=gmxstep)
    else:
        raise ValueError('Pure-QM calculations support QM_METHOD cneo or dft')

    print(f"time for this step = {time.time() - t0} seconds")
    print_memory_usage(gmxstep, "end of qmCalc")

    return energy, force_array(qmforces, len(qmatoms))


def qmCalc_cneo(qmatoms, qmnucindex, gmxstep=0):
    mol = neo.M(
        atom=qmatoms,
        basis=QM_E_BASIS,
        nuc_basis=QM_NUC_BASIS,
        quantum_nuc=qmnucindex,
        charge=QM_CHARGE,
        spin=QM_MULT-1
    )
    mf = neo.CDFT(mol, xc=DFT_E_XC)
    if DFT_DF:
        mf = mf.density_fit(auxbasis=QM_E_BASIS_AUX, df_ne=DFT_DF_NE, df_nn=DFT_DF_NN)
    if USE_GPU:
        mf = mf.to_gpu()
    # mf.mf_elec.xc = DFT_E_XC

    energy, scf_seconds = run_scf_extrapolation(mf, gmxstep)
    grad_start = wall_time()

    g = mf.Gradients()
    g_qm = g.grad()
    qmforce = -g_qm
    record_timing(mf, gmxstep, scf_seconds, wall_time() - grad_start)
    return energy, qmforce


def qmCalc_dft(qmatoms, gmxstep=0):
    mol = gto.M(atom=qmatoms,
                unit='ANG',
                basis=QM_E_BASIS,
                charge=QM_CHARGE,
                spin=QM_MULT-1)
    mf = dft.UKS(mol) if mol.spin else dft.RKS(mol)
    mf.xc = DFT_E_XC
    if DFT_DF:
        mf = mf.density_fit(auxbasis=QM_E_BASIS_AUX)
    if USE_GPU:
        mf = mf.to_gpu()

    energy, scf_seconds = run_scf_extrapolation(mf, gmxstep)
    grad_start = wall_time()
    # dm = mf.make_rdm1()

    qmgrad = mf.nuc_grad_method().kernel()
    qmforce = -qmgrad
    record_timing(mf, gmxstep, scf_seconds, wall_time() - grad_start)
    return energy, qmforce


def qmmmCalc_cneo(qmatoms, mmcoords, mmcharges, mmradii, qmnucindex, gmxstep=0):

    t0 = time.time()
    if len(mmcharges) == 0:
        energy, qmforces = qmCalc_cneo(qmatoms, qmnucindex, gmxstep)
        return energy, qmforces, np.empty((0, 3))

    mmmol = create_mm_mol(mmcoords, mmcharges, mmradii)
    mol_neo = neo.M(
        atom=qmatoms,
        basis=QM_E_BASIS,
        nuc_basis=QM_NUC_BASIS,
        # If there is no crossing covalent bond and the user
        # wishes to quantizd all H atoms, ['H'] can be used instead
        quantum_nuc=qmnucindex,
        mm_mol=mmmol,
        charge=QM_CHARGE,
        spin=QM_MULT-1
    )
    # energy
    print("mol_neo quantum_nuc", mol_neo._quantum_nuc)
    # print(qmatoms)
    mf = neo.CDFT(mol_neo, xc=DFT_E_XC)
    if DFT_DF:
        mf = mf.density_fit(auxbasis=QM_E_BASIS_AUX, df_ne=DFT_DF_NE, df_nn=DFT_DF_NN)
    if USE_GPU:
        mf = mf.to_gpu()
    # mf.mf_elec.xc = DFT_E_XC
    energy, scf_seconds = run_scf_extrapolation(mf, gmxstep)
    grad_start = wall_time()
    te = time.time()
    print(f"time for energy = {te - t0} seconds")

    # qm gradient
    g = mf.Gradients()
    g_qm = g.grad()
    qmforces = -g_qm
    tqf = time.time()
    print(f"time for qm force = {tqf - te} seconds")

    # mm gradient
    g_mm = g.grad_mm()
    mmforces = -g_mm
    record_timing(mf, gmxstep, scf_seconds, wall_time() - grad_start)
    tmf = time.time()
    print(f"time for mm force = {tmf - tqf} seconds")

    return energy, qmforces, mmforces


def qmmmCalc_dft(qmatoms, mmcoords, mmcharges, mmradii, gmxstep=0):
    t0 = time.time()
    if len(mmcharges) == 0:
        energy, qmforces = qmCalc_dft(qmatoms, gmxstep)
        return energy, qmforces, np.empty((0, 3))

    mol = gto.M(atom=qmatoms,
                unit='ANG',
                basis=QM_E_BASIS,
                charge=QM_CHARGE,
                spin=QM_MULT-1)
    mf_dft = dft.UKS(mol) if mol.spin else dft.RKS(mol)
    mf_dft.xc = DFT_E_XC

    if DFT_DF:
        mf_dft = mf_dft.density_fit(auxbasis=QM_E_BASIS_AUX)
    if USE_GPU:
        if MM_CHARGE_MODEL.lower() != 'point':
            raise NotImplementedError('Conventional GPU DFT QM/MM currently requires point charges')
        from gpu4pyscf.qmmm import itrf as gpu_itrf
        mf = gpu_itrf.mm_charge(mf_dft.to_gpu(), mmcoords, mmcharges)
    else:
        mf = itrf.mm_charge(mf_dft, mmcoords, mmcharges,
                            None if MM_CHARGE_MODEL.lower() == 'point' else mmradii)

    energy, scf_seconds = run_scf_extrapolation(mf, gmxstep)
    grad_start = wall_time()
    te = time.time()
    # print(f"time for energy = {te - t0} seconds")

    dm = mf.make_rdm1()
    mf_grad = mf.Gradients()
    qmforces = -mf_grad.kernel()
    tqf = time.time()
    # print(f"time for qm force = {tqf - te} seconds")

    mmforces_qmnuc = -mf_grad.grad_nuc_mm()
    mmforces_qme = -mf_grad.grad_hcore_mm(dm)
    mmforces = mmforces_qmnuc + mmforces_qme
    record_timing(mf, gmxstep, scf_seconds, wall_time() - grad_start)
    tmf = time.time()
    # print(f"time for mm force = {tmf - tqf} seconds")
    # print(f"time for this step = {time.time() - t0} seconds")

    return energy, qmforces, mmforces

def qmmmCalc_mulliken(qmatoms, mmcoords, mmcharges, mmradii):
    t0 = time.time()

    mol = gto.M(atom=qmatoms, unit="ANG", basis=QM_E_BASIS, charge=QM_CHARGE, spin=QM_MULT-1)
    mf_dft = dft.RKS(mol)
    mf_dft.xc = DFT_E_XC
    mf = itrf.mm_charge(mf_dft, mmcoords, mmcharges, mmradii)

    mf.kernel()

    energy = 0.000
    qmforces = np.array([[0.00,0.00,0.00] for x in qmatoms])
    mmforces = np.array([[0.00,0.00,0.00] for x in mmcharges])

    C = lo.orth_ao(mf, 'nao')

    # C is orthogonal wrt to the AO overlap matrix.  C^T S C  is an identity matrix.
    print(abs(reduce(np.dot, (C.T, mf.get_ovlp(), C)) -
            np.eye(mol.nao_nr())).max())  # should be close to 0

    # The following linear equation can also be solved using the matrix
    # multiplication reduce(np.dot (C.T, mf.get_ovlp(), mf.mo_coeff))
    mo = np.linalg.solve(C, mf.mo_coeff)

    #
    # Mulliken population analysis based on NAOs
    #
    dm = mf.make_rdm1(mo, mf.mo_occ)
    mf.mulliken_pop(mol, dm, np.eye(mol.nao_nr()))

    return energy, qmforces, mmforces

def qmmmCalc_dftb(qmcoords, mmcoords, mmcharges):
    import dftbplus
    qmcoords = np.array(qmcoords) / 0.529177249 # AA to Bohr
    mmcoords = np.array(mmcoords) / 0.529177249 # AA to Bohr
    dr = mmcoords[:,None,:] - qmcoords
    r = np.linalg.norm(dr, axis=2)
    extpot = np.einsum('R,Rr->r', mmcharges, 1/r)

    cdftb = dftbplus.DftbPlus(libpath='path to /libdftbplus.so',
                            hsdpath='dftb_in.hsd',
                            logfile='dftb_log.log')
    cdftb.set_geometry(qmcoords, latvecs=None)
    cdftb.set_external_potential(extpot, extpotgrad=np.zeros((qmcoords.shape[0], 3)))
    energy = cdftb.get_energy()
    qmforces = -cdftb.get_gradients()
    qmcharges = cdftb.get_gross_charges()
    cdftb.close()

    qmforces += -np.einsum('r,R,Rrx,Rr->rx', qmcharges, mmcharges, dr, r**-3)
    mmforces = -np.einsum('r,R,Rrx,Rr->Rx', qmcharges, mmcharges, dr, r**-3)
    return energy, qmforces, mmforces


def prop_print_xzy(xyzpropname, index, kinds, xyzprops):
    print("--------------------", xyzpropname, "--------------------")
    for i in range(len(xyzprops)):
        print(
            "%4s %4s %18.12f %18.12f %18.12f"
            % (index[i], kinds[i], xyzprops[i][0], xyzprops[i][1], xyzprops[i][2])
        )

def qmCalc_export(qmatoms, mmkinds, mmcoords, mmcharges,
                  export_qm_fname='qm.xyz', export_mm_fname='mm.xyz', **kwargs):
    # later to add export file type selection {gro, xyz, pdb, ...}
    fqm = open(export_qm_fname, 'w')
    fqm.write(str(len(qmatoms))+f"\nqm symbol {'x':14s}{'y':14s}{'z':14s}\n")
    for i in range(len(qmatoms)):
        fqm.write(f"{qmatoms[i][0]:4s} {qmatoms[i][1]:14.8f} {qmatoms[i][2]:14.8f} {qmatoms[i][3]:14.8f}\n")
    fmm = open(export_mm_fname, 'w')
    fmm.write(str(len(mmcoords))+f"\nmm symbol {'x':14s}{'y':14s}{'z':14s}{'charge':9s}\n")
    for i in range(len(mmcoords)):
        fmm.write(f"{mmkinds[i]:4s} {mmcoords[i][0]:14.8f} {mmcoords[i][1]:14.8f} {mmcoords[i][2]:14.8f} {mmcharges[i]:9.4f}\n")


def link_coord_corr(
    links,
    qmindex,
    qmkinds,
    qmcoords,
    mmindex,
    mmkinds,
    mmcoords,
    printflag=False,
    method='scale',
    rflat0=1.10,
    scale0=0.73,
):

    qmindex_link = qmindex[:]
    qmkinds_link = qmkinds[:]
    qmcoords_link = qmcoords[:]
    # Can warn user if the scale factor is out of the range 0~1
    if printflag:
        print("qm global index", qmindex)
        # print("mm global index", mmindex)
        print(
            "there are",
            len(links),
            "links to be processed:\n[QM host global index, MM host global index]:",
        )
        [print(link) for link in links]
        print(f"the method for generating link atom coordinate is {method}")

    link_scale = []
    # link_scale will be useful when partitioning forces on link
    # atoms, paritcularly when r(QMhost_link)/r(QMhost_MMhost) is
    # not the same for all links, i.e., costomized scale (to be
    # implemented) for each link or using a flat r(QMhost_link)
    for i in range((len(links))):
        link_coord = [0.000, 0.000, 0.000]
        qmindex_link.append('L' + str(i))
        qm_group_index = qmindex.index(links[i][0])
        mm_group_index = mmindex.index(links[i][1])
        qm_host_coord = np.array(qmcoords[qm_group_index])
        mm_host_coord = np.array(mmcoords[mm_group_index])
        if method.lower() == 'scale':
            # In the 'scale' method, link atom is placed along the QMhost-MMhost
            # bond, and its distance to QMhost is scale0*r_mm_qm
            # We can enable user defined (scaling factor) later, but
            # this feature is secondary compared to user-defined link atom kind
            link_coord = scale0 * mm_host_coord + (1 - scale0) * qm_host_coord
            link_scale.append(scale0)
            r_mm_qm = mm_host_coord - qm_host_coord
            d_mm_qm = np.linalg.norm(r_mm_qm)
        if method.lower() == 'flat':
            # In the 'flat' method, link atom is placed along the QMhost-MMhost
            # bond, and its distance to QMhost is set to be
            # a constant, for all links
            scale = 0
            r_mm_qm = mm_host_coord - qm_host_coord
            d_mm_qm = np.linalg.norm(r_mm_qm)
            scale = rflat0 / d_mm_qm
            link_scale.append(scale)
            link_coord = scale * mm_host_coord + (1 - scale) * qm_host_coord
        link_coord = list(link_coord)
        qmkinds_link.append('H  ')
        # We should enable user defined link atom kind later
        qmcoords_link.append(link_coord)
        if printflag:
            print(
                i + 1,
                "th link is between [QM host global index, MM host global index] =",
                links[i],
            )
            print("link's qm host in-group index", qmindex.index(links[i][0]))
            print("link's mm host in-group index", mmindex.index(links[i][1]))
            print(
                f"qm host atom is {qmkinds[qm_group_index]}, "
                f"coordinate: {qm_host_coord}"
            )
            print(
                f"mm host atom is {mmkinds[mm_group_index]}, "
                f"coordinate: {mm_host_coord}"
            )
            print(f"link atom is H , coordinate: {link_coord}")
            print(
                f"crossing bond length is {d_mm_qm} ang,\nand the scale factor for this link is {link_scale[i]}"
            )

    if printflag:
        prop_print_xzy("qm coords original", qmindex, qmkinds, qmcoords)
        prop_print_xzy("qm coords modified", qmindex_link, qmkinds_link, qmcoords_link)
        # prop_print_xzy("mm coords", mmindex, mmkinds, mmcoords)

    return qmkinds_link, qmcoords_link, qmindex_link, link_scale


def link_force_corr(
    links,
    qmindex,
    qmkinds,
    qmforces,
    mmindex,
    mmkinds,
    mmforces,
    link_scale,
    printflag=False,
):
    if printflag:
        prop_print_xzy("qm forces origin", qmindex, qmkinds, qmforces)
        # prop_print_xzy("mm forces origin", mmindex, mmkinds, mmforces)

    linkcount = 0
    for i in range((len(links))):
        linkcount = linkcount + 1
        linkindex = 'L' + str(i)

        qm_group_index = qmindex.index(links[i][0])
        mm_group_index = mmindex.index(links[i][1])
        link_group_index = qmindex.index(linkindex)
        # link_scale are calculated for each link when generating
        # link atom coordinates, and this dependence of a link atom's
        # coordinate on its hosts' coordinates is used to partition
        # that link atoms' force through chain rule
        link_force = qmforces[link_group_index]
        qm_link_force_partition = link_force * (1 - link_scale[i])
        mm_link_force_partition = link_force * link_scale[i]
        qmforces[qm_group_index] += qm_link_force_partition
        mmforces[mm_group_index] += mm_link_force_partition
        if printflag:
            print(f"the force between the {i+1} th pair")
            print(
                f"[QM host global index, MM host global index] = {links[i]} is\n{link_force}"
            )
            print(
                f"the scale factor is {link_scale[i]} for MM host and {1-link_scale[i]} for QM host"
            )
            print(
                f"QM host parition: {qm_link_force_partition}\nMM host partition: {mm_link_force_partition}"
            )

    for i in range((len(links))):
        linkindex = 'L' + str(i)
        link_group_index = qmindex.index(linkindex)
        qmindex.remove(linkindex)
        qmforces = np.delete(qmforces, link_group_index, 0)

    if printflag:
        prop_print_xzy("qm forces corrected", qmindex, qmkinds, qmforces)
        # prop_print_xzy("mm forces corrected", mmindex, mmkinds, mmforces)

    return qmindex, qmforces, mmforces


def neighborlist_gen(hostcoord, coords, index, bondthreshold=1.7, mode='radius'):

    neighbor_index = []

    if len(coords) != len(index):
        raise Exception("neighbor coords and index do not match")
    if len(coords) == 0:
        raise Exception("there is no neighbor to search for host coordinate")
    if mode.lower()[0:3] == 'rad':
        for coord in coords:
            dist = np.linalg.norm(np.array(coord) - np.array(hostcoord))
            if dist < bondthreshold and dist > 0.1:
                index[coords.index(coord)]
                neighbor_index.append(index[coords.index(coord)])

    if mode.lower()[0:4] == 'near':
        nearest_index = index[0]
        nearest_coord = coord[nearest_index]
        nearest_dist = np.linalg.norm(
            np.array(nearest_coord) - np.array(hostcoord)
        )
        for i in range(len(coords)):
            dist = np.linalg.norm(np.array(coords[i]) - np.array(hostcoord))
            if dist < nearest_dist and dist > 0.1:
                nearest_index = index[i]
                nearest_coord = coords[i]
        neighbor_index = nearest_index

    return neighbor_index


def link_charge_corr(
    mmcharges,
    mmindex,
    links,
    qm_charge,
    system_charge,
    printflag,
    charge_corr_mode='global',
    mmcoords=[],
    mmkinds=[],
    mmneighbor_thrsh=1.7,
):

    mmcharges_redist = copy.deepcopy(mmcharges)

    mmhostindex_global = [link[1] for link in links]
    mmhostindex_group = [mmindex.index(i) for i in mmhostindex_global]
    mmhostcharges = [mmcharges[i] for i in mmhostindex_group]

    if printflag:
        [
            print(
                f"mmhost global index {mmhostindex_global[i]}, in-group index {mmhostindex_group[i]}, and the mmhost charge {mmhostcharges[i]}"
            )
            for i in range(len(links))
        ]

    if charge_corr_mode.lower()[0:3] == 'glo':
        charge_total_qm_classical = system_charge - sum(mmcharges)
        # total_qm_classical_charge + total_mm_classical_charge = system_charge (interger)
        charge_total_mm_host = sum(mmhostcharges)
        # mm hosts classical charges are also spread out to the rest of mm atoms
        charge_zeroed_out = charge_total_mm_host + charge_total_qm_classical
        # residue charge to be spread out over the rest of mm atoms
        charge_corr_total = charge_zeroed_out - qm_charge

        if (len(mmcharges) - len(mmhostcharges)) < 1:
            raise RuntimeError("there is no mm atoms left to spread the charge correction")
        chargecorr = charge_corr_total / (len(mmcharges) - len(mmhostcharges))

        for i in range(len(mmcharges_redist)):
            if i in mmhostindex_group:
                mmcharges_redist[i] = 0.000
            else:
                mmcharges_redist[i] += chargecorr

        if printflag == True:
            print(f"total qm classical charge: {charge_total_qm_classical:10.4f}")
            print(f"+ total mm host charge {charge_total_mm_host:10.4f}")
            print(f"= total zeroed-out classical charge {charge_zeroed_out:10.4f}")
            print(
                f"- (qm charge : {qm_charge}) = total charge to spread over the remaining mm atoms {charge_corr_total:10.4f}"
            )
            print("remaining mm atom number", (len(mmcharges) - len(mmhostcharges)))
            print(f"each remaining mm atom adds {chargecorr:10.4f}")

    elif charge_corr_mode.lower()[0:3] == 'loc':
        print("mmhost charges will be distributed to its neighbors")
        if len(mmcoords) != len(mmcharges):
            raise Exception("mmcoords length does not match mmcharges, did you forget to input or input the wrong mmcoords for mmhost local charge distribution?")
        if len(mmkinds) != len(mmcharges):
            raise Exception("mmkinds length does not match mmcharges, did you forget to input or input the wrong mmkinds?")
        for link in links:
            mmhost_groupindex = mmindex.index(link[1])
            mmhost_coord = mmcoords[mmhost_groupindex]
            hostcharge = mmcharges_redist[mmhost_groupindex]
            if abs(hostcharge) > 0.2:
                warnings.warn("mmhost charge is higher than 0.2, is the link's choice appropriate?")
            mmcharges_redist[mmhost_groupindex] = 0
            mmneighbor_index = []
            mmneighbor_dists = []
            mmneighbor_kinds = []
            for i in range(len(mmcoords)):
                dist = np.linalg.norm(np.array(mmcoords[i]) - np.array(mmhost_coord))
                if dist < mmneighbor_thrsh and dist > 0.1:
                    mmneighbor_index.append(mmindex[i])
                    mmneighbor_dists.append(dist)
                    mmneighbor_kinds.append(mmkinds[i])
            mmneighbor_charges=[]
            mmneighbor_charges0=[]
            mmneighbor_num=len(mmneighbor_index)
            for i in mmneighbor_index:
                j = mmindex.index(i)
                mmcharges_redist[j] += hostcharge/mmneighbor_num
                mmneighbor_charges.append(mmcharges_redist[j])
                mmneighbor_charges0.append(mmcharges[j])
            if printflag == True:
                print(f"{link=}")
                print(f"{mmhost_coord=}")
                print(f"mmhost is a {mmkinds[mmhost_groupindex]} containing {mmcharges[mmhost_groupindex]}, current charge {mmcharges_redist[mmhost_groupindex]}")
                print(f'{mmneighbor_index=}')
                print(f'{mmneighbor_dists=}')
                print(f'{mmneighbor_kinds=}')
                print(f'{mmneighbor_charges=}')
                print(f'{mmneighbor_charges0=}')
    else:
        if charge_corr_mode.lower()[0:3] != 'del':
            warnings.warn("charge correction method cannot be parsed, delete all mmhost charges.")
        else:
            print("mmhost charges will be deleted.")
        for i in mmhostindex_group:
            mmcharges_redist[i] = 0
        [print(f"{mmcharges[i]=}") for i in mmhostindex_group]
        [print(f"{mmcharges_redist[i]=}") for i in mmhostindex_group]

    return mmcharges_redist

def switch(r, rs, rc):
    if r < rs and r > 0:
        return 1
    elif r >= rs and r <= rc:
        x = (r-rs)/(rc-rs)
        return (x**3 * (-6*x**2 + 15*x - 10) + 1)
    else:
        return 0
