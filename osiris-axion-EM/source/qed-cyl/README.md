# qed-cyl Algorithm

## Overview

The qed-cyl algorithm is a new simulation mode enabling to use at the same time the cylindrical mode decomposition and to account for the QED processes. Those processes are for now non-linear inverse Compton scattering and Breit-Wheeler pair production.


## Input parameters

The input parameters are the same as for a quasi-3D deck, and one then has to add the parameters from a QED deck. The main parameters of a quasi-3D input deck and all of the parameters from the QED input deck are recalled below

* The qed-cyl algorithm is currently implemented as a simulation mode. To activate it the following parameter must be added to the ``simulation`` section
  
  * `algorithm = "qed-cyl",`

* The following parameter has been added to the ``grid`` section

  * ``n_cyl_modes``, integer, default = 0

    This integer specifies how many azimuthal modes to use, ranging between 0 and 15.  Using 0 modes will assume complete azimuthal symmetry and run with only the real part of mode 0 (no imaginary part).

* As before, the number of QED groups is controlled in the ``particles`` section, using the ``num_qed`` option. A QED group consists in 3 species : an electron specie, a positron specie and a photon specie (in this exact order).

```
  particles 
  { 
    num_species = 0, 
    num_qed = 1,
  }
````

* For each QED group, a ``qed_group`` section can be included to control the parameters for photon emission and pair production within this group:
  
```
  qed_group
  {
    if_damp_classical = .true.,
    if_damp_qed = .true.,
    if_pairprod = .true.,
    
    qed_g_cutoff = 10.0,
    p_emit_cutoff = 2.0,
  }
```


  In this section, the relevant parameters are:

  * ``if_damp_classical`` and ``if_damp_qed`` - Control whether electrons and positrons lose energy due to classical radiation reaction / QED photon emission
  * ``if_pairprod`` - Controls whether photons emit electron-positron pairs
  * ``qed_g_cutoff`` - Sets the threshold energy between classical radiation reaction and QED photon emission
  * ``p_emit_cutoff`` - Sets a minimum energy for which photons emitted due to QED processes will be tracked

Following each ``qed_group`` section, a section ``qed_group_diag`` can be included to control the parameters for diagnostics relative to the QED group:

```
qed_group_diag
{
  ndump_fac_rad = 1,
  ndump_fac_pairs = 1,

  ndump_fac_radspect = 10,
  radspect_emin = 1.957d-5,
  radspect_emax = 1.957d5,
  radspect_bins = 400,

  ndump_fac_chi_emit = 10,
  chi_emit_min = 1.d-4,
  chi_emit_max = 1.d2,
  chi_emit_nbins = 240,

  ndump_fac_raddetector = 1,
  raddetector_nbins = 100
  raddetector_emin = 0
  raddetector_emax = 1.d15
  raddetector_phimin = -0.05,
  raddetector_phimax = 0.05,
  raddetector_thetamin = -0.05, 
  raddetector_thetamax = 0.05,
  detector_direction = "x1_positive",
  
  ndump_fac_radsphere = 1,
  radsphere_nbins = 100, 
  radsphere_emin = 0.0
  radsphere_emax = 1.d15
  if_add_classical_radsphere = .false.

  ndump_fac_n_emit = 1,
  ndump_fac_chi = 1,
}
```

In this section, the relevant parameters are:

  * ``ndump_fac_rad`` - Controls the dump frequency of the global radiated energy diagnostic, which includes independent columns for energy radiated by charged particles in the classical and QED regimes (for tracked and non-tracked photons); this diagnostic produces one file for each electron and positron species in the QED group (see ``HIST/par01_rad`` and ``HIST/par02_rad``)
  * ``ndump_fac_pairs`` - Controls the dump frequency of a global diagnostic that reports the total number, energy and weight of pairs emitted by the photon species in the QED group (see HIST/pho01_ene and HIST/pho01_pairs). The diagnostic is stored in ``HIST/pho01_ene/`` and ``HIST/pho01_pairs/``.
  * ``ndump_fac_radspect``, ``radspect_emin``, ``radspect_emax`` and ``radspect_bins`` - Respectively control the dump frequency, minimum and maximum energies, and number of bins of a runtime diagnostic of the energy spectrum of the photons emitted by each electron and positron species; the binning is only done logarithmically (see QED/RADSPECT/). This diagnostic is normalised in J/MeV for a 3D simulation, and in J/MeV/m for a 2D simulation. The min/max bounds and number of bins are optional and have the default values indicated in the example. The diagnostic is stored in ``QED/radspect/``.
  * ``ndump_fac_chi_emit``, ``chi_emit_min``, ``chi_emit_max`` and ``chi_emit_nbins`` - Respectively control the dump frequency, minimum and maximum chi parameters, and number of bins of a runtime diagnostic of the quantum parameter spectrum of both the electrons/positrons emitting and the photons emitted; the binning is only done logarithmically (see QED/chi_emission/). The min/max bounds and number of bins are optional and have the default values indicated in the example. The diagnostic is stored in ``QED/chi_emission/``.
  * ``ndump_fac_raddetector`` - Controls the frequency of a "detector-like" 2D diagnostic for photons emitted via NLIC, stored in ``QED/rad-detector/``. All optional parameters are given default values showed in the example. These optional parameters are the number of bins ``raddetector_nbins``, the min/max photon energies selected ``raddetector_emin``, ``raddetector_emax``, the min/max azimuthal angles selected ``raddetector_phimin``, ``raddetector_phimax``, the min/max polar angles selected ``raddetector_thetamin``, ``raddetector_thetamax``, and the direction of the detector ``detector_direction``.
  * ``ndump_fac_radsphere`` - Controls the frequency of a photon energy diagnostic on a spherical surface, stored in ``QED/radsphere/``. It is done via six 2D diagnostics. The six grids correspond to the six directions +/-x1, +/-x2 and +/-x3. They enable to reconstruct a spherical 4 pi steraian detector . Optional parameter are given default values whichcan be tuned. ``radsphere_nbins`` is the number of bins, ``radsphere_emin``, ``radsphere_emax`` are the min and max energies of photons selected and ``if_add_classical_radsphere`` is a boolean to add the contribution of photons radiated classically. 
  * ``ndump_fac_n_emit`` - Controls the dump frequency of a grid diagnostic with the number of photon emission/pair production events per species in each PIC cell.  The diagnostic is stored in ``QED/n_emit/``.
  * ``ndump_fac_chi`` - Controls the dump frequency of a grid diagnostic with the weighted-average quantum parameter per species in each PIC cell. To get this diagnostic for photons, one needs to turn on the pair production flag. The diagnostic is stored in ``QED/chi_grid/``.

After these sections, a species section (and corresponding mandatory and/or optional udist, profile, spe_bound and diag_species sections) should follow for the electron, positron and photon species in the QED group (in this order). For the moment, only the standard (Boris) particle pusher is allowed for charged particles, so those species' sections should include push_type="standard".


* The following parameters have been added to the ``species`` section

  * ``num_par_theta``, integer, default = –1

    This parameter specifies how many particles to space in rings in *ϕ*.  We recommend this number be at least 8\*``n_cyl_modes``.

  * ``rand_theta``, logical, default = .false.

    This flag specifies if a random offset should be added in theta to particle rings in each cell (one random number generated per cell if .true.).  This can help sample modes more efficiently, but may also add noise to the simulation.