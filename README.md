# Hydrogenic.vX

This repository contains the hydrogenic atom module used in the CMBSPEC
atomic-model setup. The code was developed for detailed cosmological
recombination calculations, where high-lying hydrogenic levels, dipole and
quadrupole transitions, photoionization, and recombination rates have to be
treated efficiently and consistently.

The repository is intended to hold the core hydrogenic atom classes and the
physics routines that belong directly to them, without copying the more general
Tools or Development libraries from the larger codes.

The main interface is provided by `Atom.h` / `Atom.cpp`. The module sets up
hydrogenic atoms and ions with nuclear charge `Z`, resolved angular-momentum
sublevels, bound-bound transition data, optional electric quadrupole lines,
photoionization cross sections, and recombination-related rates.

## Data

The hydrogenic module is mostly analytic and recursive, so it does not carry a
large module-specific runtime data directory. Tables and cached interpolation
data, when used by an application, are generated through the surrounding build
and runtime setup.

## Dependencies

This repository is not a standalone application. It depends on common helper
code from the surrounding CMBSPEC/CosmoSpec/CosmoTherm toolchains, including
physical constants, numerical helper routines, integration routines, file I/O,
and Voigt-profile support. Some recombination and interpolation paths also use
external numerical libraries such as GSL.

A small external demonstration project,
[Hydrogenic-demo](https://github.com/CMBSPEC/Hydrogenic-demo), shows the
minimal set of required Tools, a concrete C++ build example, and Python/Jupyter
interfaces for interactive use.

## Related Literature And Data Sources

The code is connected to several atomic-physics and cosmological recombination
references, including:

- J. Chluba and R. M. Thomas, "Towards a complete treatment of the
  cosmological recombination problem", MNRAS 412, 748, 2011.
- J. Chluba, J. A. Rubino-Martin, and R. A. Sunyaev, "Cosmological hydrogen
  recombination: populations of the high level sub-states", MNRAS 374, 1310,
  2007.
- J. Chluba and Y. Ali-Haimoud, "CosmoSpec: Fast and detailed computation of
  the cosmological recombination radiation from hydrogen and helium",
  MNRAS 456, 3494, 2016.
- P. J. Storey and D. G. Hummer, "Fast computer evaluation of radiative
  properties of hydrogenic systems", Computer Physics Communications 66, 129,
  1991.
- W. J. Karzas and R. Latter, "Electron Radiative Transitions in a Coulomb
  Field", ApJS 6, 167, 1961.
- J. D. Hey, "On the determination of radial matrix elements for high-n
  transitions in hydrogenic atoms and ions", J. Phys. B 39, 2641, 2006.
- D. Grin and C. M. Hirata, "Cosmological hydrogen recombination: The effect of
  extremely high-n states", Phys. Rev. D 81, 083005, 2010.

These references are listed to document the origin and scientific context of
the methods used here; the repository itself only provides the hydrogenic atom
setup.

These repositories were made available and documented with the help of Codex. The related release work was supported in part by a grant of access to OpenAI models through the ChatGPT for Academic Researchers program.
