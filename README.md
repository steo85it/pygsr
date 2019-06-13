# pygsr
[![DOI](https://zenodo.org/badge/191787378.svg)](https://zenodo.org/badge/latestdoi/191787378)


This code implements a simplified version of a Gaia-like Global Sphere Reconstruction code in Python3. Currently, it uses simulated observations provided by Alexey Butkevich, of the Pulkovo Observatory, Russia, which considers 1000 stars with an approximately uniform distribution on the celestial sphere. Observations are then generated for 5 years according to a nominal Gaia-like scanning law.

The general goal is to give a "hands-on" overview of the main numerical and statistical properties of the global astrometric problem.

The astrometric model is simplified with respect to:

   - The GR metric, which takes into account just the Schwarzschild-like contribution of the Sun (no planetary contributions).
   - The observation equation, which considers only the case of 2 astrometric unknowns (i.e. positions); parallaxes, proper motions, attitude, calibration and global parameters are the true ones and are not reconstructed.

On the other side, the attitude is modelled in a fully consistent relativistic way.

The code currently gets a set of observation epochs and auxiliary data (star catalog, Gaia ephemeris and attitude, Sun ephemeris), simulates observations by applying both catalog errors and observation noise and solves for a set of astrometric parameters (positions, parallaxe and proper motions).

The main classes implemented in the code are: a "star" class and an "obs_eq" class, including all methods required to set-up and solve the least-square problem from data related to each star. No global parameters (e.g., attitude corrections, ppn-parameters) are included.

Examples of how to use the code and the provided classes and methods are provided by the google colab notebook in the example/ dir as list of exercises from the Gaia Summer School 2019.
The full exercise is available at https://colab.research.google.com/drive/1chbPwxgDJmOIQ87x3Z-v3447ioSIZNSg (Google account required).
