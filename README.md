# pygsr

This code implements a simplified version of a Gaia-like Global Sphere Reconstruction code in Python3. Currently, it uses simulated observations provided by Alexey Butkevich, of the Pulkovo Observatory, Russia, which considers 1000 stars with an approximately uniform distribution on the celestial sphere. Observations are then generated for 5 years according to a nominal Gaia-like scanning law.

The general goal is to give a "hands-on" overview of the main numerical and statistical properties of the global astrometric problem.

The astrometric model is simplified with respect to:

   - The GR metric, which takes into account just the Schwarzschild-like contribution of the Sun (no planetary contributions).
   - The observation equation, which considers only the case of 2 astrometric unknowns (i.e. positions); parallaxes, proper motions, attitude, calibration and global parameters are the true ones and are not reconstructed.

On the other side, the attitude is modelled in a fully consistent relativistic way.

Examples of how to use the code and the provided classes and methods are provided by the google colab notebook in the example/ dir as list of exercises from the Gaia Summer School 2019.
