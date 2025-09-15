.. ExoComp documentation master file, created by
   sphinx-quickstart on Fri Jun 13 14:00:36 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
Welcome to ExoComp!
===================================

.. image:: _images/Exocomp_Banner_Dark_BG.png
   :alt: Project logo
   :width: 700px
   :align: center
	   
`Many different retrieval codes exist <https://ui.adsabs.harvard.edu/abs/2023RNAAS...7...54M/abstract>`__, each with different methods for quantifying atmospheric compositions. In general, chemical equilibrium retrievals quote a metallicity and C/O, but the way in which these are defined varies. For some codes, metallicity is parameterized by O/H, while C/O adjusts the C/H ratio. But in other codes, it is the opposite!

Exocomp is a toolkit with functions and utilities for normalizing measurements from different retrievals of exoplanets and brown dwarfs. Features include:

* Conversions between volume mixing ratio and mass fraction (and vice versa).
* Conversions in bulk abundances between different definitions of solar composition (or between definitions of host star abundances).
* Conversions from metallicity, C/O, and refractory-to-volatile ratios to elemental abundance ratios.
* Fitting of metallicity, C/O, and refractory-to-volatile ratios to sets of atomic and molecular abundances.
* Calculation of equilibrium molecular abundances for desired species given metallicity, C/O ratio, temperature, and pressure.

All of the above can be easily done with different solar (or stellar) abundance definitions. Available solar abundance defintions include:

* Lodders+ 2025
* Lodders+ 2020
* Lodders+ 2010
* Asplund+ 2021
* Asplund+ 2009
* Asplund+ 2005
* Caffau+ 2011

`Link to the code <https://github.com/jlothringer/exocomp>`__.
  
We also have compiled the `ExoComp Table <https://jlothringer.github.io/exoplanet_datatable.html>`__, which lists metallicity, C/O, and refractory-to-volatile measurements from retrievals of exoplanets from direct, eclipse, and transmission spectroscopy.

.. toctree::
   :maxdepth: 1
   :caption: Installation:

   installation

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules
   ExoComp Tutorial <Tutorial1>
   Population Comparison <population_comparison>
   The ExoComp Table <https://jlothringer.github.io/exoplanet_datatable.html>
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
