.. role::  raw-html(raw)
    :format: html

*PROBE* : PeRiod and OrBital Evolution Code 
==========================================

Designed to follow to evolution of both stellar rotation rate and planetary orbital period. Developed for low mass stars 0.5 < M\ :sub:`star`\/M\ :sub:`sun`\  < 1.0.

The choice of initial SMA and planetary mass is done in condinit.f90 (sma_tabl and mass_planet)

.. code-block:: bash

    perl Makefile.py
    ./condinit
    
The description of the input parameters for the rotation part can be found here https://github.com/GalletFlorian/JEVOL

    
The BGM version of PROBE include, in the rotational evolution part, the evolution of the mass of the star.
Since the mass of the star only decreases a little between the PMS and the end of the MS, this evolution is usually not included. 