---
title: 'eSCAPE: parallel global-scale landscape evolution model'
tags:
  - Python
  - landscape evolution
  - geomorphology
  - hydrology
  - surface processes
  - stratigraphy
authors:
 - name: Tristan Salles
   orcid: 0000-0001-6095-7689
   affiliation: "1"
affiliations:
 - name: School of Geosciences, The University of Sydney, Australia
   index: 1
date: 16 September 2018
bibliography: paper.bib
---

# Summary

**eSCAPE** is a parallel landscape evolution model, built to simulate Earth surface dynamics at global scale and over geological times. The model is primarily designed to address problems related to geomorphology, hydrology, and stratigraphy, but it can also be used in related fields.

**eSCAPE** accounts for both hillslope processes (_soil creep using linear diffusion_) and fluvial incision (_stream power law_). It can be forced using spatially and temporally varying tectonics (vertical displacements) and climatic conditions (precipitation changes and/or sea-level fluctuations).

The model computes flow accumulation using multiple flow direction over unstructured grids based on an adaptation of the implicit approach proposed by Richardson et al. [@Richardson:2014]. An extension of the parallel priority-flood depression-filling algorithm from [@Barnes:2016] to unstructured mesh is used to simulate sedimentation in upland areas and internally drained basins. Marine sedimentation is based on a diffusion algorithm similar to the technique proposed in [pybadlands](https://github.com/badlands-model/pyBadlands_serial) [@Salles:2018].

# References
