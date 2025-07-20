This is the repository of codes implemented in recent manuscripts by Wang Lab in the Department of Physics at the University of Tulsa:

1. Remeshing flexible membranes under the control of free energy, PLOS Computational Biology, 2022, 18(12): e1010766, written by Xinxin Wang and Gaudenz Danuser, date: 2022/07/25.
Contact email: xinxin-wang@utulsa.edu

The membrane is modularized as a matlab object @ModMembrane, and external control points as @ModSubstrate. @ModMembrane follows the physics-based remeshing algorithm in the MS to be able to represent generic morphologies, e.g. red blood cell. @ModSubstrate is mechanically coupled to @ModMembrane via @TypForce for simulating the morphologies of various case studies in the MS. All the modular objects are organized by the object @model before running the simulations. All the available simulations are written in @BiophysicsApp, and can be run from 'Scripts/runMS2022Examples.m'. Graphics are adjusted to the exact formats in the MS, and can be run from 'Scripts/runMS2022Graphics.m'.

