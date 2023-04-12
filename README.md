# DockSurf
Docking protein onto inorganic surface

DockSurf: a molecular modelling tool for protein/surface prediction

Theory overview:
One key aspect in the development of bionanotechnologies is the elucidation of the structural interfaces between proteins and inorganic surfaces[1]. Despite its wide importance, the interfacial structures between proteins and metallic surfaces remain elusive and the lack of experimental investigation refrains the development of many devices[2].In this context, molecular simulations represent an affordable alternative to understand the underlying physical processes which occur at the protein/surface interface. According to the size of proteins, molecular dynamics simulations are the reference technique but required an initial structural organization for studying the interaction with a surface. To overcome this drawback, we propose to consider the generation of protein/surface structure as a molecular docking problem for which the target is an homogenous plan.

A molecular docking software can be viewed as a conformation exploration investigation and an interaction energy computation. Our approach and originality is to consider the conformational exploration with Euler's angles, like a plane above a surface, providing this way a cartography instead of an unique structure. To date, only Au{111} surfaces are parametrized. Interaction energies were derived from either QM or MM computations for a set of small molecules that describe protein atom-types, and implemented in a DLVO potential[3]. Validation of DockSurf software has been made with MD for several corona proteins kwnown to interact with Au{111} surface.

References:
1 Ozboyaci, M.; Kokh, D. B.; Corni, S.; Wade, R. C. Quarterly Reviews of Biophysics 2016, 49, e4.
2 Dalsin, J. L.; Messersmith, P. B. Materials Today 2005, 8 (9), 38â€“46.
3 Overbeek, J. T. G. Journal of Colloid and Interface Science 1977, 58 (2), 15.

Usage:
Only few information are required for DockSurf:
- A label name which is used as a nickname for the generation of all output names.
- A pdb (protein data bank) file which contains the protein coordinates.
- An integer, called fine, which define the angle increment (typical value is 5 or 10).
- Surface type label. To date, only Au111 is available, further improvments will be made in the future.
- An integer, sel, which defines the number of output conformations to be selected (typically 3 or 5).
- Clus or conf tag to choose either the sel firt clusters or first conformations.
