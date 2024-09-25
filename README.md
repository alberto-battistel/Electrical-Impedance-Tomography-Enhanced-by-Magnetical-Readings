# Electrical Impedance Tomography Enhanced by Magnetical Readings

This repo contains several matlab scripts to investigate the addition of magnatic reading on Electrical Impedance Tomography (EIT).  
Every version corresponds to a different published work.

## v0.0.0 BMS2024 Poster

v.0.0.0 contains the data and produces the figures used for the __[BMS2024/AUTOMED2024 conference in Villingen-Schwenningen](https://www.bms-24.org/)__ Poster Abstract Proceedings titled:

__Enhance electrical impedance tomography with magnetic readings__

by:
A. Battistel 1*, A.C. Özen 2, J. Wilkie 1, R. Chen 1,3, D. von Elverfeldt 2, and K. Möller 1,3

1. Institute of Technical Medicine (ITeM), Furtwangen University, Jakob-Kienzle-Strasse 17, 78054 Villingen-
Schwenningen, Germany
2. Division of Medical Physics, Department of Diagnostic and Interventional Radiology, University Medical Center Freiburg, Faculty of Medicine, University of Freiburg, Freiburg, Germany
3. Faculty of Engineering, University of Freiburg, Freiburg, Germany
* Corresponding author, email: alberto.battistel@hs-furtwangen.de

and it won a poster prize award.

### What you need to run this code

- Matlab (it was tested with Matlab R2023b and R2024b on a linux machine)
- [EIDORS 3.11 with Netgen 5.3](https://eidors3d.sourceforge.net/)
- A matlab function called `init_eidors.mat` available to your path containing
`run /path/to/eidors3d/startup.m` (see [EIDORS first tutorial](https://eidors3d.sourceforge.net/tutorial/EIDORS_basics/one_line.shtml))

### How to use this repository

By running `sphere_test.m` you will get figure 1, and 2, of the conference paper plus additional plots that were used for the poster.
Please note that the figures you will get need formatting to appear exactly as those of the conference paper.
There are also a lot of other files around, just entropy...

## Cite this

