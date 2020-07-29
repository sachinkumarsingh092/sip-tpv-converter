# Sip-tpv converter

Sip to tpv and tpv to sip distortion convertion.
This mainly follows the algorithm as discussed in the paper:

[“More flexibility in representing geometric distortion in astronomical images,” Shupe, David L.; Laher, Russ R.; Storrie-Lombardi, Lisa; Surace, Jason; Grillmair, Carl; Levitan, David; Sesar, Branimir, 2012, in Software and Cyberinfrastructure for Astronomy II. Proceedings of the SPIE, Volume 8451, article id. 84511M.](http://web.ipac.caltech.edu/staff/shupe/reprints/SIP_to_PV_SPIE2012.pdf)

This repository contains the test files before merging them to GnuAstro.

**/scripts** contains python scripts to be used for some automation and calculations.

**/test-fits** contains fits files used to test the functions of the library.

`temp-wcs.c` is the main file containing the implementation of the algorithm.