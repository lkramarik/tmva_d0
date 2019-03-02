# TMVA Cuts method for D<sup>0</sup> meson analysis

This set of scripts works with ROOT v6-14-02. It will not work with ROOT 5 and probably not work with newer versions.

Training and application is done for different p<sub>T</sub> bins and passes in it. Pass is used in order to iterate cuts, setting 
more strict pre-cuts on the training sample. These pre-cuts are loaded from `tmvaCuts.h` (currently passes are not well worked-out).

The script that runs all is `./runTmva.sh` - create folders, copy needed files and finally run training and application.
Significance plot from training, with number of signal and background events, is calculated (TMVAGui is not needed for this)
and saved in ROOT file with training results (currently `tmvaD0.root`).


In `TMVAClassificationApplication.C`, invariant mass peak histograms could be made for cuts trained for all signal efficiencies,
that produces significance vs signal efficiency plot. This one could be compared with the one from training. In addition,
all invariant mass peaks are saved in ROOT files. Currently, it may happen, that changig training method (for example
from "Cuts" to "CutsGA") will require to change few lines in `TMVAClassificationApplication.C`.

Data file:
```sh
/gpfs01/star/pwg/lkramarik/Dmaker_ndAu/Dmaker_dAu/workDir/2019-02-08_04-15_D0/production/ntp.picoD0AnaMaker.0802.0415.root
```
Simulation (FastSim) file:
```sh
--will be created--
```
