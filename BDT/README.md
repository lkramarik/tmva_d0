# Boosted Decision Trees for D<sup>0</sup> reconstruction

This macros are made for training and application of BDT from ROOT's TMVA package.

## Running of BDT

In principal, everything can be run with `runTmva.sh`. Just set the p<sub>T</sub> interval, in which you want to train and apply BDT. It's best to comment/uncomment lines here, depending on what you want to do.

### BDT training

Macro: `TMVAClassification.C`

Input is a list of files containg background pairs.

Training results are, as always, in the folder `weights` and can be visualised with `TMVAGui.C`. For visualisation of optimal BDT cut,
expected number of signal and backgroung need to be known. This is bit hard to estimate, thus data-driven estimation is used 
(see the next sections).

FastSim output is used as signal in the training. 
Background is taken from data, currently its statistics seems to be limited and needs to be enriched in the future. Background is from sideband or wrong sign combination - need to be set.

Remember to set reasonable number of training background and signal pairs. The larger number is better, however we are limited by the statistics in the given transverse momentum range.



### BDT Application

Macro: `TMVAClassificationApplication.C`

Input is a list of files containg signal and background pairs, both from data (you are applying your training to data).
Output is a file with corresponding signal and backgroun NTuples, with BDT value for each pair.

### Data-driven significance estimation
