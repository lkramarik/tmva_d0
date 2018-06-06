# Boosted Decision Trees for D<sup>0</sup> reconstruction

This macros are made for training and application of BDT from ROOT's TMVA package.

### BDT training

Macro: `TMVAClassification.C`

Training results are, as always, in the folder `weights` and can be visualised with `TMVAGui.C`. For visualisation of optimal BDT cut,
expected number of signal and backgroung need to be known. This is bit hard to estimate, thus data-driven estimation is used 
(see the next sections).

FastSim output is used as signal in the training. 
Background is taken from data, currently its statistics seems to be limited and needs to be enriched in the future.

Remember to set reasonable number of training background and signal pairs. The larger number is better, however we are limited by
the statistics in the given transverse momentum range.



### BDT Application

Macro: `TMVAClassificationApplication.C`
