# HistPlotter
 - Header-only class to encapsulate CERN ROOT 1D/2D Histogram plotting and application of TCuts.
 - Can specify folder hierarchy while setting up fills, currently supports only one level.
 - Cuts specified using a two-column text file containing cut names, and target .root files. These .root files must contain a TCutG of name "CUTG".
 - Tested for use in macros, with TSelector design pattern and compiled code.
 - Idea inspired from MyFill() pattern created by github user gwm17
