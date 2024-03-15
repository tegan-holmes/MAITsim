# MAITsim
Repository for the Model Agnostic Isotope Tracer simulator (MAITsim)

MAITsim models isotope tracers for the outputs of hydrologic models. The basic verison of MAITsim is a set of functions in either R or MATLAB which require hydrologic storage state and flux data as inputs and outputs isotope tracer compositions for all input storages. MAITsim can link to any hydrologic model structure, provided the model tracks stored water, moves water between storage units, has only precipitation as a water input from outside the modeled area, and has a stable uni-directional flow path.

Information on the exact data needed for inputs and the formatting requirements can be found in the documentation files.

There are seperate MAITsim functions for oxygen-18 and hydrogen-2 modeling.
