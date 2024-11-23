# OpenFOAMv7-CKJac-LoadBalancedChemistryModel
The LoadBalancedChemistryModel is derived from DLBFoam (https://github.com/Aalto-CFD/DLBFoam), and integrates the recently developed CKJac (for more details in https://www.sciencedirect.com/science/article/abs/pii/S0010218024004978) chemistry model for computational speedup.

The chemical model program, CKJac, provides the pre-computed and pre-saved fully analytic molar concentration-based and more spare Jacobian and chemical reaction sources, which accelerate chemical reaction computation, which occupies 80% CPU time during CFD simulations. The speedup performance depends on chemical mechanisms.

When using the CKJac chemistry model, a library called libchemkin.so should be generated in constant/chemkin. The libchemkin.so is produced based on chemkin-type reaction mechanisms. You can give me the mechanism and I will provide the library for you.
