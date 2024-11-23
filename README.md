# OpenFOAMv7-CKJac-LoadBalancedChemistryModel
The LoadBalancedChemistryModel is derived from DLBFoam (https://github.com/Aalto-CFD/DLBFoam), and integrates the recently developed CKJac (for more details in https://www.sciencedirect.com/science/article/abs/pii/S0010218024004978) chemistry model for computational speedup.

The chemical model program, CKJac, provides the pre-computed and pre-saved fully analytic molar concentration-based and more spare Jacobian and chemical reaction sources, accelerating chemical reaction computation, which occupies 80% CPU time during CFD simulations. The speedup performance depends on chemical mechanisms.

When using the CKJac chemistry model, a library called libchemkin.so should be generated in constant/chemkin. The libchemkin.so is produced based on chemkin-type reaction mechanisms. You can give me the mechanism and I will provide the library for you.

Fig. 1 The sparsity of CKJac
![image](https://github.com/user-attachments/assets/3f202ddb-3050-4888-8d06-fc7e96d06ec2)


Fig. 2 The computational CPU time comparison.
![image](https://github.com/user-attachments/assets/0d170cf7-4d10-4ded-9267-55c165c5cbcf)

