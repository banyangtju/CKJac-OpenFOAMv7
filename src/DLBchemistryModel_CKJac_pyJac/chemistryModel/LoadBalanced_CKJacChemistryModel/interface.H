#ifndef __INTERFACE_H__
#define __INTERFACE_H__

#include <string>

#ifndef SUBROUTINE
#ifdef __linux
#define SUBROUTINE	void
#else
#define SUBROUTINE	void __stdcall 
#endif
#endif

#define CALL_OMEGA_OF_C omega_of_c_
#define CALL_JAC_OF_C jac_of_c_

SUBROUTINE CALL_OMEGA_OF_C(const double* T, double* C, const double* p, const int* KK, double* WDOT, double* DTDt);

SUBROUTINE CALL_JAC_OF_C(const int* KK, const double* p, const double* T, double* C, double* fortJAC);


#endif //__INTERFACE_H__

