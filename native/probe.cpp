#include "adapters.h"
#include "generated/beam_functions.inc"
#include "solver_adapters.h"
#include "generated/btd_functions.inc"
#include "spring_adapters.h"
#include "generated/spring_functions.inc"
int main(){Interpolation interp; ForceField force{&interp}; BeamLocalMatrices b{};b._L=2;b._A=3;b._Iy=b._Iz=1;b._J=2;b._youngM=20000;b._cPoisson=.3;b._rho=1.55e-6;force.computeStiffness(0,b);force.computeMass(0,b);std::cout<<b.m_K00[0][0]<<" "<<b.m_M00[0][0]<<"\n";}
