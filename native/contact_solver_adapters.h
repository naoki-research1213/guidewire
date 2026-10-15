#pragma once
#include "constraint_adapters.h"
#include "contact_model_adapters.h"
#include <sofa/helper/LCPcalc.h>
#define SCOPED_TIMER_VARNAME(var,name) (void)0
#define msg_info() Quiet()
namespace sofa::helper {
 struct ScopedAdvancedTimer{explicit ScopedAdvancedTimer(const char*){}};
 namespace advancedtimer {inline void stepBegin(const char*){}inline void stepEnd(const char*){}}
}
struct ContactSolver {
 using PersistentID=long long;
 struct ConstraintBlockBuf {int nbLines=0;std::map<PersistentID,int>persistentToConstraintIdMap;};
 std::map<core::behavior::BaseLagrangianConstraint*,ConstraintBlockBuf>_previousConstraints;
 VecConstraintBlockInfo m_constraintBlockInfo;VecPersistentID m_constraintIds;std::vector<double>_previousForces;std::vector<double>*_result=nullptr;
 unsigned _numConstraints=0;
 Data<double>d_mu{parameters().mu},d_tol{parameters().tolerance};Data<int>d_maxIt{static_cast<int>(parameters().maxIterations)};
 std::vector<Correction<RigidData>*>l_constraintCorrections,_cclist_elem1,_cclist_elem2;
 ConstraintMatrix _Wdiag;ConstraintMatrix*_W=&_Wdiag;
 std::vector<double>unbuilt_d;std::vector<sofa::helper::LocalBlock33>unbuilt_W33;
 void addRegularization(ConstraintMatrix&){}
 int nlcp_gaussseidel_unbuilt(SReal*,SReal*,std::vector<SReal>*);
 void computeInitialGuess();void keepContactForcesValue();
};
