#include "implementation.h"
int main(){
 Kernel kernel;for(int i=0;i<100;i++)kernel.advanceWithoutContact(1.,0.);
 kernel.freeMotion(0.,0.);
 ConstraintState state{&kernel.state};state.freePosition.value=kernel.freePositions;state.freeVelocity.value=kernel.freeVelocities;
 for(int i=0;i<3;i++){Vec6 v;v[i]=1.;state.jacobian.value.rows[i][60]=v;}
 ConstraintLinearSystem system{&kernel.system,&state};CorrectionSolver linear{&kernel.solver,&system};OdeFactors ode;
 Correction<RigidData>correction{&state,&linear,{&ode}};
 ContactSolver contact;contact._numConstraints=3;contact._Wdiag.resize(3,3);contact.l_constraintCorrections={&correction};
 std::vector<double>d={-.2,.1,0.},force(3,0.),residual;
 int solved=contact.nlcp_gaussseidel_unbuilt(d.data(),force.data(),&residual);
 std::cout<<"converged="<<solved<<" iterations="<<residual.size()<<" force="<<force[0]<<","<<force[1]<<","<<force[2]<<" residual="<<residual.back()<<"\n";
 if(!solved||force[0]<0||std::hypot(force[1],force[2])>.15*force[0]+1e-12)return 2;
 Vector lambda(1);for(int i=0;i<3;i++)lambda.set(i,force[i]);
 correction.applyContactForce(&lambda);
 return 0;
}
