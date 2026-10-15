#pragma once
#include "solver_adapters.h"
#include "spring_adapters.h"
struct RigidDerivativeView {
 const Vec6& v;
 Vec3 getVCenter()const{return Vec3(v[0],v[1],v[2]);}
 Vec3 getVOrientation()const{return Vec3(v[3],v[4],v[5]);}
};
#include "generated/rigid_integrate.inc"
inline Pose operator+(const Pose& p,const Vec6& dx){Pose out=p;integratePose(out,RigidDerivativeView{dx});return out;}
// 現行設定の接続。計算関数の本体は generated/ にある原文を使用する。
// 接触検出・摩擦拘束の求解・状態補正は scene.h が接続する。
struct Kernel {
 const double h=parameters().dt;
 MechanicalState state;FixedConstraint fixed;Interpolation interpolation;
 Controller controller{&state,&fixed,{&interpolation}};
 ForceField forceField{&interpolation};Spring<RigidData>spring{&state.restPositions};
 Matrix system;BTD solver;DataVecDeriv forces;VecCoord freePositions;VecDeriv freeVelocities;
 double time=0;
 Kernel(){
  // InterventionalRadiologyController::init normalizes startingPos.
  controller.d_startingPos.value.q.normalize();
  spring.d_stiffness.value={parameters().spring};spring.d_angularStiffness.value={parameters().angularSpring};
  state.positions.value.assign(61,controller.d_startingPos.value);
  state.restPositions.value=state.positions.value;state.velocities.value.resize(61);
  controller.applyInterventionalRadiologyController();system.resize(366,366);
 }
 void freeMotion(double insert,double rotate){
  // interactive_xray_scene.py::onAnimateBeginEvent (same rates and clamp).
  controller.d_xTip.value[0]=std::clamp(controller.d_xTip.value[0]+insert*parameters().insertRate*h,0.,parameters().total());
  controller.d_rotationInstrument.value[0]+=rotate*parameters().rotateRate*h;
  controller.applyInterventionalRadiologyController();
  unsigned root=controller.d_indexFirstNode.value;
  spring.m_indices={root};spring.m_ext_indices={root};
  forces.value.assign(61,Vec6());
  forceField.addForce(nullptr,forces,state.positions,state.velocities);
  spring.addForce(nullptr,forces,state.positions,state.velocities);
  freePositions=state.positions.value;freeVelocities=state.velocities.value;
  if(interpolation.beams.empty()){return;}
  // EulerImplicitSolver::solve lines 145–175: compute RHTerm then MBK.
  DataVecDeriv rhs=forces;
  forceField.addMDx(nullptr,rhs,state.velocities,-parameters().alpha);
  for(unsigned b=0;b<interpolation.beams.size();b++){
   auto& e=interpolation.beams[b];
   forceField.applyStiffnessLarge(rhs.value,state.velocities.value,b,e.node0,e.node1,h+parameters().beta);
  }
  for(int j=0;j<6;j++)rhs.value[root][j]-=(j<3?parameters().spring:parameters().angularSpring)*(h+parameters().beta)*state.velocities.value[root][j];
  for(auto& v:rhs.value)v*=h;
  system.resize(366,366);
  MechanicalParams params{-h*(h+parameters().beta),1+h*parameters().alpha};
  sofa::core::behavior::MultiMatrixAccessor matrices{&system};
  forceField.addMBKToMatrix(&params,&matrices);
  for(int j=0;j<6;j++)system.add(root*6+j,root*6+j,h*(h+parameters().beta)*(j<3?parameters().spring:parameters().angularSpring));
  // FixedProjectiveConstraint: project response and matrix; active root remains spring-supported.
  for(unsigned node:fixed.indices){
   rhs.value[node].clear();
   for(int i=0;i<system.n;i++){system.blocks[node*system.n+i].clear();system.blocks[i*system.n+node].clear();}
   system.blocks[node*system.n+node].identity();
  }
  solver.invert(system);Vector b,x;
  for(int i=0;i<61;i++)b.blocks[i]=rhs.value[i];
  solver.solve(system,x,b);
  // EulerImplicitSolver and FreeMotionAnimationLoop: v_free=v+dv, q_free=q⊕h*v_free.
  for(int i=0;i<61;i++){
   freeVelocities[i]+=x.blocks[i];Vec6 displacement=freeVelocities[i]*h;
   integratePose(freePositions[i],RigidDerivativeView{displacement});
  }
 }
 // Only used by the internal no-contact conformance probe.
 void advanceWithoutContact(double insert,double rotate){
  freeMotion(insert,rotate);state.positions.value=freePositions;state.velocities.value=freeVelocities;time+=h;
 }
};
