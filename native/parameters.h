#pragma once
#include <cmath>
#include <stdexcept>
// シーン設定のみ。抽出した原関数の計算方式は変更しない。
struct Parameters {
 double dt=.005,radius=.9,young=20000.,density=.00000155,poisson=.3;
 double straight=980.,tip=20.,spireDiameter=25.,insertRate=10.,rotateRate=.5;
 double alpha=.1,beta=.2,spring=1e8,angularSpring=1e8,alarm=2.,contact=.1;
 double tolerance=.0001,maxIterations=1000.,gravity=-9.81,threshold=.01,angleCone=.2,coneFactor=.5,mu=.15;
 double total()const{return straight+tip;}
};
inline const Parameters defaultParameters;
inline thread_local const Parameters* activeParameters=&defaultParameters;
inline const Parameters& parameters(){return *activeParameters;}
struct ParameterScope {
 const Parameters* previous;
 explicit ParameterScope(const Parameters& p):previous(activeParameters){activeParameters=&p;}
 ~ParameterScope(){activeParameters=previous;}
};
