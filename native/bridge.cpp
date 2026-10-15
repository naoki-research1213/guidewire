#include "scene.h"
#include <string>
static thread_local std::string lastError;
#define API extern "C" __declspec(dllexport)
API const char*gw_error(){return lastError.c_str();}
API void*gw_create(const char*mesh){try{return new Scene(mesh);}catch(const std::exception&e){lastError=e.what();return nullptr;}}
API void*gw_create_parameters(const char*mesh,const double*values,int count){try{
 if(!values||count!=23)throw std::invalid_argument("Expected 23 scene parameters");
 Parameters p;double*fields[]={&p.dt,&p.radius,&p.young,&p.density,&p.poisson,&p.straight,&p.tip,&p.spireDiameter,&p.insertRate,&p.rotateRate,&p.alpha,&p.beta,&p.spring,&p.angularSpring,&p.alarm,&p.contact,&p.tolerance,&p.maxIterations,&p.gravity,&p.threshold,&p.angleCone,&p.coneFactor,&p.mu};
 for(int i=0;i<count;i++){if(!std::isfinite(values[i]))throw std::invalid_argument("Non-finite scene parameter");*fields[i]=values[i];}
 if(p.dt<=0||p.radius<=0||p.young<=0||p.density<=0||p.straight<=0||p.tip<=0||p.spireDiameter<=0||p.spring<=0||p.angularSpring<=0||p.tolerance<=0||p.threshold<=0||p.contact<0||p.alarm<=p.contact||p.poisson<0||p.poisson>=.5||p.alpha<0||p.beta<0||p.mu<0||p.maxIterations<1||p.maxIterations>100000||std::floor(p.maxIterations)!=p.maxIterations)throw std::invalid_argument("Invalid scene parameter range");
 ParameterScope scope(p);return new Scene(mesh);
 }catch(const std::exception&e){lastError=e.what();return nullptr;}}
API void gw_destroy(void*p){delete static_cast<Scene*>(p);}
API int gw_step(void*p,double insert,double rotate){try{auto&s=*static_cast<Scene*>(p);ParameterScope scope(s.settings);s.advance(insert,rotate);for(auto&q:s.kernel.state.positions.value)for(unsigned i=0;i<3;i++)if(!std::isfinite(q.p[i]))throw std::runtime_error("non-finite mechanical state");return 0;}catch(const std::exception&e){lastError=e.what();return -1;}}
// 0=current Rigid3 pose, 1=rest pose, 2=velocity, 3=collision points,
// 4=visual centerline, 5=constraint forces. All buffers are caller-owned.
API int gw_read(void*p,int kind,double*out,int capacity){
 auto&s=*static_cast<Scene*>(p);int n=0;
 auto put=[&](double v){if(n<capacity)out[n]=v;n++;};
 if(kind==0||kind==1){const auto&x=kind==0?s.kernel.state.positions.value:s.kernel.state.restPositions.value;for(auto&q:x){for(unsigned i=0;i<3;i++)put(q.p[i]);for(unsigned i=0;i<4;i++)put(q.q[i]);}}
 else if(kind==2){for(auto&v:s.kernel.state.velocities.value)for(unsigned i=0;i<6;i++)put(v[i]);}
 else if(kind==3||kind==4){const auto&x=kind==3?s.wireState.positions.value:s.visualOutput.positions;for(auto&v:x)for(unsigned i=0;i<3;i++)put(v[i]);}
 else if(kind==5){for(auto f:s.forces)put(f);}
 else if(kind==6){put(s.kernel.time);put(s.kernel.controller.d_xTip.value[0]);put(s.kernel.controller.d_rotationInstrument.value[0]);put(s.kernel.controller.d_indexFirstNode.value);put(s.contactCount);put(s.residuals.size());put(s.converged);}
 else if(kind==7){put(s.collisionJacobian.rows.size());for(auto&[r,cols]:s.collisionJacobian.rows){put(r);put(cols.size());for(auto&[p,v]:cols){put(p);for(unsigned i=0;i<3;i++)put(v[i]);}}}
 return n;
}
