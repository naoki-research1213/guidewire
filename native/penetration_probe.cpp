// 調査専用。通常DLL・Scene・抽出した原関数は変更しない。
#include "scene.h"
#include <iomanip>
double dot3(const Vec3&a,const Vec3&b){double d=0;for(int i=0;i<3;i++)d+=a[i]*b[i];return d;}
Vec3 cross3(const Vec3&a,const Vec3&b){return Vec3(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]);}
unsigned crossings(const Scene&s){
 unsigned result=0;const auto&w=s.wireState.positions.value;const auto&v=s.vesselState.positions.value;
 for(unsigned i=1;i<w.size();i++){
  Vec3 d=w[i]-w[i-1];if(dot3(d,d)<1e-20)continue;
  for(auto&t:s.vesselTopology.triangles){
   Vec3 e1=v[t[1]]-v[t[0]],e2=v[t[2]]-v[t[0]],p=cross3(d,e2);double det=dot3(e1,p);
   if(std::abs(det)<1e-12)continue;
   Vec3 delta=w[i-1]-v[t[0]];double u=dot3(delta,p)/det;if(u<0||u>1)continue;
   Vec3 q=cross3(delta,e1);double b=dot3(d,q)/det;if(b<0||u+b>1)continue;
   double a=dot3(e2,q)/det;if(a>1e-8&&a<1-1e-8){result++;break;}
  }
 }
 return result;
}
int main(int argc,char**argv){
 if(argc!=5)return 2;
 Parameters p;p.dt=.002;p.alpha=0;p.beta=.2;p.radius=1.08;
 p.young=7.5e10/6000.;p.density=6500./(6000.*6000.*6000.);
 p.gravity=-9.81*6000.;p.contact=1.08;p.maxIterations=5000;
 p.insertRate=std::stod(argv[2]);p.mu=std::stod(argv[3]);
 ParameterScope scope(p);Scene s(argv[1]);int steps=std::stoi(argv[4]);
 std::cout<<std::setprecision(17)<<"step,time,input,xtip,contacts,converged,iterations,min_pred_gap,min_actual_gap,max_gap_mismatch,max_normal_impulse,max_speed,crossing_edges\n";
 for(int i=0;i<steps;i++){
  double input=i<steps*3/4?1.:0.;s.advance(input,0.);
  double predMin=0,actualMin=0,mismatch=0,impulse=0,speed=0;
  for(auto&v:s.kernel.state.velocities.value){double n=0;for(int k=0;k<3;k++)n+=v[k]*v[k];speed=std::max(speed,std::sqrt(n));}
  for(auto&[key,owned]:s.responses){
   auto&r=*owned;VectorVec3 q=r.mapper1.output.positions.value,pnew=r.mapper2.output.positions.value;
   r.mapper1.apply(q,r.model1->state->positions.value);
   r.mapper2.apply(pnew,r.model2->state->positions.value);
   for(auto&c:r.constraint.contacts){
    double predicted=c.dfree;
    auto row=s.constraintState.jacobian.value.rows.find(c.id);
    if(row!=s.constraintState.jacobian.value.rows.end())for(auto&[node,j]:row->second){
     for(int k=0;k<6;k++)predicted+=s.kernel.h*j[k]*(s.kernel.state.velocities.value[node][k]-s.kernel.freeVelocities[node][k]);
    }
    double actual=dot3(c.norm,pnew[c.m2]-q[c.m1])-c.contactDistance;
    predMin=std::min(predMin,predicted);actualMin=std::min(actualMin,actual);
    mismatch=std::max(mismatch,std::abs(actual-predicted));
    if(c.id<s.forces.size())impulse=std::max(impulse,s.forces[c.id]);
   }
  }
  int cross=i%25==24||i==steps-1?static_cast<int>(crossings(s)):-1;
  std::cout<<i+1<<','<<s.kernel.time<<','<<input<<','<<s.kernel.controller.d_xTip.value[0]<<','<<s.contactCount<<','<<s.converged<<','<<s.residuals.size()<<','<<predMin<<','<<actualMin<<','<<mismatch<<','<<impulse<<','<<speed<<','<<cross<<'\n';
 }
}
