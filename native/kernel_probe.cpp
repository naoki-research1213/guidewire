#include "kernel.h"
#include "generated/beam_functions.inc"
#include "generated/btd_functions.inc"
#include "generated/spring_functions.inc"
#include <iomanip>
int main(){
 Kernel kernel;
 for(int i=0;i<400;i++){
  kernel.advanceWithoutContact(1.,i>200?1.:0.);
  for(auto& p:kernel.state.positions.value)for(int j=0;j<3;j++)if(!std::isfinite(p.p[j]))return 2;
  if(i%40==39){auto& p=kernel.state.positions.value.back();std::cout<<std::setprecision(17)<<i+1<<","<<kernel.controller.d_xTip.value[0]<<","<<p.p[0]<<","<<p.p[1]<<","<<p.p[2]<<"\n";}
 }
}
