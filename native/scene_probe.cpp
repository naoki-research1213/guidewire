#include "scene.h"
#include <iomanip>
#include <chrono>
int main(int argc,char**argv){
 if(argc!=2)return 2;Scene scene(argv[1]);auto start=std::chrono::steady_clock::now();std::cout<<std::setprecision(17);
 for(unsigned i=0;i<400;i++){scene.advance(1.,i>200?1.:0.);if(i%40==39){auto&p=scene.kernel.state.positions.value.back().p;std::cout<<i+1<<","<<scene.kernel.controller.d_xTip.value[0]<<","<<p[0]<<","<<p[1]<<","<<p[2]<<","<<scene.contactCount<<","<<scene.residuals.size()<<"\n";}}
 std::cerr<<"elapsed="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
}
