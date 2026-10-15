#include "implementation.h"
CollisionModel* first(CollisionModel*m){while(m->getPrevious())m=m->getPrevious();return m;}
int main(int argc,char**argv){
 if(argc!=2)return 2;
 MeshSTLLoader loader;if(!loader.readBinarySTL(argv[1]))return 3;
 MeshTopology vessel;vessel.triangles=loader.d_triangles.value;vessel.nbPoints=loader.d_positions.value.size();vessel.connect();
 CollisionState vs;vs.positions.value=loader.d_positions.value;
 sofa::simulation::Node vn,wn;
 TriangleModel tri;LineModel vl;PointModel wp;LineModel wl;
 vn.line=&vl;wn.line=&wl;
 for(CollisionModel*m:{static_cast<CollisionModel*>(&tri),static_cast<CollisionModel*>(&vl)}){m->state=&vs;m->topology=&vessel;m->node=&vn;m->simulated=false;m->moving=false;}
 tri.m_mstate=&vs;tri.m_topology=&vessel;tri.updateFromTopology();vl.mstate=&vs;
 MeshTopology wire;wire.nbPoints=2;wire.d_seqEdges.value={{0,1}};wire.connect();
 CollisionState ws;auto t=vessel.triangles[0];Vec3 p=(vs.positions.value[t[0]]+vs.positions.value[t[1]]+vs.positions.value[t[2]])/3.;
 Vec3 n=cross(vs.positions.value[t[1]]-vs.positions.value[t[0]],vs.positions.value[t[2]]-vs.positions.value[t[0]]);n.normalize();
 ws.positions.value={p+n*.5,p+n*.5+Vec3(.1,0,0)};
 for(CollisionModel*m:{static_cast<CollisionModel*>(&wp),static_cast<CollisionModel*>(&wl)}){m->state=&ws;m->topology=&wire;m->node=&wn;}
 wp.mstate=&ws;wl.mstate=&ws;
 tri.computeBoundingTree(6);vl.computeBoundingTree(6);wp.computeBoundingTree(6);wl.computeBoundingTree(6);
 LocalMinDistance lmd;IntersectionMethod intersection(&lmd);BVHNarrowPhase narrow;narrow.intersectionMethod=&intersection;
 narrow.addCollisionPair({first(&tri),first(&wp)});narrow.addCollisionPair({first(&vl),first(&wp)});narrow.addCollisionPair({first(&vl),first(&wl)});
 size_t count=0;for(auto&[key,out]:narrow.outputs)if(out)count+=out->size();
 std::cout<<"vertices="<<vs.getSize()<<" triangles="<<tri.size<<" edges="<<vl.size<<" contacts="<<count<<"\n";
 return count?0:4;
}
