#pragma once
#include "collision_adapters.h"
#include <queue>
#include <stack>
struct CubeCollisionModel;
struct Cube:core::CollisionElementIterator {
 Cube(CubeCollisionModel*model=nullptr,unsigned index=0);
 explicit Cube(const core::CollisionElementIterator& i):core::CollisionElementIterator(i){}
 const Vec3&minVect()const;const Vec3&maxVect()const;
 const std::pair<Cube,Cube>&subcells()const;double getConeAngle()const;const Vec3&getConeAxis()const;
 double getContactDistance(){return 0.;}
};
struct CubeCollisionModel:CollisionModel {
 struct CubeData {Vec3 minBBox,maxBBox;std::pair<Cube,Cube>subcells;std::pair<core::CollisionElementIterator,core::CollisionElementIterator>children;Vec3 coneAxis;double coneAngle=0;};
 struct CubeSortPredicate {int axis;explicit CubeSortPredicate(int a):axis(a){}bool operator()(const CubeData&a,const CubeData&b)const{return a.minBBox[axis]+a.maxBBox[axis]<b.minBBox[axis]+b.maxBBox[axis];}};
 std::vector<CubeData>elems;std::vector<unsigned>parentOf;
 CubeCollisionModel(){kind=3;}
 void resize(unsigned)override;void setParentOf(unsigned,const Vec3&,const Vec3&);
 void setParentOf(unsigned,const Vec3&,const Vec3&,const Vec3&,const SReal angle=0);
 int addCube(Cube,Cube);void updateCube(unsigned);void updateCubes();void computeBoundingTree(int);
 std::pair<core::CollisionElementIterator,core::CollisionElementIterator>getInternalChildren(unsigned)const override;
 std::pair<core::CollisionElementIterator,core::CollisionElementIterator>getExternalChildren(unsigned)const override;
};
inline Cube::Cube(CubeCollisionModel*m,unsigned i):core::CollisionElementIterator(m,i){}
inline const Vec3&Cube::minVect()const{return static_cast<CubeCollisionModel*>(model)->elems[index].minBBox;}
inline const Vec3&Cube::maxVect()const{return static_cast<CubeCollisionModel*>(model)->elems[index].maxBBox;}
inline const std::pair<Cube,Cube>&Cube::subcells()const{return static_cast<CubeCollisionModel*>(model)->elems[index].subcells;}
inline double Cube::getConeAngle()const{return static_cast<CubeCollisionModel*>(model)->elems[index].coneAngle;}
inline const Vec3&Cube::getConeAxis()const{return static_cast<CubeCollisionModel*>(model)->elems[index].coneAxis;}
namespace sofa::core::collision {
 using DetectionOutputVector=std::vector<DetectionOutput>;
 struct ElementIntersector {
  virtual ~ElementIntersector()=default;
  virtual bool canIntersect(core::CollisionElementIterator,core::CollisionElementIterator,const Intersection*)=0;
  virtual void intersect(core::CollisionElementIterator,core::CollisionElementIterator,DetectionOutputVector*,const Intersection*)=0;
  void beginIntersect(CollisionModel*,CollisionModel*,DetectionOutputVector*& out){if(!out)out=new DetectionOutputVector();}
 };
}
struct MirrorIntersector:core::collision::ElementIntersector {
 core::collision::ElementIntersector*intersector;
 bool canIntersect(core::CollisionElementIterator a,core::CollisionElementIterator b,const core::collision::Intersection*i)override{return intersector->canIntersect(b,a,i);}
 void intersect(core::CollisionElementIterator a,core::CollisionElementIterator b,core::collision::DetectionOutputVector*out,const core::collision::Intersection*i)override{intersector->intersect(b,a,out,i);}
};
struct LMDIntersector:core::collision::ElementIntersector {
 LocalMinDistance*lmd;
 bool canIntersect(core::CollisionElementIterator,core::CollisionElementIterator,const core::collision::Intersection*)override;
 void intersect(core::CollisionElementIterator,core::CollisionElementIterator,core::collision::DetectionOutputVector*,const core::collision::Intersection*)override;
 bool testCubes(Cube&,Cube&,const core::collision::Intersection*);
};
struct IntersectionMethod:core::collision::Intersection {
 LMDIntersector intersector;explicit IntersectionMethod(LocalMinDistance*lmd){intersector.lmd=lmd;}
 core::collision::ElementIntersector*findIntersector(CollisionModel*a,CollisionModel*b,bool&swap){
  swap=false;if(a->kind==3&&b->kind==3)return &intersector;
  if((a->kind==2&&b->kind==0)||(a->kind==1&&b->kind==0)||(a->kind==1&&b->kind==1))return &intersector;
  if((b->kind==2&&a->kind==0)||(b->kind==1&&a->kind==0)){swap=true;return &intersector;}return nullptr;
 }
 const char*getName(){return "LocalMinDistance";}
};
struct BVHNarrowPhase {
 using CollisionIteratorRange=std::pair<core::CollisionElementIterator,core::CollisionElementIterator>;
 using TestPair=std::pair<CollisionIteratorRange,CollisionIteratorRange>;
 struct FinestCollision {CollisionModel*cm1,*cm2;core::collision::ElementIntersector*intersector;bool selfCollision;};
 IntersectionMethod*intersectionMethod;
 std::map<std::pair<CollisionModel*,CollisionModel*>,core::collision::DetectionOutputVector*>outputs;
 ~BVHNarrowPhase(){for(auto& [key,p]:outputs)delete p;}
 auto&getDetectionOutputs(CollisionModel*a,CollisionModel*b){return outputs[{a,b}];}
 static bool isSelfCollision(CollisionModel*,CollisionModel*);
 void addCollisionPair(const std::pair<CollisionModel*,CollisionModel*>&);
 static void initializeExternalCells(CollisionModel*,CollisionModel*,std::queue<TestPair>&);
 void processExternalCell(const TestPair&,CollisionModel*&,CollisionModel*&,core::collision::ElementIntersector*,const FinestCollision&,MirrorIntersector*,std::queue<TestPair>&,core::collision::DetectionOutputVector*&)const;
 static void processInternalCell(const TestPair&,core::collision::ElementIntersector*,const FinestCollision&,std::queue<TestPair>&,std::stack<TestPair>&,core::collision::DetectionOutputVector*&,const core::collision::Intersection*);
 static void visitCollisionElements(const TestPair&,core::collision::ElementIntersector*,const FinestCollision&,std::queue<TestPair>&,std::stack<TestPair>&,core::collision::DetectionOutputVector*&,const core::collision::Intersection*);
 static void visitExternalChildren(const core::CollisionElementIterator&,const core::CollisionElementIterator&,core::collision::ElementIntersector*,const FinestCollision&,std::queue<TestPair>&,core::collision::DetectionOutputVector*&,const core::collision::Intersection*);
 static void finalCollisionPairs(const TestPair&,bool,core::collision::ElementIntersector*,core::collision::DetectionOutputVector*&,const core::collision::Intersection*);
 static std::pair<CollisionModel*,CollisionModel*>getCollisionModelsFromTestPair(const TestPair&);
 static bool isRangeEmpty(const CollisionIteratorRange&);
};
