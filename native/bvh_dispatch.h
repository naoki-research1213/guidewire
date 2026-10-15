#pragma once
#include "bvh_adapters.h"
inline bool LMDIntersector::canIntersect(core::CollisionElementIterator a,core::CollisionElementIterator b,const core::collision::Intersection*i){
 if(a.model->kind!=3||b.model->kind!=3)throw std::logic_error("Unexpected coarse intersector pair");
 Cube x(a),y(b);return testCubes(x,y,i);
}
inline void LMDIntersector::intersect(core::CollisionElementIterator a,core::CollisionElementIterator b,core::collision::DetectionOutputVector*out,const core::collision::Intersection*i){
 if(a.model->kind==2&&b.model->kind==0){Triangle t{static_cast<TriangleModel*>(a.model),a.index};Point p{static_cast<PointModel*>(b.model),b.index};lmd->computeIntersection(t,p,out,i);}
 else if(a.model->kind==1&&b.model->kind==0){Line l{static_cast<LineModel*>(a.model),a.index};Point p{static_cast<PointModel*>(b.model),b.index};lmd->computeIntersection(l,p,out,i);}
 else if(a.model->kind==1&&b.model->kind==1){Line l{static_cast<LineModel*>(a.model),a.index},r{static_cast<LineModel*>(b.model),b.index};lmd->computeIntersection(l,r,out,i);}
 else throw std::logic_error("Unexpected finest intersector pair");
}
