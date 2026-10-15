#pragma once
#include "constraint_adapters.h"
#include <sofa/type/Mat.h>
#define EMIT_EXTRA_DEBUG_MESSAGE false
using Matrix2=sofa::type::Mat<2,2,double>;
using sofa::type::MatNoInit;using sofa::type::VecNoInit;
namespace sofa::defaulttype {struct Vec3Types{};}
namespace sofa::simulation {
 struct Context {virtual ~Context()=default;};
 struct Node:Context {void*line=nullptr;template<class T>T*get()const{return static_cast<T*>(line);}};
}
namespace sofa::core::vec_id::read_access {inline auto position=sofa::core::vec_id::write_access::position;}
struct CollisionState {Data<VectorVec3>positions;unsigned getSize(){return positions.value.size();}Data<VectorVec3>*read(core::vec_id::write_access::Position){return &positions;}};
namespace sofa::core {struct CollisionElementIterator;}
struct CollisionModel {
 virtual ~CollisionModel()=default;
 struct SPtr {CollisionModel*p;SPtr(CollisionModel*x):p(x){}CollisionModel*operator->(){return p;}operator CollisionModel*(){return p;}void reset(){p=nullptr;}};
 CollisionState*state;BaseMeshTopology*topology;sofa::simulation::Node*node;
 Data<bool>d_bothSide{false};Data<double>d_contactDistance{0.};unsigned size=0;int kind=-1;bool simulated=true,moving=true;
 std::unique_ptr<CollisionModel>previous;CollisionModel*next=nullptr;
 sofa::simulation::Context*getContext(){return node;}BaseMeshTopology*getCollisionTopology(){return topology;}
 CollisionState*getMechanicalState(){return state;}unsigned getSize(){return size;}
 bool empty(){return size==0;}bool isSimulated(){return simulated;}
 double getContactStiffness(unsigned){return 10.;}double getContactDistance(){return d_contactDistance.value;}
 bool isMoving(){return moving;}
 virtual void resize(unsigned n){size=n;}
 CollisionModel*getPrevious(){return previous.get();}CollisionModel*getNext(){return next;}
 CollisionModel*getLast(){auto*m=this;while(m->next)m=m->next;return m;}
 void setPrevious(CollisionModel*){throw std::logic_error("dynamic BVH depth changes not extracted");}
 CollisionModel*getMaster(){return nullptr;}void removeSlave(SPtr){}
 template<class T>T*createPrevious(){if(!previous){previous=std::make_unique<T>();previous->next=this;previous->simulated=simulated;previous->node=node;}return static_cast<T*>(previous.get());}
 const char*getClassName(){return "ExtractedCollisionModel";}
 core::CollisionElementIterator begin();core::CollisionElementIterator end();
 virtual std::pair<core::CollisionElementIterator,core::CollisionElementIterator>getInternalChildren(unsigned)const;
 virtual std::pair<core::CollisionElementIterator,core::CollisionElementIterator>getExternalChildren(unsigned)const;
};
struct Triangle;
template<class T>struct PointCollisionModel:CollisionModel {CollisionState*mstate=nullptr;Data<bool>d_computeNormals{false};PointCollisionModel(){kind=0;}void updateNormals(){}void computeBoundingTree(int);};
template<class T>struct LineCollisionModel:CollisionModel {CollisionState*mstate=nullptr;struct LineData{unsigned p[2];};std::vector<LineData>elems;bool needsUpdate=true;LineCollisionModel(){kind=1;}
 void updateFromTopology(){if(elems.size()!=topology->edges.size()){size=topology->edges.size();elems.resize(size);for(unsigned i=0;i<size;i++){elems[i].p[0]=topology->edges[i][0];elems[i].p[1]=topology->edges[i][1];}needsUpdate=true;}}
 void computeBoundingTree(int);};
template<class T>struct TriangleCollisionModel:CollisionModel {using VecCoord=VectorVec3;using Element=Triangle;CollisionState*m_mstate=nullptr;BaseMeshTopology*m_topology=nullptr;int m_topologyRevision=0;bool m_needsUpdate=true;Data<bool>d_computeNormals{true},d_useCurvature{false};std::vector<Vec3>normals;
 TriangleCollisionModel(){kind=2;}void updateFromTopology(){size=topology->triangles.size();normals.resize(size);m_needsUpdate=true;}
 void computeBoundingTree(int);};
using PointModel=PointCollisionModel<sofa::defaulttype::Vec3Types>;
using LineModel=LineCollisionModel<sofa::defaulttype::Vec3Types>;
using TriangleModel=TriangleCollisionModel<sofa::defaulttype::Vec3Types>;
struct Point {
 using Coord=Vec3;PointModel*model;unsigned index;
 PointModel*getCollisionModel(){return model;}unsigned getIndex(){return index;}
 const Vec3&p(){return model->state->positions.value[index];}
 bool isActive(CollisionModel*){return true;}double getContactDistance(){return 0.;}
};
template<class T>using TPoint=Point;
struct Line {
 using Coord=Vec3;LineModel*model;unsigned index;
 LineModel*getCollisionModel(){return model;}unsigned getIndex(){return index;}
 unsigned i1(){return model->topology->edges[index][0];}unsigned i2(){return model->topology->edges[index][1];}
 const Vec3&p1(){return model->state->positions.value[i1()];}const Vec3&p2(){return model->state->positions.value[i2()];}
 bool isActive(CollisionModel*){return true;}double getContactDistance(){return 0.;}
};
struct Triangle {
 using Coord=Vec3;TriangleModel*model;unsigned index;
 TriangleModel*getCollisionModel(){return model;}unsigned getIndex(){return index;}
 unsigned p1Index(){return model->topology->triangles[index][0];}unsigned p2Index(){return model->topology->triangles[index][1];}unsigned p3Index(){return model->topology->triangles[index][2];}
 Vec3&n(){return model->normals[index];}
 const Vec3&p1(){return model->state->positions.value[model->topology->triangles[index][0]];}
 const Vec3&p2(){return model->state->positions.value[model->topology->triangles[index][1]];}
 const Vec3&p3(){return model->state->positions.value[model->topology->triangles[index][2]];}
 double getContactDistance(){return 0.;}
};
namespace sofa::core {
 using CollisionModel=::CollisionModel;
 struct CollisionElementIterator {CollisionModel*model=nullptr;unsigned index=0;
 CollisionElementIterator()=default;template<class Element>CollisionElementIterator(Element e):model(e.model),index(e.index){}
 CollisionElementIterator(CollisionModel*m,unsigned i):model(m),index(i){}
 unsigned getIndex()const{return index;}CollisionModel*getCollisionModel()const{return model;}
 CollisionElementIterator&operator++(){++index;return *this;}CollisionElementIterator operator+(int n)const{return {model,index+n};}
 bool operator==(const CollisionElementIterator& b)const{return model==b.model&&index==b.index;}bool operator!=(const CollisionElementIterator& b)const{return !(*this==b);}
 bool valid()const{return model!=nullptr;}bool canCollideWith(const CollisionElementIterator&)const{return false;}
 std::pair<CollisionElementIterator,CollisionElementIterator>getInternalChildren()const{return model?model->getInternalChildren(index):std::pair<CollisionElementIterator,CollisionElementIterator>{};}
 std::pair<CollisionElementIterator,CollisionElementIterator>getExternalChildren()const{return model?model->getExternalChildren(index):std::pair<CollisionElementIterator,CollisionElementIterator>{};}
 };
 namespace collision {struct Intersection {double getAlarmDistance()const{return parameters().alarm;}double getContactDistance()const{return parameters().contact;}};}
}
inline core::CollisionElementIterator CollisionModel::begin(){return {this,0};}
inline core::CollisionElementIterator CollisionModel::end(){return {this,size};}
inline std::pair<core::CollisionElementIterator,core::CollisionElementIterator>CollisionModel::getInternalChildren(unsigned)const{return {};}
inline std::pair<core::CollisionElementIterator,core::CollisionElementIterator>CollisionModel::getExternalChildren(unsigned)const{return {};}
struct DetectionOutput {using ContactId=int64_t;std::pair<core::CollisionElementIterator,core::CollisionElementIterator>elem;ContactId id;Vec3 point[2],normal;double value;};
namespace sofa::core::collision {using DetectionOutput=::DetectionOutput;}
struct LocalMinDistance {
 using OutputVector=std::vector<DetectionOutput>;
 Data<bool>d_useLMDFilters{false},d_filterIntersection{true};Data<double>d_angleCone{parameters().angleCone},d_coneFactor{parameters().coneFactor};
 int computeIntersection(Triangle&,Point&,OutputVector*,const core::collision::Intersection*);
 int computeIntersection(Line&,Point&,OutputVector*,const core::collision::Intersection*);
 int computeIntersection(Line&,Line&,OutputVector*,const core::collision::Intersection*);
 bool testValidity(Point&,const Vec3&)const;bool testValidity(Line&,const Vec3&)const;bool testValidity(Triangle&,const Vec3&)const;
};
