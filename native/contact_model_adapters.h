#pragma once
#include "barycentric_adapters.h"
#include "generated/contact_identifier.inc"
namespace sofa::helper {
 template<class T>using WriteAccessor=::WriteAccessor<T>;
 inline double rabs(double x){return std::abs(x);}
}
struct ContactPointState {
 Data<VectorVec3>positions,freePositions;void resize(unsigned n){positions.value.resize(n);freePositions.value.resize(n);}
 unsigned getSize(){return positions.value.size();}
 Data<VectorVec3>*read(core::vec_id::read_access::FreePosition){return &freePositions;}
 Data<VectorVec3>*read(core::vec_id::write_access::Position){return &positions;}
 Data<VectorVec3>*read(core::vec_id::write_access::RestPosition){return &positions;}
};
namespace sofa::core::vec_id::read_access {inline auto restPosition=sofa::core::vec_id::write_access::restPosition;}
struct SubsetMapper {
 using IndexArray=std::vector<unsigned>;using InVecCoord=VectorVec3;using OutVecCoord=VectorVec3;
 using InDataVecCoord=Data<VectorVec3>;using OutDataVecCoord=Data<VectorVec3>;
 using InMatrixDeriv=SparsePoints;using OutMatrixDeriv=SparsePoints;using InDataMatrixDeriv=Data<SparsePoints>;using OutDataMatrixDeriv=Data<SparsePoints>;
 using In=PointData;using Out=PointData;
 Data<IndexArray>d_indices;Data<bool>d_resizeToModel{false};CollisionState*fromModel;ContactPointState*toModel;
 void apply(const MappingParams*,OutDataVecCoord&,const InDataVecCoord&);
 void applyJT(const core::ConstraintParams*,InDataMatrixDeriv&,const OutDataMatrixDeriv&);
};
struct FrictionParameters {double mu=.8;bool hasTangentialComponent()const{return mu>0.;}};
struct ContactConstraint;
struct ConstraintBlockInfo {ContactConstraint*parent=nullptr;int const0=0,nbLines=0,offsetId=0,nbGroups=0;bool hasId=false;};
using VecConstraintBlockInfo=std::vector<ConstraintBlockInfo>;using VecPersistentID=std::vector<long long>;
namespace sofa::core::behavior {using BaseLagrangianConstraint=::ContactConstraint;}
struct ContactConstraint {
 using ContactParams=FrictionParameters;using Deriv=Vec3;using Coord=Vec3;using VecCoord=VectorVec3;using PersistentID=long long;
 using DataVecCoord=Data<VectorVec3>;using DataMatrixDeriv=Data<SparsePoints>;using MatrixDerivRowIterator=SparsePoints::RowIterator;
 struct Contact {int m1,m2;Deriv norm,t,s;Real contactDistance;unsigned id;long contactId;PersistentID localId;ContactParams parameters;Coord P,Q;mutable Real dfree;};
 std::vector<Contact>contacts;ContactPointState*mstate1,*mstate2;bool yetIntegrated=false;
 ContactPointState*getMState1(){return mstate1;}ContactPointState*getMState2(){return mstate2;}
 void clear(int n){contacts.clear();contacts.reserve(n);}void setName(std::string){}void setCustomTolerance(double){}
 void addContact(const ContactParams&,Deriv,Coord,Coord,Real,int,int,long,PersistentID);
 void addContact(const ContactParams&,Deriv,Real,int,int,long,PersistentID);
 void addContact(const ContactParams&,Deriv,Coord,Coord,Real,int,int,Coord,Coord,long,PersistentID);
 void buildConstraintMatrix(const core::ConstraintParams*,DataMatrixDeriv&,DataMatrixDeriv&,unsigned&,const DataVecCoord&,const DataVecCoord&);
 void getPositionViolation(linearalgebra::BaseVector*);
 void getConstraintInfo(const core::ConstraintParams*,VecConstraintBlockInfo&,VecPersistentID&);
};
// Scene-graph allocation is replaced by owned data. Barycentric calculations
// below call the extracted BarycentricMapperMeshTopology bodies.
struct ContactMapper {
 CollisionModel*model;VectorVec3*freeInput;ContactPointState output;BaryMapper bary;SubsetMapper subset;
 ContactPointState*createMapping(const char*){return &output;}
 void resize(unsigned n){output.resize(n);bary.m_map1d.clear();bary.m_map2d.clear();bary.m_map3d.clear();subset.d_indices.value.clear();subset.fromModel=model->state;subset.toModel=&output;}
 unsigned addPointB(const Vec3&p,unsigned index,double&){
  if(model->kind==0){subset.d_indices.value.push_back(index);return subset.d_indices.value.size()-1;}
  if(model->kind==1)return bary.createPointInLine(p,index,&model->state->positions.value);
  if(model->kind==2)return bary.createPointInTriangle(p,index,&model->state->positions.value);
  throw std::logic_error("unsupported contact geometry");
 }
 void apply(VectorVec3&out,const VectorVec3&in){if(model->kind==0){Data<VectorVec3>dataOut{out},dataIn{in};subset.apply(nullptr,dataOut,dataIn);out=std::move(dataOut.value);}else bary.apply(out,in);}
 void update(){apply(output.positions.value,model->state->positions.value);}void updateXfree(){apply(output.freePositions.value,*freeInput);}
 void applyJT(SparsePoints&out,const SparsePoints&in){if(model->kind==0){Data<SparsePoints>dataOut{out},dataIn{in};subset.applyJT(nullptr,dataOut,dataIn);out=std::move(dataOut.value);}else bary.applyJT(out,in);}
};
namespace sofa::core::objectmodel {struct BaseContext {template<class T>void removeObject(T*){}template<class T>void addObject(T*){}};}
struct ContactResponse:sofa::component::collision::response::contact::ContactIdentifier {
 using TOutputVector=std::vector<DetectionOutput>;using OutputVector=TOutputVector;
 using MechanicalState1=ContactPointState;using MechanicalState2=ContactPointState;using CollisionElement1=core::CollisionElementIterator;using CollisionElement2=core::CollisionElementIterator;
 using DataTypes1=PointData;using DataTypes2=PointData;using ConstraintParameters=FrictionParameters;
 CollisionModel*model1,*model2;core::collision::Intersection*intersectionMethod;
 std::vector<DetectionOutput*>contacts;std::vector<std::pair<std::pair<int,int>,double>>mappedContacts;
 ContactMapper mapper1,mapper2;ContactConstraint constraint;ContactConstraint*m_constraint=&constraint;
 bool selfCollision=false;Data<double>d_tol{0.};core::objectmodel::BaseContext*parent=nullptr;
 std::string getName(){return "FrictionContactConstraint";}
 void setupConstraint(ContactPointState*a,ContactPointState*b){constraint.mstate1=a;constraint.mstate2=b;}
 void setInteractionTags(ContactPointState*,ContactPointState*){}
 FrictionParameters getParameterFromDatas(){return {};}
 static long cantorPolynomia(int64_t x,int64_t y){return sofa::component::collision::response::contact::cantorPolynomia(x,y);}
 void setDetectionOutputs(OutputVector*);void activateMappers();void createResponse(core::objectmodel::BaseContext*);
};
