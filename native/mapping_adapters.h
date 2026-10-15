#pragma once
#include "collision_adapters.h"
#include <cassert>
using InReal=double;
using sofa::type::Vec;using sofa::helper::AdvancedTimer;
namespace sofa::core {
 struct ConstraintParams{};
 using ConstVecCoordId=int;using VecCoordId=int;
 struct MultiVecCoordId {int id;MultiVecCoordId(int x):id(x){}MultiVecCoordId(vec_id::write_access::Position):id(0){}int getId(const void*)const{return id;}};
 using ConstMultiVecCoordId=MultiVecCoordId;
}
using MultiVecCoordId=core::MultiVecCoordId;using ConstMultiVecCoordId=core::ConstMultiVecCoordId;using VecCoordId=int;
struct MappingParams {core::MultiVecCoordId id{0};const core::MultiVecCoordId&x()const{return id;}};
namespace sofa::simulation {
 enum class ForEachExecutionPolicy{PARALLEL,SEQUENTIAL};struct TaskScheduler{};
 struct MainTaskSchedulerFactory {static TaskScheduler*createInRegistry(){static TaskScheduler task;return &task;}};
 template<class Iter,class Fn>void forEachRange(ForEachExecutionPolicy,TaskScheduler&,Iter a,Iter b,Fn f){struct Range{Iter start,end;};f(Range{a,b});}
}
namespace simulation=sofa::simulation;
using MappingTopology=BaseMeshTopology;
struct MappingContext {MappingTopology*topology;MappingTopology*getMeshTopology(){return topology;}};
struct MappingOutput {VectorVec3 positions;void resize(unsigned n){positions.resize(n);}};
struct Mapping {
 using InVecCoord=::VecCoord;using VecCoord=VectorVec3;using Deriv=Vec3;using InVecDeriv=::VecDeriv;using VecDeriv=VectorVec3;
 using InMatrixDeriv=SparseRows;using OutMatrixDeriv=SparsePoints;using InDeriv=Vec6;
 struct In {using MatrixDeriv=SparseRows;static Vec3 getDPos(const Vec6&v){return Vec3(v[0],v[1],v[2]);}static Vec3 getDRot(const Vec6&v){return Vec3(v[3],v[4],v[5]);}static void setDPos(Vec6&v,const Vec3&p){for(int i=0;i<3;i++)v[i]=p[i];}static void setDRot(Vec6&v,const Vec3&p){for(int i=0;i<3;i++)v[i+3]=p[i];}};
 struct Out {using MatrixDeriv=SparsePoints;};
 struct PosPointDefinition {unsigned beamId;Vec3 baryPoint;PosPointDefinition(unsigned b,const Vec3&p):beamId(b),baryPoint(p){}};
 Interpolation*l_adaptativebeamInterpolation;MechanicalState*fromModel=nullptr;MappingOutput*toModel;MechanicalState*getFromModel(){return fromModel;}
 Data<VectorVec3>d_points;Data<std::vector<double>>d_segmentsCurvAbs;
 Data<bool>d_useCurvAbs{true},d_contactDuplicate{false},d_parallelMapping{false},f_printLog{false};Data<unsigned>d_nbPointsPerBeam{0};
 Data<double>d_proximity{0.};
 bool m_isSubMapping=false,m_isBarycentricMapping=false;
 std::vector<PosPointDefinition>m_pointBeamDistribution;std::vector<unsigned>m_idPointSubMap;
 MappingContext context{nullptr};MappingContext*getContext(){return &context;}
 void computeDistribution();void computeIdxAndBaryCoordsForAbs(unsigned&,Real&,const Real&);
 void apply(const MappingParams*,Data<VecCoord>&,const Data<InVecCoord>&);
 void applyJTonPoint(unsigned,const Deriv&,SpatialVector&,SpatialVector&,const InVecCoord&);
 void applyJonPoint(unsigned,SpatialVector&,SpatialVector&,Deriv&,const InVecCoord&);
 void applyJ(const MappingParams*,Data<VecDeriv>&,const Data<InVecDeriv>&);
 void applyJT(const core::ConstraintParams*,Data<InMatrixDeriv>&,const Data<OutMatrixDeriv>&);
};
