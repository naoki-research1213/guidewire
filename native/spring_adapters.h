#pragma once
#include "adapters.h"
#include <array>
namespace sofa::core::objectmodel {enum class ComponentState{Valid,Invalid};}
namespace sofa::type {template<class T>inline constexpr bool isRigidType=true;}
struct RigidData {using VecCoord=::VecCoord;using VecDeriv=::VecDeriv;using Deriv=Vec6;};
template<class D>struct ReadAccessor;
template<class T>struct ReadAccessor<Data<T>>{const T& value;ReadAccessor(const Data<T>& d):value(d.value){}const auto& operator[](size_t i)const{return value[i];}size_t size()const{return value.size();}};
struct ForceSlice {Vec6& v;int offset;void operator-=(const Vec3& x){for(int i=0;i<3;i++)v[i+offset]-=x[i];}};
inline ForceSlice getVCenter(Vec6& v){return {v,0};}
inline ForceSlice getVOrientation(Vec6& v){return {v,3};}
inline Vec3 getVCenter(const Vec6& v){return Vec3(v[0],v[1],v[2]);}
inline Vec3 getVOrientation(const Vec6& v){return Vec3(v[3],v[4],v[5]);}
template<class DataTypes>struct Spring {
 using VecCoord=typename DataTypes::VecCoord;using VecDeriv=typename DataTypes::VecDeriv;
 using DataVecCoord=Data<VecCoord>;using DataVecDeriv=Data<VecDeriv>;using Deriv=typename DataTypes::Deriv;
 using VecReal=std::vector<Real>;using CPos=Vec3;static constexpr size_t spatial_dimensions=3;
 DataVecCoord* rest;
 Data<sofa::core::objectmodel::ComponentState>d_componentState{sofa::core::objectmodel::ComponentState::Valid};
 Data<bool>d_recompute_indices{false};bool useRestMState=true;
 std::vector<unsigned>m_indices{60},m_ext_indices{60};std::vector<Vec3>m_pivots;
 Data<std::array<bool,6>>d_activeDirections{{true,true,true,true,true,true}};
 Data<VecReal>d_stiffness{{1e8}},d_angularStiffness{{1e8}};
 const DataVecCoord* getExtPosition(){return rest;}void recomputeIndices(){}
 void addForce(const MechanicalParams*,DataVecDeriv&,const DataVecCoord&,const DataVecDeriv&);
};
