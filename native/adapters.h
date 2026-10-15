#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
// Data/container adapters only; original SOFA math remains in upstream/include.
#include <vector>
#include <limits>
#include <iostream>
#include <algorithm>
#include <array>
#include "parameters.h"
#include <sofa/type/Quat.inl>
#include <sofa/type/SpatialVector.inl>
#include <sofa/type/Transform.inl>
#define SOFA_UNUSED(x) (void)(x)
#define SCOPED_TIMER(x) (void)0
struct Quiet { template<class T> Quiet& operator<<(const T&) {return *this;} };
#define dmsg_info() Quiet()
#define msg_warning() Quiet()
#define msg_error(...) Quiet()
#define msg_warning_when(x) Quiet()
#define dmsg_error() Quiet()
using Real=double;
using Vec3=sofa::type::Vec<3,double>;
using Vec2=sofa::type::Vec<2,double>;
using Vec6=sofa::type::Vec<6,double>;
using Vec3NoInit=sofa::type::VecNoInit<3,double>;
using Vec6NoInit=sofa::type::VecNoInit<6,double>;
using Matrix6x6=sofa::type::Mat<6,6,double>;
using Quat=sofa::type::Quat<double>;
using Transform=sofa::type::Transform<double>;
using SpatialVector=sofa::type::SpatialVector<double>;
namespace type=sofa::type;
namespace sofa::type {template<class T>using vector=std::vector<T>;}
namespace sofa::helper {struct AdvancedTimer{struct StepVar{explicit StepVar(const char*){}};static void stepBegin(const char*){}static void stepEnd(const char*){}static void valSet(const char*,int){}};}
struct Pose {Vec3 p;Quat q;Vec3& getCenter(){return p;} const Vec3& getCenter()const{return p;}
 Quat& getOrientation(){return q;} const Quat& getOrientation()const{return q;} };
using VecCoord=std::vector<Pose>; using VecDeriv=std::vector<Vec6>;
template<class T> struct Data {T value; const T& getValue()const{return value;} T& getValue(){return value;}
 void setValue(const T& x){value=x;} Data& operator=(const T& x){value=x;return *this;} T*beginEdit(){return &value;}void endEdit(){} };
using DataVecDeriv=Data<VecDeriv>;using DataVecCoord=Data<VecCoord>;
template<class T>struct Access {T* ptr;auto& operator[](size_t i){return (*ptr)[i];} auto size()const{return ptr->size();}
 void resize(size_t n){ptr->resize(n);}void clear(){ptr->clear();}void reserve(size_t n){ptr->reserve(n);}template<class V>void push_back(V&& v){ptr->push_back(std::forward<V>(v));}T*operator->(){return ptr;}T& wref(){return *ptr;}T& ref(){return *ptr;} };
namespace sofa::helper {template<class T>Access<T> getWriteOnlyAccessor(Data<T>& d){return {&d.value};}
 template<class T>const T& getReadAccessor(const Data<T>& d){return d.value;} }
namespace helper=sofa::helper;
namespace sofa::core {using ConstVecCoordId=int;}
using VectorVec3=std::vector<Vec3>;
struct BezierState {Data<VectorVec3>positions[2];const Data<VectorVec3>*read(int id){return &positions[id];}};
struct MechanicalParams {double k=0,m=0;double kFactor()const{return k;}double mFactor()const{return m;}};
struct ScalarMatrix {virtual ~ScalarMatrix()=default;virtual void add(int,int,double)=0;virtual int rowSize()const=0;virtual int colSize()const=0;virtual double element(int,int)const=0;};
namespace sofa::core::behavior {struct MultiMatrixAccessor {
 struct MatrixRef {ScalarMatrix*matrix;unsigned offset=0;};
 ScalarMatrix*matrix;MatrixRef getMatrix(const void*)const{return {matrix,0};}
};}
struct Element {unsigned node0,node1;double length,x0,x1;Transform local0,local1;};
struct RestShape {
 static constexpr Real EPSILON=.0000000001;
 struct Section {RestShape*owner;bool curved;void getRestTransformOnX(Transform&t,Real x,Real start){if(curved)owner->spire(t,x,start);else owner->straight(t,x,start);}};
 struct Sections {std::array<Section,2>values;size_t size(){return values.size();}Section*get(unsigned i){return &values[i];}}l_sectionMaterials{{{{this,false},{this,true}}}};
 Data<std::vector<Real>>d_keyPoints{{0.,parameters().straight,parameters().total()}};
 Data<double> d_spireDiameter{parameters().spireDiameter},d_spireHeight{0.};
 void spire(Transform&,Real,Real);
 void straight(Transform&,Real,Real);Real getLength(){return parameters().total();}
 void getRestTransformOnX(Transform&,Real);
};
struct Interpolation {
 static constexpr Real EPSILON=.0000000001;
 struct Section {unsigned edges;unsigned getNbCollisionEdges(){return edges;}};
 struct Sections {std::array<Section,2>values{{{50},{10}}};size_t size(){return values.size();}Section*get(unsigned i){return &values[i];}}l_sectionMaterials;
 Data<std::vector<Real>>d_keyPoints{{0.,parameters().straight,parameters().total()}};
 std::vector<Element> beams;RestShape rest;
 Data<std::vector<Real>>d_lengthList;Data<std::vector<unsigned>>d_edgeList;
 BezierState cache;BezierState*m_StateNodes=&cache;
 std::vector<unsigned> collisionBeams;
 void clear(){beams.clear();d_lengthList.value.clear();d_edgeList.value.clear();}void clearCollisionOnBeam(){collisionBeams.clear();}
 void addCollisionOnBeam(unsigned b){collisionBeams.push_back(b);}
 Real getRestTotalLength(){return parameters().total();}
 void getSamplingParameters(std::vector<Real>& x,std::vector<unsigned>& density){x={0.,parameters().straight,parameters().total()};density={50,10};}
 void getCollisionSampling(Real&,Real);
 void addBeam(unsigned edge,Real length,Real x0,Real x1,Real angle){Quat q;q.axisToQuat(Vec3(1,0,0),angle);q.normalize();beams.push_back({edge,edge+1,length,x0,x1,Transform(Vec3(),q),Transform(Vec3(),q)});d_lengthList.value.push_back(length);d_edgeList.value.push_back(edge);}
 unsigned getNumBeams()const{return beams.size();} unsigned getStateSize()const{return 61;}
 int getNodeIndices(unsigned id,unsigned& n0,unsigned& n1){n0=beams[id].node0;n1=beams[id].node1;return 1;}
 Real getLength(unsigned id){return beams[id].length;}
 Real getLength(){return parameters().total();}
 void getAbsCurvXFromBeam(unsigned id,Real&a,Real&b){a=beams[id].x0;b=beams[id].x1;}
 const char*getName(){return "WireBeamInterpolation";}
 void getControlPointsFromFrame(const Transform&,const Transform&,const Real,Vec3&,Vec3&,Vec3&,Vec3&);
 void RotateFrameForAlignNormalizedX(const Quat&,const Vec3&,Quat&);
 void InterpolateTransformUsingSpline(Transform&,const Real,const Transform&,const Transform&,const Real);
 void InterpolateTransformUsingSpline(unsigned,const Real,const Vec3&,const VecCoord&,Transform&);
 void interpolatePointUsingSpline(unsigned,const Real,const Vec3&,const VecCoord&,Vec3&,bool,const sofa::core::ConstVecCoordId&);
 void getBeamAtCurvAbs(const Real,unsigned&,Real&,unsigned start=0);
 void updateBezierPoints(const VecCoord& x,int id){auto& points=cache.positions[id].value;points.resize(beams.size()*4);for(unsigned i=0;i<beams.size();i++){Transform a,b;computeTransform(i,a,b,x);getControlPointsFromFrame(a,b,getLength(i),points[4*i],points[4*i+1],points[4*i+2],points[4*i+3]);}}
 void computeActualLength(Real&,const Vec3&,const Vec3&,const Vec3&,const Vec3&);
 void computeStrechAndTwist(unsigned,const VecCoord&,Vec3&,Vec3&);
 void MapForceOnNodeUsingSpline(unsigned,const Real,const Vec3&,const VecCoord&,const Vec3&,SpatialVector&,SpatialVector&);
 void getDOFtoLocalTransform(unsigned id,Transform& a,Transform& b){a=beams[id].local0;b=beams[id].local1;}
 int computeTransform(unsigned id,unsigned n0,unsigned n1,Transform& a,Transform& b,const VecCoord& x){
  a=Transform(x[n0].p,x[n0].q)*beams[id].local0;b=Transform(x[n1].p,x[n1].q)*beams[id].local1;return 1;}
 int computeTransform(unsigned id,Transform& a,Transform& b,const VecCoord& x){return computeTransform(id,beams[id].node0,beams[id].node1,a,b,x);}
 void getDOFtoLocalTransformInGlobalFrame(unsigned id,Transform& a,Transform& b,const VecCoord& x){
  a=Transform(Vec3(),x[beams[id].node0].q)*beams[id].local0;b=Transform(Vec3(),x[beams[id].node1].q)*beams[id].local1;}
 void getSplinePoints(unsigned id,const VecCoord& x,Vec3& a,Vec3& b,Vec3& c,Vec3& d){Transform t0,t1;computeTransform(id,t0,t1,x);getControlPointsFromFrame(t0,t1,getLength(id),a,b,c,d);}
 void getSplineRestTransform(unsigned id,Transform& a,Transform& b){Transform middle,t0,t1;rest.getRestTransformOnX(middle,(beams[id].x0+beams[id].x1)/2);rest.getRestTransformOnX(t0,beams[id].x0);rest.getRestTransformOnX(t1,beams[id].x1);a=middle.inversed()*t0;b=middle.inversed()*t1;}
 void getInterpolationParameters(unsigned id,Real& l,Real& A,Real& iy,Real& iz,Real& sy,Real& sz,Real& J){l=beams[id].length;A=M_PI*parameters().radius*parameters().radius;iy=iz=M_PI*pow(parameters().radius,4)/4;J=iy+iz;sy=sz=0.;}
 void getMechanicalParameters(unsigned,Real& E,Real& nu,Real& rho){E=parameters().young;nu=parameters().poisson;rho=parameters().density;}
};
namespace sofa::core::vec_id::write_access {struct Position{operator int()const{return 0;}};struct RestPosition{};struct Velocity{};inline Position position;inline int freePosition=1;inline RestPosition restPosition;inline Velocity velocity;}
struct MechanicalState {
 DataVecCoord positions,restPositions;DataVecDeriv velocities;
 DataVecCoord* write(sofa::core::vec_id::write_access::Position){return &positions;}
 DataVecCoord* write(sofa::core::vec_id::write_access::RestPosition){return &restPositions;}
 DataVecDeriv* write(sofa::core::vec_id::write_access::Velocity){return &velocities;}
 DataVecCoord*read(sofa::core::vec_id::write_access::Position){return &positions;}
};
template<class D>struct WriteAccessor;
template<class T>struct WriteAccessor<Data<T>>:Access<T>{WriteAccessor(Data<T>& d):Access<T>{&d.value}{}};
struct FixedConstraint {std::vector<unsigned>indices;void clearConstraints(){indices.clear();}void addConstraint(unsigned i){indices.push_back(i);}};
struct BaseMeshTopology {
 using EdgeID=unsigned;using Edge=sofa::type::Vec<2,unsigned>;using Triangle=sofa::type::Vec<3,unsigned>;
 std::vector<Edge>edges;std::vector<Triangle>triangles;
 std::vector<std::vector<unsigned>>vertexEdges,vertexTriangles,edgeTriangles;
 void clear(){edges.clear();triangles.clear();}void addEdge(unsigned a,unsigned b){edges.emplace_back(a,b);}
 const Edge&getEdge(unsigned i){return edges[i];}const Triangle&getTriangle(unsigned i){return triangles[i];}
 const std::vector<unsigned>&getEdgesAroundVertex(unsigned i){return vertexEdges[i];}
 const std::vector<unsigned>&getTrianglesAroundVertex(unsigned i){return vertexTriangles[i];}
 const std::vector<unsigned>&getTrianglesAroundEdge(unsigned i){return edgeTriangles[i];}
 int getRevision(){return 0;}
};
struct Controller {
 using RealConstIterator=std::vector<Real>::const_iterator;
 MechanicalState* state;FixedConstraint* l_fixedConstraint;
 std::vector<Interpolation*>m_instrumentsList;
 Data<Real> d_threshold{parameters().threshold},d_step{.05};Data<unsigned>d_indexFirstNode{60};
 Data<std::vector<Real>> d_rigidCurvAbs,d_xTip{{0.}},d_rotationInstrument{{0.}};
 Data<Pose>d_startingPos{Pose{Vec3(),Quat(0,.7071068,0,.7071068)}};
 std::vector<Real>m_nodeCurvAbs{0.};std::vector<std::vector<int>>m_idInstrumentCurvAbsTable{{0}};
 std::vector<bool>m_activatedPointsBuf;
 MechanicalState*getMechanicalState(){return state;}
 void computeInstrumentsCurvAbs(std::vector<Real>&,const std::vector<Real>&,const Real&);
 void totalLengthIsChanging(const std::vector<Real>&,std::vector<Real>&,const std::vector<std::vector<int>>&);
 void sortCurvAbs(std::vector<Real>&);
 void fillInstrumentCurvAbsTable(const std::vector<Real>&,const std::vector<Real>&,const std::vector<Real>&,std::vector<std::vector<int>>&);
 void fixFirstNodesWithUntil(unsigned);
 void activateBeamListForCollision(std::vector<Real>&,std::vector<std::vector<int>>&);
 void interventionalRadiologyCollisionControls(std::vector<Real>&,std::vector<int>&,std::vector<int>&);
 void applyInterventionalRadiologyController();
};
struct BeamLocalMatrices {double _L,_A,_Iy,_Iz,_Asy,_Asz,_J,_youngM,_cPoisson,_rho;
 Matrix6x6 m_K00,m_K01,m_K10,m_K11,m_M00,m_M01,m_M10,m_M11,m_A0Ref,m_A1Ref;};
struct ForceField {
 Interpolation* l_interpolation; std::vector<BeamLocalMatrices> m_localBeamMatrices;
 void* mstate=nullptr;
 Data<bool> d_computeMass{true},d_reinforceLength{false},d_useShearStressComputation{true};
 Data<std::vector<Real>> d_massDensity{{1.}}; Real m_defaultMassDensity=1.;Vec6 m_gravity;
 void computeGravityVector(){m_gravity.clear();m_gravity[1]=parameters().gravity;}
 void computeStiffness(unsigned,BeamLocalMatrices&);void computeMass(unsigned,BeamLocalMatrices&);
 void applyMassLarge(VecDeriv&,unsigned,unsigned,unsigned,SReal);
 void applyStiffnessLarge(VecDeriv&,const VecDeriv&,unsigned,unsigned,unsigned,SReal);
 void addForce(const MechanicalParams*,DataVecDeriv&,const DataVecCoord&,const DataVecDeriv&);
 void addMDx(const MechanicalParams*,DataVecDeriv&,const DataVecDeriv&,SReal);
 void addMBKToMatrix(const MechanicalParams*,const sofa::core::behavior::MultiMatrixAccessor*);
};
