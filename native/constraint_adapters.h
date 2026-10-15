#pragma once
#include "kernel.h"
#include <memory>
inline constexpr char msgendl='\n';
using Size=sofa::Size;
namespace core=sofa::core;
template<class V>struct SparseMatrix {
 using Row=std::map<unsigned,V>;using Rows=std::map<unsigned,Row>;Rows rows;
 struct ColIterator {typename Row::const_iterator it;unsigned index()const{return it->first;}const V&val()const{return it->second;}void operator++(){++it;}bool operator!=(const ColIterator& b)const{return it!=b.it;}};
 struct RowConstIterator {typename Rows::const_iterator it;unsigned index()const{return it->first;}ColIterator begin()const{return {it->second.begin()};}ColIterator end()const{return {it->second.end()};}void operator++(){++it;}bool operator!=(const RowConstIterator& b)const{return it!=b.it;}};
 using ColConstIterator=ColIterator;
 struct RowIterator {Row*row;void addCol(unsigned i,const V&v){(*row)[i]+=v;}void setCol(unsigned i,const V&v){(*row)[i]=v;}};
 RowIterator writeLine(unsigned i){return {&rows[i]};}
 RowConstIterator begin()const{return {rows.begin()};}RowConstIterator end()const{return {rows.end()};}RowConstIterator readLine(int i)const{return {rows.find(i)};}size_t size()const{return rows.size();}
};
using SparseRows=SparseMatrix<Vec6>;using SparsePoints=SparseMatrix<Vec3>;
struct ScalarSparse {
 using Row=std::map<int,double>;using Rows=std::map<int,Row>;
 using LineConstIterator=Rows::const_iterator;using LElementConstIterator=Row::const_iterator;
 Rows rows;int n=0,m=0;
 void resize(int r,int c){rows.clear();n=r;m=c;}
 int rowSize()const{return n;}int colSize()const{return m;}
 double element(int i,int j)const{auto r=rows.find(i);if(r==rows.end())return 0;auto c=r->second.find(j);return c==r->second.end()?0:c->second;}
 void add(int i,int j,double v){rows[i][j]+=v;}
 LineConstIterator begin()const{return rows.begin();}LineConstIterator end()const{return rows.end();}
};
struct ConstraintMatrix:ScalarMatrix {
 int n=0;std::vector<double>values;
 void resize(int r,int){n=r;values.assign(n*n,0.);}
 int rowSize()const override{return n;}int colSize()const override{return n;}
 void add(int i,int j,double v)override{values[i*n+j]+=v;}
 double element(int i,int j)const override{return values[i*n+j];}
};
namespace sofa::linearalgebra {using BaseMatrix=ScalarMatrix;using BaseVector=Vector;template<class T>using FullVector=Vector;}
namespace linearalgebra=sofa::linearalgebra;
namespace sofa::core {
 struct VecDerivId {static constexpr int V_FIRST_DYNAMIC_INDEX=1;int value;explicit VecDerivId(int v=0):value(v){}};
 namespace vec_id::write_access {inline VecDerivId dx{2},force{1};}
 namespace vec_id::read_access {struct Jacobian{};struct FreePosition{};struct FreeVelocity{};inline Jacobian constraintJacobian;inline FreePosition freePosition;inline FreeVelocity freeVelocity;}
 namespace execparams {inline void*defaultInstance(){return nullptr;}}
}
struct ConstraintState {
 MechanicalState*state;DataVecCoord freePosition;DataVecDeriv freeVelocity,force,dx;Data<SparseRows>jacobian;
 unsigned getSize(){return state->positions.value.size();}
 Data<SparseRows>*read(core::vec_id::read_access::Jacobian){return &jacobian;}
 DataVecCoord*read(core::vec_id::read_access::FreePosition){return &freePosition;}
 DataVecDeriv*read(core::vec_id::read_access::FreeVelocity){return &freeVelocity;}
 DataVecCoord*write(core::vec_id::write_access::Position){return &state->positions;}
 DataVecDeriv*write(core::vec_id::write_access::Velocity){return &state->velocities;}
 DataVecDeriv*write(core::VecDerivId id){return id.value==1?&force:&dx;}
 void vFree(void*,core::VecDerivId){}
};
struct OdeFactors {double getPositionIntegrationFactor(){return parameters().dt;}double getVelocityIntegrationFactor(){return 1.;}};
template<class T>struct Link {T*ptr;T*get(){return ptr;}T*operator->(){return ptr;}};
// BTD buffers are exposed through the same RHS/solution interface used by correction.
struct ConstraintLinearSystem:LinearSystem {
 ConstraintState*state;Vector rhsBuffer,solutionBuffer;
 ConstraintLinearSystem(Matrix* m,ConstraintState* s):state(s){matrix=m;rhs=&rhsBuffer;solution=&solutionBuffer;}
 void setRHS(core::VecDerivId id){auto& data=state->write(id)->value;for(size_t i=0;i<data.size();i++)rhs->blocks[i]=data[i];}
 void setSystemSolution(core::VecDerivId){}
 ScalarMatrix*getSystemBaseMatrix(){return matrix;}Vector*getSystemRHSBaseVector(){return rhs;}Vector*getSystemSolutionBaseVector(){return solution;}
 void dispatchSystemSolution(core::VecDerivId id){auto& data=state->write(id)->value;data.resize(solution->blocks.size());for(size_t i=0;i<data.size();i++)data[i]=solution->blocks[i];}
};
struct CorrectionSolver {
 BTD*solver;ConstraintLinearSystem*system;
 ConstraintLinearSystem*getLinearSystem(){return system;}
 void solveSystem(){solver->solve(*system->matrix,*system->solution,*system->rhs);}
 void init_partial_solve(){solver->l_linearSystem=system;solver->init_partial_solve();}
 void partial_solve(std::list<int>& a,std::list<int>& b,bool update){solver->partial_solve(a,b,update);}
 void addJMInvJt(ScalarMatrix* w,ScalarSparse* j,double f){solver->addJMInvJt(*w,*j,f);}
};
template<class DataTypes>struct Correction {
 using VecCoord=typename DataTypes::VecCoord;using VecDeriv=typename DataTypes::VecDeriv;using Deriv=typename DataTypes::Deriv;
 using MatrixDeriv=SparseRows;using MatrixDerivRowConstIterator=SparseRows::RowConstIterator;using MatrixDerivColConstIterator=SparseRows::ColIterator;
 using ListIndex=std::list<int>;using ComponentState=sofa::core::objectmodel::ComponentState;
 ConstraintState*mstate;CorrectionSolver*l_linearSolver;Link<OdeFactors>l_ODESolver;
 Data<bool>wire_optimization{true};Data<ComponentState>d_componentState{ComponentState::Valid};
 VecDeriv constraint_force;ListIndex constraint_dofs;std::vector<ListIndex>Vec_I_list_dof;
 ScalarSparse m_constraintJacobian;ScalarMatrix*systemMatrix_buf=nullptr;
 Vector*systemRHVector_buf=nullptr,*systemLHVector_buf=nullptr,*systemLHVector_buf_fullvector=nullptr;
 int last_disp=0,last_force=0;bool _new_force=false;
 void verify_constraints(){}void addRegularization(ScalarMatrix*){}const char*getName(){return "ExtractedCorrection";}
 bool hasConstraintNumber(int index){return mstate->jacobian.value.rows.contains(index);}
 void resetForUnbuiltResolution(SReal*,std::list<unsigned>&);
 void addConstraintDisplacement(SReal*,int,int);
 void setConstraintDForce(SReal*,int,int,bool);
 void getBlockDiagonalCompliance(ScalarMatrix*,int,int);
 void applyContactForce(const Vector*);
};
