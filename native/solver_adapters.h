#pragma once
#include "adapters.h"
#include <map>
#include <list>
#include <sstream>
#include <stdexcept>
namespace sofa {using SignedIndex=int;}
using Index=int;using SignedIndex=int;
#define msg_info_when(x) Quiet()
#define msg_error_when(x) Quiet()
#define dmsg_info_when(x) Quiet()
#define dmsg_warning() Quiet()
using std::make_pair;
// Original BTDMatrix::Block is inserted by the extraction script below.
#include "generated/btd_block.inc"
#include "generated/btd_vector.inc"
struct Matrix:ScalarMatrix {
 std::vector<Block>blocks;int n=0;
 static constexpr int getSubMatrixDim(){return 6;}
 void resize(int rows,int){n=rows/6;blocks.assign(n*n,Block());}
 int rowSize()const override{return n*6;}int colSize()const override{return n*6;}
 Block& asub(int i,int j,int,int){return blocks[i*n+j];}
 void getAlignedSubMatrix(int i,int j,int,int,Block& b){b=blocks[i*n+j];}
 void setAlignedSubMatrix(int i,int j,int,int,const Block& b){blocks[i*n+j]=b;}
 double element(int i,int j)const override{return blocks[(i/6)*n+j/6][i%6][j%6];}
 void add(int i,int j,double v)override{blocks[(i/6)*n+j/6][i%6][j%6]+=v;}
};
struct Vector {
 std::vector<SubVector>blocks;
 explicit Vector(int n=61):blocks(n){}
 int size()const{return blocks.size()*6;}
 SubVector& asub(int i,int){return blocks[i];}
 const SubVector& asub(int i,int)const{return blocks[i];}
 void clear(){blocks.clear();}void resize(int size){blocks.resize(size/6);}
 Vector&operator-=(const Vector& x){for(size_t i=0;i<blocks.size();i++)blocks[i]-=x.blocks[i];return *this;}
 double element(int i)const{return blocks[i/6][i%6];}void set(int i,double v){blocks[i/6][i%6]=v;}
 double&operator[](int i){return blocks[i/6][i%6];}
};
struct LinearSystem {Matrix*matrix;Vector*rhs,*solution;Matrix*getSystemMatrix(){return matrix;}Vector*getRHSVector(){return rhs;}Vector*getSolutionVector(){return solution;}};
struct BTD {
 using SubMatrix=Block;using BlocType=Matrix6x6;using IndexPair=std::pair<int,int>;
 Data<bool>d_verbose{false},d_subpartSolve{false};
 Data<bool>d_problem{false},d_verification{false};
 using ListIndex=std::list<Index>;
 LinearSystem* l_linearSystem=nullptr;
 Vector bwdContributionOnLH,fwdContributionOnRH,_rh_buf;SubVector _acc_rh_bloc,_acc_lh_bloc;
 Index current_bloc=0,first_block=0,_indMaxNonNullForce=0,_indMaxFwdLHComputed=0;
 std::vector<SubVector>Vec_dRH;
 std::vector<Block>alpha_inv,lambda,B;std::vector<int>nBlockComputedMinv;
 Matrix Minv;std::map<IndexPair,Block>H;std::map<IndexPair,Block>::iterator H_it;
 void init_partial_inverse(const Index&,const Index&);
 void init_partial_solve();void bwdAccumulateRHinBloc(Index);void bwdAccumulateLHGlobal();
 void fwdAccumulateRHGlobal(Index);void fwdComputeLHinBloc(Index);void partial_solve(ListIndex&,ListIndex&,bool);
 void my_identity(SubMatrix&,const Index);
 void invert(SubMatrix&,const BlocType&);
 void invert(Matrix&);
 void solve(Matrix&,Vector&,Vector&);
 void computeMinvBlock(Index,Index);
 double getMinvElement(Index,Index);
 bool notMuted(){return false;}
 template<class RMatrix,class JMatrix>bool addJMInvJt(RMatrix&,JMatrix&,double);
};
