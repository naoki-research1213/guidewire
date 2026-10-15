#pragma once
#include "implementation.h"
// interactive_xray_scene.py の現行構成を固定配列・リンクへ展開する。
// 計算本体は generated/ の原関数。ここは状態の所有と呼び出し順の接続。
// FreeMotionAnimationLoop の自由運動→写像→検出→拘束求解→補正に対応する。
struct Scene {
 Parameters settings=parameters();
 Kernel kernel;MeshTopology vesselTopology,wireTopology;CollisionState vesselState,wireState;
 VectorVec3 wireFree;sofa::simulation::Node vesselNode,wireNode;
 // 原 CollisionResponse はモデル対のポインタ順に反復する。現行で観測した
 // 三角形/点→辺/点→辺/辺の順序を保持する（履歴・反復解に影響する）。
 TriangleModel vesselTriangles;LineModel vesselLines;PointModel wirePoints;LineModel wireLines;
 LocalMinDistance lmd;IntersectionMethod intersection{&lmd};BVHNarrowPhase narrow;
 MappingOutput collisionOutput,visualOutput;Mapping collisionMapping{&kernel.interpolation},visualMapping{&kernel.interpolation};MappingTopology visualTopology;
 std::vector<double>xPoints=std::vector<double>(61);std::vector<int>instrumentIds=std::vector<int>(61);
 ConstraintState constraintState{&kernel.state};ConstraintLinearSystem constraintSystem{&kernel.system,&constraintState};CorrectionSolver linear{&kernel.solver,&constraintSystem};OdeFactors ode;
 Correction<RigidData>correction{&constraintState,&linear,{&ode}};ContactSolver contactSolver;
 std::map<std::pair<CollisionModel*,CollisionModel*>,std::unique_ptr<ContactResponse>>responses;
 std::vector<double>forces,residuals;unsigned contactCount=0;int converged=0;
 SparsePoints collisionJacobian;
 static CollisionModel*first(CollisionModel*m){while(m->getPrevious())m=m->getPrevious();return m;}
 Scene(const char*filename){
  MeshSTLLoader loader;if(!loader.readBinarySTL(filename)){std::ifstream file(filename);if(!file||!loader.readSTL(file))throw std::runtime_error("STL loading failed");}
  vesselState.positions.value=loader.d_positions.value;vesselTopology.triangles=loader.d_triangles.value;vesselTopology.nbPoints=vesselState.getSize();vesselTopology.connect();
  vesselNode.line=&vesselLines;wireNode.line=&wireLines;
  for(CollisionModel*m:{static_cast<CollisionModel*>(&vesselTriangles),static_cast<CollisionModel*>(&vesselLines)}){m->state=&vesselState;m->topology=&vesselTopology;m->node=&vesselNode;m->moving=false;m->simulated=false;}
  vesselTriangles.m_mstate=&vesselState;vesselTriangles.m_topology=&vesselTopology;vesselTriangles.updateFromTopology();vesselLines.mstate=&vesselState;
  wireTopology.nbPoints=61;for(unsigned i=0;i<60;i++)wireTopology.d_seqEdges.value.emplace_back(i,i+1);wireTopology.connect();wireState.positions.value.resize(61);
  for(CollisionModel*m:{static_cast<CollisionModel*>(&wireLines),static_cast<CollisionModel*>(&wirePoints)}){m->state=&wireState;m->topology=&wireTopology;m->node=&wireNode;}
  wireLines.mstate=&wireState;wirePoints.mstate=&wireState;
  collisionMapping.fromModel=&kernel.state;collisionMapping.toModel=&collisionOutput;collisionMapping.m_isSubMapping=true;collisionMapping.d_points.value.resize(61);collisionMapping.m_idPointSubMap.resize(61);
  for(unsigned i=0;i<61;i++)collisionMapping.m_idPointSubMap[i]=i;
  visualMapping.fromModel=&kernel.state;visualMapping.toModel=&visualOutput;visualMapping.context.topology=&visualTopology;visualMapping.d_nbPointsPerBeam.value=8;
  narrow.intersectionMethod=&intersection;contactSolver.l_constraintCorrections={&correction};
  vesselTriangles.computeBoundingTree(6);vesselLines.computeBoundingTree(6);
  updateCollisionPoints();updateDisplay();
 }
 void updateCollisionPoints(){
  std::vector<int>remove;kernel.controller.interventionalRadiologyCollisionControls(xPoints,instrumentIds,remove);
  if(!remove.empty())throw std::logic_error("dropped instrument topology not extracted");
  for(unsigned i=0;i<61;i++)collisionMapping.d_points.value[i]=Vec3(xPoints[i],0,0);
  Data<VectorVec3>mapped{wireState.positions.value};MappingParams params;collisionMapping.apply(&params,mapped,kernel.state.positions);wireState.positions.value=std::move(mapped.value);
 }
 void updateDisplay(){Data<VectorVec3>out;MappingParams params;visualMapping.apply(&params,out,kernel.state.positions);visualOutput.positions=std::move(out.value);}
 void advance(double insert,double rotate){
  kernel.freeMotion(insert,rotate);updateCollisionPoints();
  // 原 FreeMotionAnimationLoop は自由速度を J で写像した後、
  // 接触点の自由位置を現在の写像位置 + dt * 写像された自由速度で更新する。
  Data<VecDeriv>v{kernel.freeVelocities};Data<VectorVec3>mappedVelocity{VectorVec3(61)};MappingParams params;collisionMapping.applyJ(&params,mappedVelocity,v);
  wireFree=wireState.positions.value;for(unsigned i=0;i<wireFree.size();i++)wireFree[i]+=mappedVelocity.value[i]*kernel.h;
  wireLines.computeBoundingTree(6);wirePoints.computeBoundingTree(6);
  for(auto&[key,out]:narrow.outputs)if(out)out->clear();
  narrow.addCollisionPair({first(&vesselTriangles),first(&wirePoints)});
  narrow.addCollisionPair({first(&vesselLines),first(&wirePoints)});
  narrow.addCollisionPair({first(&vesselLines),first(&wireLines)});
  constraintState.freePosition.value=kernel.freePositions;constraintState.freeVelocity.value=kernel.freeVelocities;constraintState.jacobian.value.rows.clear();
  contactSolver.m_constraintBlockInfo.clear();contactSolver.m_constraintIds.clear();
  collisionJacobian.rows.clear();
  unsigned row=0;Vector violation(1);std::vector<ContactResponse*>active;
  for(auto&[key,out]:narrow.outputs){
   if(!out||out->empty()){responses.erase(key);continue;}
   auto&owned=responses[key];if(!owned){owned=std::make_unique<ContactResponse>();owned->model1=key.first;owned->model2=key.second;owned->intersectionMethod=&intersection;
    owned->mapper1.model=key.first;owned->mapper2.model=key.second;owned->mapper1.bary.m_fromTopology=&vesselTopology;owned->mapper2.bary.m_fromTopology=&wireTopology;
    owned->mapper1.freeInput=&vesselState.positions.value;owned->mapper2.freeInput=&wireFree;owned->setupConstraint(&owned->mapper1.output,&owned->mapper2.output);
   }
   auto&r=*owned;r.setDetectionOutputs(out);r.createResponse(nullptr);
   Data<SparsePoints>j1,j2;r.constraint.buildConstraintMatrix(nullptr,j1,j2,row,r.mapper1.output.positions,r.mapper2.output.positions);
   Data<SparsePoints>wireJ;r.mapper2.applyJT(wireJ.value,j2.value);for(auto&[r,cols]:wireJ.value.rows)for(auto&[p,v]:cols)collisionJacobian.rows[r][p]+=v;collisionMapping.applyJT(nullptr,constraintState.jacobian,wireJ);
   r.constraint.getConstraintInfo(nullptr,contactSolver.m_constraintBlockInfo,contactSolver.m_constraintIds);active.push_back(&r);
  }
  contactCount=row/3;forces.assign(row,0.);residuals.clear();converged=1;
  if(row){
   violation.blocks.resize((row+5)/6);for(auto*r:active)r->constraint.getPositionViolation(&violation);
   std::vector<double>d(row);for(unsigned i=0;i<row;i++)d[i]=violation.element(i);
   contactSolver._numConstraints=row;contactSolver._Wdiag.resize(row,row);contactSolver._result=&forces;contactSolver.computeInitialGuess();
   converged=contactSolver.nlcp_gaussseidel_unbuilt(d.data(),forces.data(),&residuals);contactSolver.keepContactForcesValue();
   Vector lambda((row+5)/6);for(unsigned i=0;i<row;i++)lambda.set(i,forces[i]);correction.applyContactForce(&lambda);
  }else{kernel.state.positions.value=kernel.freePositions;kernel.state.velocities.value=kernel.freeVelocities;}
  // 拘束補正後に固定節点を射影し、現在の接触点・表示中心線へ写像する。
  for(unsigned i:kernel.fixed.indices){kernel.state.positions.value[i]=kernel.state.restPositions.value[i];kernel.state.velocities.value[i].clear();}
  kernel.time+=kernel.h;Data<VectorVec3>mapped{wireState.positions.value};collisionMapping.apply(&params,mapped,kernel.state.positions);wireState.positions.value=std::move(mapped.value);updateDisplay();
 }
};
