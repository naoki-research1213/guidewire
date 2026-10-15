#pragma once
#include "bvh_adapters.h"
#include <fstream>
using sofa::helper::getWriteOnlyAccessor;
using std::vector;
namespace topology {using Triangle=BaseMeshTopology::Triangle;}
struct MeshSTLLoader {
 using Triangle=topology::Triangle;
 Data<VectorVec3>d_positions,d_normals;Data<std::vector<Triangle>>d_triangles;
 Data<bool>d_mergePositionUsingMap{true},d_flipNormals{false};Data<int>d_headerSize{80};
 bool readBinarySTL(const char*);bool readSTL(std::ifstream&);void addTriangle(std::vector<Triangle>&,const Triangle&);
};
struct MeshTopology:BaseMeshTopology {
 using SeqEdges=std::vector<Edge>;using SeqTriangles=std::vector<Triangle>;using SeqQuads=std::vector<sofa::type::Vec<4,unsigned>>;
 using Quad=sofa::type::Vec<4,unsigned>;using Tetra=sofa::type::Vec<4,unsigned>;using Hexa=sofa::type::Vec<8,unsigned>;using SeqLines=SeqEdges;using SeqTetrahedra=std::vector<Tetra>;using SeqHexahedra=std::vector<Hexa>;
 SeqTetrahedra tetrahedra;SeqHexahedra hexas;
 using EdgesInTriangle=Triangle;static constexpr unsigned InvalidID=std::numeric_limits<unsigned>::max();
 Data<SeqEdges>d_seqEdges;SeqQuads quads;unsigned nbPoints=0;
 std::vector<EdgesInTriangle>m_edgesInTriangle;
 std::vector<std::vector<unsigned>>&m_edgesAroundVertex=vertexEdges;
 std::vector<std::vector<unsigned>>&m_trianglesAroundVertex=vertexTriangles;
 std::vector<std::vector<unsigned>>&m_trianglesAroundEdge=edgeTriangles;
 std::map<Edge,unsigned>edgeIds;
 const SeqEdges&getEdges(){return d_seqEdges.value;}const SeqTriangles&getTriangles(){return triangles;}const SeqQuads&getQuads(){return quads;}
 unsigned getNbEdges(){return d_seqEdges.value.size();}
 const SeqLines&getLines(){return d_seqEdges.value;}const Edge&getLine(unsigned i){return d_seqEdges.value[i];}const SeqTetrahedra&getTetrahedra(){return tetrahedra;}const SeqHexahedra&getHexahedra(){return hexas;}
 unsigned getEdgeIndex(unsigned a,unsigned b){auto it=edgeIds.find(Edge(std::min(a,b),std::max(a,b)));return it==edgeIds.end()?InvalidID:it->second;}
 struct EdgeUpdate {MeshTopology*topology;using Edge=MeshTopology::Edge;void updateFromSurface();};
 void createEdgesAroundVertexArray();void createEdgesInTriangleArray();void createTrianglesAroundVertexArray();void createTrianglesAroundEdgeArray();
 void connect(){
  if(d_seqEdges.value.empty()&&!triangles.empty()){EdgeUpdate update{this};update.updateFromSurface();}
  edges=d_seqEdges.value;edgeIds.clear();for(unsigned i=0;i<edges.size();i++)edgeIds[Edge(std::min(edges[i][0],edges[i][1]),std::max(edges[i][0],edges[i][1]))]=i;
  createEdgesAroundVertexArray();createTrianglesAroundVertexArray();createTrianglesAroundEdgeArray();
 }
};
