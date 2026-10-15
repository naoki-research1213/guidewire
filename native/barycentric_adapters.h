#pragma once
#include "mesh_adapters.h"
#include "mapping_adapters.h"
struct PointData {
 using Real=double;using Coord=Vec3;using Deriv=Vec3;using VecCoord=VectorVec3;using MatrixDeriv=SparsePoints;
 static const Vec3&getCPos(const Vec3&p){return p;}static const Vec3&getDPos(const Vec3&p){return p;}static void setCPos(Vec3&p,const Vec3&q){p=q;}
};
struct BaryMapper {
 using In=PointData;using Out=PointData;using Index=unsigned;using OutReal=double;using InDeriv=Vec3;
 using Edge=MeshTopology::Edge;using Triangle=MeshTopology::Triangle;using Quad=MeshTopology::Quad;using Tetra=MeshTopology::Tetra;using Hexa=MeshTopology::Hexa;
 using SeqLines=MeshTopology::SeqLines;using SeqTriangles=MeshTopology::SeqTriangles;using SeqQuads=MeshTopology::SeqQuads;using SeqTetrahedra=MeshTopology::SeqTetrahedra;using SeqHexahedra=MeshTopology::SeqHexahedra;
 struct MappingData1D{unsigned in_index;double baryCoords[1];};struct MappingData2D{unsigned in_index;double baryCoords[2];};struct MappingData3D{unsigned in_index;double baryCoords[3];};
 std::vector<MappingData1D>m_map1d;std::vector<MappingData2D>m_map2d;std::vector<MappingData3D>m_map3d;MeshTopology*m_fromTopology;
 unsigned addPointInLine(unsigned,const SReal*);unsigned addPointInTriangle(unsigned,const SReal*);
 unsigned createPointInLine(const Vec3&,unsigned,const VectorVec3*);unsigned createPointInTriangle(const Vec3&,unsigned,const VectorVec3*);
 void apply(VectorVec3&,const VectorVec3&);void applyJT(SparsePoints&,const SparsePoints&);
};
