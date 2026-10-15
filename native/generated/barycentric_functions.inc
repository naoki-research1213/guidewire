/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#pragma once

// BarycentricMapperMeshTopology.inl:268
unsigned 
BaryMapper::addPointInLine ( const Index lineIndex, const SReal* baryCoords )
{
    m_map1d.resize ( m_map1d.size() +1 );
    MappingData1D& data = *m_map1d.rbegin();
    data.in_index = lineIndex;
    data.baryCoords[0] = ( Real ) baryCoords[0];
    return Index(m_map1d.size()-1);
}

// BarycentricMapperMeshTopology.inl:280
unsigned 
BaryMapper::addPointInTriangle ( const Index triangleIndex, const SReal* baryCoords )
{
    m_map2d.resize ( m_map2d.size() +1 );
    MappingData2D& data = *m_map2d.rbegin();
    data.in_index = triangleIndex;
    data.baryCoords[0] = ( Real ) baryCoords[0];
    data.baryCoords[1] = ( Real ) baryCoords[1];
    return Index(m_map2d.size()-1);
}

// BarycentricMapperMeshTopology.inl:334
unsigned 
BaryMapper::createPointInLine ( const typename Out::Coord& p, Index lineIndex, const typename In::VecCoord* points )
{
    SReal baryCoords[1];
    const Edge& elem = this->m_fromTopology->getLine ( lineIndex );
    const typename In::Coord p0 = ( *points ) [elem[0]];
    const typename In::Coord pA = ( *points ) [elem[1]] - p0;
    typename In::Coord pos = Out::getCPos(p) - p0;

    const SReal L2 = pA.norm2(); 
    if (L2 < std::numeric_limits<SReal>::epsilon()) // in case of null length edge, avoid division by 0
        baryCoords[0] = 0.0;
    else
        baryCoords[0] = ((pos * pA) / L2);

    return this->addPointInLine ( lineIndex, baryCoords );
}

// BarycentricMapperMeshTopology.inl:354
unsigned 
BaryMapper::createPointInTriangle ( const typename Out::Coord& p, Index triangleIndex, const typename In::VecCoord* points )
{
    SReal baryCoords[2];
    const Triangle& elem = this->m_fromTopology->getTriangle ( triangleIndex );

    const typename In::Coord & p1 = ( *points ) [elem[0]];
    const typename In::Coord & p2 = ( *points ) [elem[1]];
    const typename In::Coord & p3 = ( *points ) [elem[2]];
    const typename In::Coord & to_be_projected = Out::getCPos(p);

    const typename In::Coord AB = p2-p1;
    const typename In::Coord AC = p3-p1;
    const typename In::Coord AQ = to_be_projected -p1;
    sofa::type::Mat<2,2,typename In::Real> A;
    sofa::type::Vec<2,typename In::Real> b;
    A(0,0) = AB*AB;
    A(1,1) = AC*AC;
    A(0,1) = A(1,0) = AB*AC;
    b[0] = AQ*AB;
    b[1] = AQ*AC;
    const typename In::Real det = sofa::type::determinant(A);

    baryCoords[0] = (b[0]*A(1,1) - b[1]*A(0,1))/det;
    baryCoords[1]  = (b[1]*A(0,0) - b[0]*A(1,0))/det;

    if (baryCoords[0] < 0 || baryCoords[1] < 0 || baryCoords[0] + baryCoords[1] > 1)
    {
        // nearest point is on an edge or corner
        // barycentric coordinate on AB
        const SReal pAB = b[0] / A(0,0); // AQ*AB / AB*AB
        // barycentric coordinate on AC
        const SReal pAC = b[1] / A(1,1); // AQ*AC / AB*AB
        if (pAB < 0 && pAC < 0)
        {
            // closest point is A
            baryCoords[0] = 0.0;
            baryCoords[1] = 0.0;
        }
        else if (pAB < 1 && baryCoords[1] < 0)
        {
            // closest point is on AB
            baryCoords[0] = pAB;
            baryCoords[1] = 0.0;
        }
        else if (pAC < 1 && baryCoords[0] < 0)
        {
            // closest point is on AC
            baryCoords[0] = 0.0;
            baryCoords[1] = pAC;
        }
        else
        {
            // barycentric coordinate on BC
            // BQ*BC / BC*BC = (AQ-AB)*(AC-AB) / (AC-AB)*(AC-AB) = (AQ*AC-AQ*AB + AB*AB-AB*AC) / (AB*AB+AC*AC-2AB*AC)
            const SReal pBC = (b[1] - b[0] + A(0,0) - A(0,1)) / (A(0,0) + A(1,1) - 2*A(0,1)); // BQ*BC / BC*BC
            if (pBC < 0)
            {
                // closest point is B
                baryCoords[0] = 1.0;
                baryCoords[1] = 0.0;
            }
            else if (pBC > 1)
            {
                // closest point is C
                baryCoords[0] = 0.0;
                baryCoords[1] = 1.0;
            }
            else
            {
                // closest point is on BC
                baryCoords[0] = 1.0-pBC;
                baryCoords[1] = pBC;
            }
        }
    }

    return this->addPointInTriangle ( triangleIndex, baryCoords );
}

// BarycentricMapperMeshTopology.inl:1003
void BaryMapper::apply ( typename Out::VecCoord& out, const typename In::VecCoord& in )
{
    out.resize( m_map1d.size() +m_map2d.size() +m_map3d.size() );

    const SeqLines& lines = this->m_fromTopology->getLines();
    const SeqTriangles& triangles = this->m_fromTopology->getTriangles();
    const SeqQuads& quads = this->m_fromTopology->getQuads();
    const SeqTetrahedra& tetrahedra = this->m_fromTopology->getTetrahedra();
    const SeqHexahedra& cubes = this->m_fromTopology->getHexahedra();

    // 1D elements
    {
        for ( std::size_t i=0; i<m_map1d.size(); i++ )
        {
            const Real fx = m_map1d[i].baryCoords[0];
            const Index index = m_map1d[i].in_index;
            {
                const Edge& line = lines[index];
                Out::setCPos(out[i] , in[line[0]] * ( 1-fx )
                        + in[line[1]] * fx );
            }
        }
    }
    // 2D elements
    {
        const std::size_t i0 = m_map1d.size();
        const std::size_t c0 = triangles.size();
        for ( std::size_t i=0; i<m_map2d.size(); i++ )
        {
            const Real fx = m_map2d[i].baryCoords[0];
            const Real fy = m_map2d[i].baryCoords[1];
            const Index index = m_map2d[i].in_index;
            if ( index<c0 )
            {
                const Triangle& triangle = triangles[index];
                Out::setCPos(out[i+i0] , in[triangle[0]] * ( 1-fx-fy )
                        + in[triangle[1]] * fx
                        + in[triangle[2]] * fy );
            }
            else
            {
                if (quads.size())
                {
                    const Quad& quad = quads[index-c0];
                    Out::setCPos(out[i+i0] , in[quad[0]] * ( ( 1-fx ) * ( 1-fy ) )
                            + in[quad[1]] * ( ( fx ) * ( 1-fy ) )
                            + in[quad[3]] * ( ( 1-fx ) * ( fy ) )
                            + in[quad[2]] * ( ( fx ) * ( fy ) ) );
                }
            }
        }
    }
    // 3D elements
    {
        const std::size_t i0 = m_map1d.size() + m_map2d.size();
        const std::size_t c0 = tetrahedra.size();
        for ( std::size_t i=0; i<m_map3d.size(); i++ )
        {
            const Real fx = m_map3d[i].baryCoords[0];
            const Real fy = m_map3d[i].baryCoords[1];
            const Real fz = m_map3d[i].baryCoords[2];
            const Index index = m_map3d[i].in_index;
            if ( index<c0 )
            {
                const Tetra& tetra = tetrahedra[index];
                Out::setCPos(out[i+i0] , in[tetra[0]] * ( 1-fx-fy-fz )
                        + in[tetra[1]] * fx
                        + in[tetra[2]] * fy
                        + in[tetra[3]] * fz );
            }
            else
            {
                const Hexa& cube = cubes[index-c0];

                Out::setCPos(out[i+i0] , in[cube[0]] * ( ( 1-fx ) * ( 1-fy ) * ( 1-fz ) )
                        + in[cube[1]] * ( ( fx ) * ( 1-fy ) * ( 1-fz ) )
                        + in[cube[3]] * ( ( 1-fx ) * ( fy ) * ( 1-fz ) )
                        + in[cube[2]] * ( ( fx ) * ( fy ) * ( 1-fz ) )
                        + in[cube[4]] * ( ( 1-fx ) * ( 1-fy ) * ( fz ) )
                        + in[cube[5]] * ( ( fx ) * ( 1-fy ) * ( fz ) )
                        + in[cube[7]] * ( ( 1-fx ) * ( fy ) * ( fz ) )
                        + in[cube[6]] * ( ( fx ) * ( fy ) * ( fz ) ) );
            }
        }
    }
}

// BarycentricMapperMeshTopology.inl:462
void BaryMapper::applyJT ( typename In::MatrixDeriv& out, const typename Out::MatrixDeriv& in )
{
    const SeqLines& lines = this->m_fromTopology->getLines();
    const SeqTriangles& triangles = this->m_fromTopology->getTriangles();
    const SeqQuads& quads = this->m_fromTopology->getQuads();
    const SeqTetrahedra& tetrahedra = this->m_fromTopology->getTetrahedra();
    const SeqHexahedra& hexas = this->m_fromTopology->getHexahedra();

    const size_t iTri = triangles.size();
    const size_t iTetra= tetrahedra.size();

    const size_t i1d = m_map1d.size();
    const size_t i2d = m_map2d.size();
    const size_t i3d = m_map3d.size();

    size_t indexIn;

    typename Out::MatrixDeriv::RowConstIterator rowItEnd = in.end();

    for (typename Out::MatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        typename Out::MatrixDeriv::ColConstIterator colItEnd = rowIt.end();
        typename Out::MatrixDeriv::ColConstIterator colIt = rowIt.begin();

        if (colIt != colItEnd)
        {
            typename In::MatrixDeriv::RowIterator o = out.writeLine(rowIt.index());

            for ( ; colIt != colItEnd; ++colIt)
            {
                indexIn = colIt.index();
                InDeriv data = (InDeriv) Out::getDPos(colIt.val());

                // 1D elements
                if ( indexIn < i1d )
                {
                    const OutReal fx = ( OutReal ) m_map1d[indexIn].baryCoords[0];
                    const Index index = m_map1d[indexIn].in_index;
                    {
                        const Edge& line = lines[index];
                        o.addCol( line[0], data * ( 1-fx ) );
                        o.addCol( line[1], data * fx );
                    }
                }
                // 2D elements : triangle or quad
                else if ( indexIn < i2d )
                {
                    const OutReal fx = ( OutReal ) m_map2d[indexIn].baryCoords[0];
                    const OutReal fy = ( OutReal ) m_map2d[indexIn].baryCoords[1];
                    const Index index = m_map2d[indexIn].in_index;
                    if ( index < iTri ) // triangle
                    {
                        const Triangle& triangle = triangles[index];
                        o.addCol( triangle[0], data * ( 1-fx-fy ) );
                        o.addCol( triangle[1], data * fx );
                        o.addCol( triangle[2], data * fy );
                    }
                    else // quad
                    {
                        const Quad& quad = quads[index - iTri];
                        o.addCol( quad[0], data * ( ( 1-fx ) * ( 1-fy ) ) );
                        o.addCol( quad[1], data * ( ( fx ) * ( 1-fy ) ) );
                        o.addCol( quad[3], data * ( ( 1-fx ) * ( fy ) ) );
                        o.addCol( quad[2], data * ( ( fx ) * ( fy ) ) );
                    }
                }
                // 3D elements : tetra or hexa
                else if ( indexIn < i3d )
                {
                    const OutReal fx = ( OutReal ) m_map3d[indexIn].baryCoords[0];
                    const OutReal fy = ( OutReal ) m_map3d[indexIn].baryCoords[1];
                    const OutReal fz = ( OutReal ) m_map3d[indexIn].baryCoords[2];
                    const size_t index = m_map3d[indexIn].in_index;
                    if ( index < iTetra ) // tetra
                    {
                        const Tetra& tetra = tetrahedra[index];
                        o.addCol ( tetra[0], data * ( 1-fx-fy-fz ) );
                        o.addCol ( tetra[1], data * fx );
                        o.addCol ( tetra[2], data * fy );
                        o.addCol ( tetra[3], data * fz );
                    }
                    else // hexa
                    {
                        const Hexa& hexa = hexas[index-iTetra];

                        o.addCol ( hexa[0],data * ( ( 1-fx ) * ( 1-fy ) * ( 1-fz ) ) ) ;
                        o.addCol ( hexa[1],data * ( ( fx ) * ( 1-fy ) * ( 1-fz ) ) ) ;

                        o.addCol ( hexa[3],data * ( ( 1-fx ) * ( fy ) * ( 1-fz ) ) ) ;
                        o.addCol ( hexa[2],data * ( ( fx ) * ( fy ) * ( 1-fz ) ) ) ;

                        o.addCol ( hexa[4],data * ( ( 1-fx ) * ( 1-fy ) * ( fz ) ) ) ;
                        o.addCol ( hexa[5],data * ( ( fx ) * ( 1-fy ) * ( fz ) ) ) ;

                        o.addCol ( hexa[7],data * ( ( 1-fx ) * ( fy ) * ( fz ) ) ) ;
                        o.addCol ( hexa[6],data * ( ( fx ) * ( fy ) * ( fz ) ) ) ;
                    }
                }
            }
        }
    }
}
