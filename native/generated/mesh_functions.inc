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

// upstream MeshSTLLoader.cpp:93; 7c18e95d5c5f2839079892c69e7d89a313c79603
bool isBinarySTLValid(const char* filename, const MeshSTLLoader* _this)
{
    // Binary STL files have 80-bytes headers. The following 4-bytes is the number of triangular d_facets in the file
    // Each facet is described with a 50-bytes field, so a valid binary STL file verifies the following condition:
    // nFacets * 50 + 84-bytes header == filesize

    long filesize;
    std::ifstream f(filename, std::ifstream::ate | std::ifstream::binary);
    filesize = f.tellg();
    if (filesize < 84)
    {
        msg_error(_this) << "Can't read binary STL file: " << filename;
        return false;
    }
    f.seekg(0);
    char buffer[80];
    f.read(buffer, 80);
    uint32_t ntriangles;
    f.read(reinterpret_cast<char*>(&ntriangles), 4);
    const uint32_t expectedFileSize = ntriangles * 50 + 84;
    if (filesize != expectedFileSize)
    {
        msg_error(_this) << filename << " isn't binary STL file. File size expected to be "
            << expectedFileSize << " (with " << ntriangles << " triangles) but it is " << filesize;
        return false;
    }
    return true;
}

// upstream MeshSTLLoader.cpp:122; 7c18e95d5c5f2839079892c69e7d89a313c79603
bool MeshSTLLoader::readBinarySTL(const char *filename)
{
    dmsg_info() << "Reading binary STL file..." ;
    if (!isBinarySTLValid(filename, this))
        return false;

    auto my_positions = getWriteOnlyAccessor(d_positions);
    auto my_normals = getWriteOnlyAccessor(d_normals);
    auto my_triangles = getWriteOnlyAccessor(this->d_triangles);

    std::map< sofa::type::Vec3, sofa::Index > my_map;
    sofa::Index positionCounter = 0;
    const bool useMap = d_mergePositionUsingMap.getValue();

    std::ifstream dataFile(filename, std::ios::in | std::ifstream::binary);

    // Skipping header file
    char buffer[256];
    dataFile.read(buffer, d_headerSize.getValue());

    uint32_t nbrFacet;
    dataFile.read(reinterpret_cast<char*>(&nbrFacet), 4);

    my_normals.resize( nbrFacet ); // exact size
    my_positions.reserve( nbrFacet * 3 ); // max size

#ifndef NDEBUG
    {
    // checking that the file is large enough to contain the given nb of d_facets
    // store current pos in file
    std::streampos pos = dataFile.tellg();
    // get length of file
    dataFile.seekg(0, std::ios::end);
    std::streampos length = dataFile.tellg();
    // restore pos in file
    dataFile.seekg(pos);
    // check for length
    assert(length >= d_headerSize.getValue() + 4 + nbrFacet * (12 /*normal*/ + 3 * 12 /*points*/ + 2 /*attribute*/ ) );
    }
#endif

    // temporaries
    sofa::type::Vec3f vertexf, normalf;

    // reserve vector before filling it
    my_triangles.reserve( nbrFacet );

    unsigned int nbDegeneratedTriangles = 0;

    // Parsing d_facets
    for (uint32_t i = 0; i<nbrFacet; ++i)
    {
        topology::Triangle the_tri;

        // Normal:
        dataFile.read((char*)&normalf[0], 4);
        dataFile.read((char*)&normalf[1], 4);
        dataFile.read((char*)&normalf[2], 4);
        my_normals[i] = type::toVec3(normalf);

        // Vertices:
        for (size_t j = 0; j<3; ++j)
        {
            dataFile.read((char*)&vertexf[0], 4);
            dataFile.read((char*)&vertexf[1], 4);
            dataFile.read((char*)&vertexf[2], 4);
            const auto vertex = type::toVec3(vertexf);


            if( useMap )
            {
                auto it = my_map.find( vertex );
                if( it == my_map.end() )
                {
                    the_tri[j] = positionCounter;
                    my_map[vertex] = positionCounter++;
                    my_positions.push_back(vertex);
                }
                else
                {
                    the_tri[j] = it->second;
                }
            }
            else
            {
                bool find = false;
                for (size_t k=0; k<my_positions.size(); ++k)
                    if ( (vertex[0] == my_positions[k][0]) && (vertex[1] == my_positions[k][1])  && (vertex[2] == my_positions[k][2]))
                    {
                        find = true;
                        the_tri[j] = static_cast<sofa::Index>(k);
                        break;
                    }

                if (!find)
                {
                    my_positions.push_back(vertex);
                    the_tri[j] = my_positions.size()-1;
                }
            }
        }

        if (the_tri[0] == the_tri[1] || the_tri[1] == the_tri[2] || the_tri[0] == the_tri[2])
        {
            ++nbDegeneratedTriangles;
        }
        else
        {
            this->addTriangle(my_triangles.wref(), the_tri);
        }

        // Attribute byte count
        uint16_t count;
        dataFile.read((char*)&count, 2);
    }

    if(my_triangles.size() != (size_t)(nbrFacet - nbDegeneratedTriangles))
    {
        msg_error() << "Size mismatch between triangle vector and facetSize";
        return false;
    }

    msg_warning_when(nbDegeneratedTriangles > 0) << "Found " << nbDegeneratedTriangles << " degenerated triangles ("
        "triangles which indices are not all different). Those triangles have not been added to the list of triangles";

    dmsg_info() << "done!" ;
    return true;
}

// upstream MeshSTLLoader.cpp:252; 7c18e95d5c5f2839079892c69e7d89a313c79603
bool MeshSTLLoader::readSTL(std::ifstream& dataFile)
{
    Vec3 result;
    std::string line;

    auto my_positions = getWriteOnlyAccessor(d_positions);
    auto my_normals = getWriteOnlyAccessor(d_normals);
    auto my_triangles = getWriteOnlyAccessor(d_triangles);

    std::map< sofa::type::Vec3, sofa::Index > my_map;
    sofa::Index positionCounter = 0, vertexCounter = 0;
    const bool useMap = d_mergePositionUsingMap.getValue();

    topology::Triangle the_tri;

    while (std::getline(dataFile, line))
    {
        if (line.empty()) continue;
        std::istringstream values(line);

        std::string bufferWord;
        values >> bufferWord;

        if (bufferWord == "facet")
        {
            // Normal
            values >> bufferWord >> result[0] >> result[1] >> result[2];
            my_normals.push_back(result);
        }
        else if (bufferWord == "vertex")
        {
            // Vertex
            values >> result[0] >> result[1] >> result[2];

            if( useMap )
            {
                auto it = my_map.find(result);
                if( it == my_map.end() )
                {
                    the_tri[vertexCounter] = positionCounter;
                    my_map[result] = positionCounter++;
                    my_positions.push_back(result);
                }
                else
                {
                    the_tri[vertexCounter] = it->second;
                }
            }
            else
            {

                bool find = false;
                for (size_t i=0; i<my_positions.size(); ++i)
                    if ( (result[0] == my_positions[i][0]) && (result[1] == my_positions[i][1])  && (result[2] == my_positions[i][2]))
                    {
                        find = true;
                        the_tri[vertexCounter] = static_cast<sofa::Index>(i);
                        break;
                    }

                if (!find)
                {
                    my_positions.push_back(result);
                    the_tri[vertexCounter] = static_cast<sofa::Index>(my_positions.size()-1);
                }
            }
            vertexCounter++;
        }
        else if (bufferWord == "endfacet")
        {
            this->addTriangle(my_triangles.wref(), the_tri);
            vertexCounter = 0;
        }
        else if (bufferWord == "endsolid" || bufferWord == "end")
        {
            break;
        }
    }

    dataFile.close();

    dmsg_info() << "done!" ;

    return true;
}
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
void MeshSTLLoader::addTriangle(type::vector<Triangle >& pTriangles, const Triangle& p)
{
    if (d_flipNormals.getValue())
    {
        Triangle revertP;
        std::reverse_copy(p.begin(), p.end(), revertP.begin());
        pTriangles.push_back(revertP);
    }
    else
    {
        pTriangles.push_back(p);
    }
}
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

// upstream MeshTopology.cpp:145; 7c18e95d5c5f2839079892c69e7d89a313c79603
void MeshTopology::EdgeUpdate::updateFromSurface()
{
    typedef MeshTopology::SeqTriangles SeqTriangles;
    typedef MeshTopology::SeqQuads     SeqQuads;
    typedef MeshTopology::SeqEdges     SeqEdges;

    std::map<Edge,unsigned int> edgeMap;
    unsigned int edgeIndex;
    SeqEdges& seqEdges = *topology->d_seqEdges.beginEdit();
    seqEdges.clear();
    const SeqTriangles& triangles = topology->getTriangles(); // do not use d_seqTriangles directly as it might not be up-to-date
    for (const auto& t : triangles)
    {
        std::map<Edge,unsigned int>::iterator ite;
        Edge e;
        for (unsigned int j=0; j<3; ++j)
        {
            unsigned int v1=t[(j+1)%3];
            unsigned int v2=t[(j+2)%3];
            // sort vertices in lexicographics order
            if (v1<v2)
                e=Edge(v1,v2);
            else
                e=Edge(v2,v1);
            ite=edgeMap.find(e);
            if (ite==edgeMap.end())
            {
                // edge not in edgeMap so create a new one
                edgeIndex=(unsigned int)seqEdges.size();
                edgeMap[e]=edgeIndex;
                // To be similar to TriangleSetTopologyContainer::createEdgeSetArray
                //d_seqEdges.push_back(e); Changed to have oriented edges on the border of the triangulation.
                seqEdges.push_back(Edge(v1, v2));
            }
//            else
//            {
//                edgeIndex=(*ite).second;
//            }
            //m_edgesInTriangle[i][j]=edgeIndex;
        }
    }

    const SeqQuads& quads = topology->getQuads(); // do not use d_seqQuads directly as it might not be up-to-date
    for (const auto& t : quads)
    {
        std::map<Edge,unsigned int>::iterator ite;
        Edge e;
        for (unsigned int j=0; j<4; ++j)
        {
            unsigned int v1=t[(j+1)%4];
            unsigned int v2=t[(j+2)%4];
            // sort vertices in lexicographics order
            if (v1<v2)
                e=Edge(v1,v2);
            else
                e=Edge(v2,v1);
            ite=edgeMap.find(e);
            if (ite==edgeMap.end())
            {
                // edge not in edgeMap so create a new one
                edgeIndex=(unsigned int)seqEdges.size();
                edgeMap[e]=edgeIndex;
                seqEdges.push_back(e);
            }
//            else
//            {
//                edgeIndex=(*ite).second;
//            }
            //m_edgesInQuad[i][j]=edgeIndex;
        }
    }

    topology->d_seqEdges.endEdit();
}

// upstream MeshTopology.cpp:820; 7c18e95d5c5f2839079892c69e7d89a313c79603
void MeshTopology::createEdgesAroundVertexArray ()
{
    const SeqEdges& edges = getEdges();
    m_edgesAroundVertex.clear();
    m_edgesAroundVertex.resize( nbPoints );
    /*....
        if (getNbTetrahedra() || getNbHexahedra())
        { // Unordered shells if the mesh is volumic
            for (unsigned int i = 0; i < edges.size(); ++i)
            {
                // adding edge i in the edge shell of both points
                m_edgesAroundVertex[ edges[i][0] ].push_back( i );
                m_edgesAroundVertex[ edges[i][1] ].push_back( i );
            }
        }
        else if (getNbTriangles() || getNbQuads())
        { // order edges in consistent order if possible (i.e. on manifold meshes)
            bool createdTriangleShell = getNbTriangles() && m_edgesInTriangle.empty();
            bool createdQuadShell = getNbQuads() && m_edgesInQuad.empty();
            if (createdTriangleShell) createTrianglesAroundVertexArray();
            if (createdQuadShell) createQuadsAroundVertexArray();
            const SeqTriangles& triangles = getTriangles();
            const SeqQuads& quads = getQuads();

        }
        else*/
    {
        // 1D mesh : put inbound edges before outbound edges
        for (unsigned int i = 0; i < edges.size(); ++i)
        {
            // adding edge i in the edge shell of both points
            m_edgesAroundVertex[ edges[i][0] ].push_back( i );
            m_edgesAroundVertex[ edges[i][1] ].push_back( i );
        }
    }
}

// upstream MeshTopology.cpp:857; 7c18e95d5c5f2839079892c69e7d89a313c79603
void MeshTopology::createEdgesInTriangleArray ()
{
    //const SeqEdges& edges = getEdges(); // do not use d_seqEdges directly as it might not be up-to-date
    const SeqTriangles& triangles = getTriangles(); // do not use d_seqTriangles directly as it might not be up-to-date
    m_edgesInTriangle.clear();
    m_edgesInTriangle.resize(triangles.size());
    for (unsigned int i = 0; i < triangles.size(); ++i)
    {
        const Triangle &t=triangles[i];
        // adding edge i in the edge shell of both points
        for (unsigned int j=0; j<3; ++j)
        {
            const EdgeID edgeIndex=getEdgeIndex(t[(j+1)%3],t[(j+2)%3]);
            assert(edgeIndex != InvalidID);
            m_edgesInTriangle[i][j]=edgeIndex;
        }
    }
}

// upstream MeshTopology.cpp:936; 7c18e95d5c5f2839079892c69e7d89a313c79603
void MeshTopology::createTrianglesAroundVertexArray ()
{
    const SeqTriangles& triangles = getTriangles(); // do not use d_seqTriangles directly as it might not be up-to-date
    m_trianglesAroundVertex.clear();
    m_trianglesAroundVertex.resize( nbPoints );

    for (unsigned int i = 0; i < triangles.size(); ++i)
    {
        // adding triangle i in the triangle shell of all points
        for (unsigned int j=0; j<3; ++j)
            m_trianglesAroundVertex[ triangles[i][j]  ].push_back( i );
    }
}

// upstream MeshTopology.cpp:1097; 7c18e95d5c5f2839079892c69e7d89a313c79603
void MeshTopology::createTrianglesAroundEdgeArray ()
{
    const SeqTriangles& triangles = getTriangles(); // do not use d_seqTriangles directly as it might not be up-to-date
    if (m_edgesInTriangle.empty())
        createEdgesInTriangleArray();
    m_trianglesAroundEdge.clear();
    m_trianglesAroundEdge.resize( getNbEdges());
    const vector< EdgesInTriangle > &tea=m_edgesInTriangle;
    unsigned int j;

    for (unsigned int i = 0; i < triangles.size(); ++i)
    {
        const Triangle &t=triangles[i];
        // adding triangle i in the triangle shell of all edges
        for (j=0; j<3; ++j)
        {
            if (d_seqEdges.getValue()[tea[i][j]][0] == t[(j + 1) % 3])
                m_trianglesAroundEdge[ tea[i][j] ].insert(m_trianglesAroundEdge[ tea[i][j] ].begin(), i); // triangle is on the left of the edge
            else
                m_trianglesAroundEdge[ tea[i][j] ].push_back( i ); // triangle is on the right of the edge
        }
    }
}
