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

// SubsetMapping.inl:182
void SubsetMapper::apply ( const MappingParams* /*mparams*/, OutDataVecCoord& dOut, const InDataVecCoord& dIn )
{
    const IndexArray& indices = d_indices.getValue();
    
    if (d_resizeToModel.getValue() || this->toModel->getSize() < indices.size())
    { 
        if (this->toModel->getSize() != indices.size()) 
        { 
            this->toModel->resize(Size(indices.size())); 
        } 
    }
    
    const InVecCoord& in = dIn.getValue();
    const OutVecCoord& out0 = this->toModel->read(core::vec_id::read_access::restPosition)->getValue();
    OutVecCoord& out = *dOut.beginEdit();
    const auto fromSize = in.size();

    out.resize(indices.size());
    for(std::size_t i = 0; i < out.size(); ++i)
    {
        if(indices[i] < fromSize)
            out[i] = in[ indices[i] ];
        else
            out[i] = out0[i];
    }

    dOut.endEdit();
}

// SubsetMapping.inl:257
void SubsetMapper::applyJT ( const core::ConstraintParams * /*cparams*/, InDataMatrixDeriv& dOut, const OutDataMatrixDeriv& dIn)
{
    const IndexArray& indices = d_indices.getValue();

    InMatrixDeriv& out = *dOut.beginEdit();
    const OutMatrixDeriv& in = dIn.getValue();
    const unsigned int fromSize = this->fromModel->getSize();

    typename Out::MatrixDeriv::RowConstIterator rowItEnd = in.end();

    for (typename Out::MatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        typename Out::MatrixDeriv::ColConstIterator colIt = rowIt.begin();
        typename Out::MatrixDeriv::ColConstIterator colItEnd = rowIt.end();

        // Creates a constraints if the input constraint is not empty.
        if (colIt != colItEnd)
        {
            typename In::MatrixDeriv::RowIterator o = out.writeLine(rowIt.index());

            while (colIt != colItEnd)
            {
                if(indices[colIt.index()] < fromSize)
                    o.addCol(indices[colIt.index()], colIt.val());
                ++colIt;
            }
        }
    }
    dOut.endEdit();
}
