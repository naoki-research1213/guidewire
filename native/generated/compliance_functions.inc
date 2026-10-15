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

template<class RMatrix,class JMatrix>
bool BTD::addJMInvJt(RMatrix& result, JMatrix& J, double fact)
{
    const Index Jcols = J.colSize();
    if (Jcols != Minv.rowSize())
    {
        msg_error() << "AddJMInvJt: incompatible J matrix size.";
        return false;
    }

    if (this->notMuted())
    {
        std::stringstream tmpStr;
        tmpStr<< "C = ["<<msgendl;
        for  (Index mr=0; mr<Minv.rowSize(); mr++)
        {
            tmpStr<<" "<<msgendl;
            for (Index mc=0; mc<Minv.colSize(); mc++)
            {
                tmpStr<<" "<< getMinvElement(mr,mc);
            }
        }
        tmpStr << "];"<<msgendl;

        tmpStr<< "J = ["<<msgendl;
        for  (Index jr=0; jr<J.rowSize(); jr++)
        {
            tmpStr<<" "<<msgendl;
            for (Index jc=0; jc<J.colSize(); jc++)
            {
                tmpStr<<" "<< J.element(jr, jc) ;
            }
        }
        tmpStr << "];";
        msg_info() << tmpStr.str();
    }

    std::stringstream tmpStr2;
    const typename JMatrix::LineConstIterator jitend = J.end();
    for (typename JMatrix::LineConstIterator jit1 = J.begin(); jit1 != jitend; ++jit1)
    {
        Index row1 = jit1->first;
        for (typename JMatrix::LineConstIterator jit2 = jit1; jit2 != jitend; ++jit2)
        {
            Index row2 = jit2->first;
            double acc = 0.0;
            for (typename JMatrix::LElementConstIterator i1 = jit1->second.begin(), i1end = jit1->second.end(); i1 != i1end; ++i1)
            {
                const Index col1 = i1->first;
                const double val1 = i1->second;
                for (typename JMatrix::LElementConstIterator i2 = jit2->second.begin(), i2end = jit2->second.end(); i2 != i2end; ++i2)
                {
                    const Index col2 = i2->first;
                    const double val2 = i2->second;
                    acc += val1 * getMinvElement(col1,col2) * val2;
                }
            }

            if (this->notMuted())
            {
                tmpStr2 << "W("<<row1<<","<<row2<<") += "<<acc<<" * "<<fact<<msgendl;
            }

            acc *= fact;
            result.add(row1,row2,acc);
            if (row1!=row2)
                result.add(row2,row1,acc);
        }
    }

    return true;
}
