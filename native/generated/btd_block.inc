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

constexpr size_t BSIZE=6;
using T=double;
    class TransposedBlock
    {
    public:
        const type::Mat<BSIZE,BSIZE,Real>& m;
        TransposedBlock(const type::Mat<BSIZE,BSIZE,Real>& m) : m(m) {}
        type::Vec<BSIZE,Real> operator*(const type::Vec<BSIZE,Real>& v)
        {
            return m.multTranspose(v);
        }
        type::Mat<BSIZE,BSIZE,Real> operator-() const
        {
            type::Mat<BSIZE,BSIZE,Real> r;
            for (Index i=0; i<BSIZE; i++)
                for (Index j=0; j<BSIZE; j++)
                    r(i,j)=-m(j,i);
            return r;
        }
    };

    class Block : public type::Mat<BSIZE,BSIZE,Real>
    {
    public:
        Index Nrows() const { return BSIZE; }
        Index Ncols() const { return BSIZE; }
        void resize(Index, Index)
        {
            this->clear();
        }
        const T& element(Index i, Index j) const { return (*this)(i,j); }
        void set(Index i, Index j, const T& v) { (*this)(i,j) = v; }
        void add(Index i, Index j, const T& v) { (*this)(i,j) += v; }
        void operator=(const type::Mat<BSIZE,BSIZE,Real>& v)
        {
            type::Mat<BSIZE,BSIZE,Real>::operator=(v);
        }
        type::Mat<BSIZE,BSIZE,Real> operator-() const
        {
            type::Mat<BSIZE,BSIZE,Real> r;
            for (Index i=0; i<BSIZE; i++)
                for (Index j=0; j<BSIZE; j++)
                    r(i,j)=-(*this)(i,j);
            return r;
        }
        type::Mat<BSIZE,BSIZE,Real> operator-(const type::Mat<BSIZE,BSIZE,Real>& m) const
        {
            return type::Mat<BSIZE,BSIZE,Real>::operator-(m);
        }
        type::Vec<BSIZE,Real> operator*(const type::Vec<BSIZE,Real>& v)
        {
            return type::Mat<BSIZE,BSIZE,Real>::operator*(v);
        }
        type::Mat<BSIZE,BSIZE,Real> operator*(const type::Mat<BSIZE,BSIZE,Real>& m)
        {
            return sofa::type::operator*(*this, m);
        }
        type::Mat<BSIZE,BSIZE,Real> operator*(const TransposedBlock& mt)
        {
            return operator*(mt.m.transposed());
        }
        TransposedBlock t() const
        {
            return TransposedBlock(*this);
        }
        Block i() const
        {
            Block r;
            const bool canInvert = r.invert(*this);
            assert(canInvert);
            SOFA_UNUSED(canInvert);
            return r;
        }
    };

