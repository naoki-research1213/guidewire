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

// upstream:333; 7c18e95d5c5f2839079892c69e7d89a313c79603
template<class DataTypes>
void Spring<DataTypes>::addForce(const MechanicalParams*  mparams , DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv&  v )
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(v);

    WriteAccessor< DataVecDeriv > f1 = f;
    ReadAccessor< DataVecCoord > p1 = x;

    const DataVecCoord* extPosition = getExtPosition();
    if (!extPosition)
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }

    ReadAccessor< DataVecCoord > p0 = *extPosition;

    const VecReal& k = d_stiffness.getValue();
    const VecReal& k_a = d_angularStiffness.getValue();

    f1.resize(p1.size());

    if (d_recompute_indices.getValue())
    {
        recomputeIndices();
    }

    for (sofa::Index i = 0; i < m_indices.size(); i++)
    {
        const sofa::Index index = m_indices[i];
        sofa::Index ext_index = m_indices[i];
        if (useRestMState)
            ext_index = m_ext_indices[i];

        const auto stiffness = k[static_cast<std::size_t>(i < k.size()) * i];

        const auto activeDirections = d_activeDirections.getValue();

        // rigid case
        if constexpr (sofa::type::isRigidType<DataTypes>)
        {
            // translation
            if (i >= m_pivots.size())
            {
                CPos dx = p1[index].getCenter() - p0[ext_index].getCenter();
                // We filter the difference dx by setting to 0 the entries corresponding
                // to 0 values in d_activeDirections
                for (sofa::Size entryId = 0; entryId < spatial_dimensions; ++entryId)
                {
                    if (!activeDirections[entryId])
                        dx[entryId] = 0;
                }
                getVCenter(f1[index]) -= dx * stiffness;
            }
            else
            {
                CPos localPivot = p0[ext_index].getOrientation().inverseRotate(m_pivots[i] - p0[ext_index].getCenter());
                CPos rotatedPivot = p1[index].getOrientation().rotate(localPivot);
                CPos pivot2 = p1[index].getCenter() + rotatedPivot;
                CPos dx = pivot2 - m_pivots[i];
                getVCenter(f1[index]) -= dx * stiffness;
            }

            // rotation
            type::Quat<Real> dq = p1[index].getOrientation() * p0[ext_index].getOrientation().inverse();
            dq.normalize();

            type::Vec<3, Real> dir{type::NOINIT};
            Real angle = 0.;

            if (dq[3] < 0.)
            {
                dq = dq * -1.0;
            }

            if (dq[3] < 1.0)
                dq.quatToAxis(dir, angle);

            // We change the direction of the axis of rotation based on
            // the 0 values in d_activeDirections. This is equivalent
            // to senting to 0 the rotation axis components along x, y
            // and/or z, depending on the rotations we want to take into
            // account.
            for (sofa::Size entryId = spatial_dimensions; entryId < Deriv::total_size; ++entryId)
            {
                if (!activeDirections[entryId])
                    dir[entryId-spatial_dimensions] = 0;
            }

            const auto angularStiffness = k_a[static_cast<std::size_t>(i < k_a.size()) * i];
            getVOrientation(f1[index]) -= dir * angle * angularStiffness;
        }
        else // non-rigid implementation
        {
            Deriv dx = p1[index] - p0[ext_index];
            // We filter the difference dx by setting to 0 the entries corresponding
            // to 0 values in d_activeDirections
            for (sofa::Size entryId = 0; entryId < spatial_dimensions; ++entryId)
            {
                if (!activeDirections[entryId])
                    dx[entryId] = 0;
            }
            f1[index] -= dx * stiffness;
        }
    }
}
