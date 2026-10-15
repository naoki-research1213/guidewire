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


// LCPConstraintSolver.cpp:314
void ContactSolver::computeInitialGuess()
{
    sofa::helper::AdvancedTimer::StepVar vtimer("InitialGuess");

    const auto _mu = d_mu.getValue();
    const VecConstraintBlockInfo& constraintBlockInfo = m_constraintBlockInfo;
    const VecPersistentID& constraintIds = m_constraintIds;
    const int numContact = (_mu > 0.0) ? _numConstraints/3 : _numConstraints;

    for (int c=0; c<numContact; c++)
    {
        if (_mu>0.0)
        {
            (*_result)[3*c  ] = 0.0;
            (*_result)[3*c+1] = 0.0;
            (*_result)[3*c+2] = 0.0;
        }
        else
        {
            (*_result)[c] =  0.0;
        }
    }
    for (const ConstraintBlockInfo& info : constraintBlockInfo)
    {
        if (!info.hasId) continue;
        std::map<core::behavior::BaseLagrangianConstraint*, ConstraintBlockBuf>::const_iterator previt = _previousConstraints.find(info.parent);
        if (previt == _previousConstraints.end()) continue;
        const ConstraintBlockBuf& buf = previt->second;
        const int c0 = info.const0;
        const int nbl = (info.nbLines < buf.nbLines) ? info.nbLines : buf.nbLines;
        for (int c = 0; c < info.nbGroups; ++c)
        {
            std::map<PersistentID,int>::const_iterator it = buf.persistentToConstraintIdMap.find(constraintIds[info.offsetId + c]);
            if (it == buf.persistentToConstraintIdMap.end()) continue;
            const int prevIndex = it->second;
            if (prevIndex >= 0 && prevIndex+nbl <= (int) _previousForces.size())
            {
                for (int l=0; l<nbl; ++l)
                    (*_result)[c0 + c*nbl + l] = _previousForces[prevIndex + l];
            }
        }
    }
}

// LCPConstraintSolver.cpp:358
void ContactSolver::keepContactForcesValue()
{
    sofa::helper::AdvancedTimer::StepVar vtimer("KeepForces");
    const VecConstraintBlockInfo& constraintBlockInfo = m_constraintBlockInfo;
    const VecPersistentID& constraintIds = m_constraintIds;
    // store current force
    _previousForces.resize(_numConstraints);
    for (unsigned int c=0; c<_numConstraints; ++c)
        _previousForces[c] = (*_result)[c];
    // clear previous history
    for (auto& previousConstraint : _previousConstraints)
    {
        ConstraintBlockBuf& buf = previousConstraint.second;
        for (auto& it2 : buf.persistentToConstraintIdMap)
            it2.second = -1;
    }
    // fill info from current ids
    for (const ConstraintBlockInfo& info : constraintBlockInfo)
    {
        if (!info.parent) continue;
        if (!info.hasId) continue;
        ConstraintBlockBuf& buf = _previousConstraints[info.parent];
        const int c0 = info.const0;
        const int nbl = info.nbLines;
        buf.nbLines = nbl;
        for (int c = 0; c < info.nbGroups; ++c)
        {
            buf.persistentToConstraintIdMap[constraintIds[info.offsetId + c]] = c0 + c*nbl;
        }
    }
}
