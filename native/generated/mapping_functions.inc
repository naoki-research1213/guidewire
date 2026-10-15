/******************************************************************************
*                              BeamAdapter plugin                             *
*                  (c) 2006 Inria, University of Lille, CNRS                  *
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
* Authors: see Authors.md                                                     *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
//
// C++ Implementation : AdaptiveBeamMapping
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#pragma once

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////

// AdaptiveBeamMapping.inl:567
void Mapping::computeIdxAndBaryCoordsForAbs(unsigned int &b, Real &xBary, const Real &xAbs )
{
    const InReal xAbsInput = (InReal) xAbs;
    InReal xBaryOutput = (InReal) xBary;

    l_adaptativebeamInterpolation->getBeamAtCurvAbs(xAbsInput,b,xBaryOutput);
    xBary = (Real) xBaryOutput;
}

// AdaptiveBeamMapping.inl:578
void Mapping::computeDistribution()
{
    //The normal procedure is givng "points", computeDistribution comppute theses points and put
    //on the baryPoints which are "pointBeamDistribution"
    //If the mapping is a barycentric one, that mean a list of baryPoints are already setted
    //Then no need to recompute this
    if(!m_isBarycentricMapping)
    {
        bool curvAbs = d_useCurvAbs.getValue();
        const type::vector<Vec3>& points = d_points.getValue();
        m_pointBeamDistribution.clear();

        const unsigned int numBeams = l_adaptativebeamInterpolation->getNumBeams();
        if(numBeams==0)
        {
            if (this->f_printLog.getValue())
                msg_info() <<"No beams found in adaptBeamInterpolation in BeamInterpolation named"<<l_adaptativebeamInterpolation->getName();
            return;
        }

        if (curvAbs)
        {
            double ptsPerBeam = d_nbPointsPerBeam.getValue();
            if(ptsPerBeam)
            {	// Recreating the distribution based on the current sampling of the beams
                WriteAccessor<Data<type::vector<Real>>> waSegmentsCurvAbs = d_segmentsCurvAbs;
                waSegmentsCurvAbs.clear();

                double step = 1.0 / ptsPerBeam;
                double posInBeam = 0;
                InReal segStart, segEnd, segLength;
                for(unsigned int b=0; b< numBeams; ++b)
                {
                    l_adaptativebeamInterpolation->getAbsCurvXFromBeam(b, segStart, segEnd);
                    segLength = segEnd - segStart;

                    for (; posInBeam <= 1.0; posInBeam += step)
                    {
                        waSegmentsCurvAbs.push_back(segStart + segLength * posInBeam);
                        m_pointBeamDistribution.emplace_back( b, Vec3{ posInBeam, 0.0, 0.0} );
                    }

                    posInBeam -= 1.0;
                }

                if (numBeams && fabs(posInBeam - step) > 0.01)
                {
                    // Last point
                    waSegmentsCurvAbs.push_back(segEnd);
                    m_pointBeamDistribution.emplace_back( numBeams - 1, Vec3{ 1.0, 0.0, 0.0} );
                }

                BaseMeshTopology* topo = this->getContext()->getMeshTopology();
                if(topo)
                {
                    topo->clear();
                    int nbEdges = static_cast<int>(m_pointBeamDistribution.size()) - 1;
                    for(int i=0; i<nbEdges; ++i)
                        topo->addEdge(i, i+1);
                }
            }
            else
            {
                // We use the points Data
                m_pointBeamDistribution.reserve(points.size());
                for (unsigned int i=0; i<points.size(); i++)
                {
                    unsigned int b=0;
                    Real xAbs = points[i][0];
                    Real xBary = 0.0;

                    computeIdxAndBaryCoordsForAbs(b, xBary,  xAbs );

                    m_pointBeamDistribution.emplace_back( b, Vec3{ xBary, points[i][1], points[i][2]});
                }
            }
        }
        else
        {
            m_pointBeamDistribution.reserve(points.size());
            for (unsigned int i=0; i<points.size(); i++)
            {
                const unsigned int beamId = (int)floor(points[i][0]);
                msg_warning_when(beamId > numBeams-1) << "Points[" << i << "][0] = " << beamId << " is defined outside of the beam length";
                m_pointBeamDistribution.emplace_back( beamId, Vec3{ points[i][0] - floor(points[i][0]), points[i][1], points[i][2]});
            }
        }
    }
}

// AdaptiveBeamMapping.inl:783
void Mapping::applyJTonPoint(unsigned int i, const Deriv& finput, SpatialVector& FNode0output, SpatialVector& FNode1output, const InVecCoord& x)
{
    const PosPointDefinition&  pointBeamDistribution = m_pointBeamDistribution[i];
    const Vec3 localPos(0.,pointBeamDistribution.baryPoint[1],pointBeamDistribution.baryPoint[2]);
    const Vec3 Fin(finput[0], finput[1], finput[2]);
    l_adaptativebeamInterpolation->MapForceOnNodeUsingSpline(pointBeamDistribution.beamId, pointBeamDistribution.baryPoint[0], localPos, x, Fin, FNode0output, FNode1output );
}

// AdaptiveBeamMapping.inl:720
void Mapping::applyJonPoint(unsigned int i, SpatialVector& VNode0input, SpatialVector& VNode1input, Deriv& vOutput, const InVecCoord& x)
{

    // 1. get the curvilinear abs;
    const PosPointDefinition&  ppd = m_pointBeamDistribution[i];

    // 2. get the indices
    unsigned int IdxNode0, IdxNode1;
    l_adaptativebeamInterpolation->getNodeIndices(ppd.beamId,IdxNode0,IdxNode1);

    // 3. get the transform to DOF in global frame from local frame
    Transform DOF0Global_H_local0, DOF1Global_H_local1;
    l_adaptativebeamInterpolation->getDOFtoLocalTransformInGlobalFrame(ppd.beamId, DOF0Global_H_local0, DOF1Global_H_local1, x);

    // 4. project the velocities in local frame:
    SpatialVector v_local0, v_local1;
    v_local0 = DOF0Global_H_local0.inversed()*VNode0input;
    v_local1 = DOF1Global_H_local1.inversed()*VNode1input;

    // 5. Computes the local velocities of the 4 points of the spline
    Real L = l_adaptativebeamInterpolation->getLength(ppd.beamId);
    Vec3 lever(L/3,0,0);
    Vec3 V0, V1, V2, V3;
    V0 = v_local0.getLinearVelocity();
    V1 = V0 - lever.cross(v_local0.getAngularVelocity());
    V3 = v_local1.getLinearVelocity();
    V2 = V3 + lever.cross(v_local1.getAngularVelocity());

    // 6. Rotate back the vectors in the global frame
    V0 = DOF0Global_H_local0.getOrientation().rotate(V0);
    V1 = DOF0Global_H_local0.getOrientation().rotate(V1);
    V2 = DOF1Global_H_local1.getOrientation().rotate(V2);
    V3 = DOF1Global_H_local1.getOrientation().rotate(V3);

    const Vec3 localPos(0.,ppd.baryPoint[1],ppd.baryPoint[2]);
    if(localPos.norm() > L*1e-10)
    {
        lever = localPos;

        // 7. compute the effect of the angular velocities on the points that are not aligned along the center line
        Vec3 DV0, DV3;
        DV0 = DOF0Global_H_local0.getOrientation().rotate( - lever.cross( v_local0.getAngularVelocity() ) );
        DV3 = DOF1Global_H_local1.getOrientation().rotate( - lever.cross( v_local1.getAngularVelocity() ) );

        // uses spline to interpolate:
        Real bx = ppd.baryPoint[0];
        Real a0=(1-bx)*(1-bx)*(1-bx);
        Real a1=3*bx*(1-bx)*(1-bx);
        Real a2=3*bx*bx*(1-bx);
        Real a3=bx*bx*bx;
        vOutput = V0*a0 + V1*a1 + V2*a2 + V3*a3 + DV0*(a0+a1) + DV3*(a2+a3);

    }
    else
    {
        // uses spline to interpolate:
        Real bx = ppd.baryPoint[0];
        vOutput = V0*(1-bx)*(1-bx)*(1-bx) + V1*3*bx*(1-bx)*(1-bx) + V2*3*bx*bx*(1-bx) + V3*bx*bx*bx;
    }
}

// AdaptiveBeamMapping.inl:201
void Mapping::apply(const MappingParams* mparams, Data<VecCoord>& dOut, const Data<InVecCoord>& dIn)
{
    SCOPED_TIMER("AdaptiveBeamMapping_Apply");

    auto out = sofa::helper::getWriteOnlyAccessor(dOut);
    const InVecCoord& in = dIn.getValue();

    // When using an adaptatif controller, one need to redistribute the points at each time step
    {
        SCOPED_TIMER("computeDistribution");
        if (d_useCurvAbs.getValue() && !d_contactDuplicate.getValue())
            computeDistribution();
    }

    if (!m_isSubMapping)
    {
        SCOPED_TIMER("resizeToModel&Out");
        if (d_nbPointsPerBeam.getValue() > 0)
        {
            this->toModel->resize(m_pointBeamDistribution.size());
            out.resize(m_pointBeamDistribution.size());
        }
        else
        {
            this->toModel->resize(d_points.getValue().size());
            out.resize(d_points.getValue().size());
        }
    }

    MultiVecCoordId x = sofa::core::vec_id::write_access::position;
    MultiVecCoordId xfree = sofa::core::vec_id::write_access::freePosition;

    const ConstMultiVecCoordId &xId = mparams->x();
    sofa::core::ConstVecCoordId xtest = xId.getId(this->fromModel);

    if(xtest == xfree.getId(this->fromModel))
    {
        VecCoordId xfreeIn = sofa::core::vec_id::write_access::freePosition;
        l_adaptativebeamInterpolation->updateBezierPoints(in, xfreeIn);
    }
    else if(xtest == x.getId(this->fromModel))
    {
        VecCoordId positionIn = sofa::core::vec_id::write_access::position;
        l_adaptativebeamInterpolation->updateBezierPoints(in, positionIn);
    }

    AdvancedTimer::stepBegin("computeNewInterpolation");

    // effective computation for each element in pointBeamDistribution
    auto apply_impl = [&](const auto& range)
    {
        auto elementID = std::distance(m_pointBeamDistribution.begin(), range.start);
        for (auto it = range.start; it != range.end; ++it, ++elementID)
        {
            const auto& pointBeamDistribution = m_pointBeamDistribution[elementID];

            Vec<3, InReal> pos(sofa::type::NOINIT);
            const Vec3 localPos(0., pointBeamDistribution.baryPoint[1], pointBeamDistribution.baryPoint[2]);
            l_adaptativebeamInterpolation->interpolatePointUsingSpline(pointBeamDistribution.beamId, pointBeamDistribution.baryPoint[0], localPos, in, pos, false, xtest);

            if (m_isSubMapping)
            {
                if (!m_idPointSubMap.empty())
                    out[m_idPointSubMap[elementID]] = pos;
            }
            else
                out[elementID] = pos;
        }

    };

    
    // HACK for init: In case the number of output points is bigger to the number of distribution, set all points to 0
    for (std::size_t i = m_pointBeamDistribution.size(); i < d_points.getValue().size(); i++)
    {
        out[i].clear();
    }


    const bool multithreading = d_parallelMapping.getValue();

    const simulation::ForEachExecutionPolicy execution = multithreading ?
        simulation::ForEachExecutionPolicy::PARALLEL :
        simulation::ForEachExecutionPolicy::SEQUENTIAL;

    simulation::TaskScheduler* taskScheduler = simulation::MainTaskSchedulerFactory::createInRegistry();
    assert(taskScheduler);

    simulation::forEachRange(execution, *taskScheduler, m_pointBeamDistribution.begin(), m_pointBeamDistribution.end(), apply_impl);


    AdvancedTimer::stepEnd("computeNewInterpolation");
}

// AdaptiveBeamMapping.inl:296
void Mapping::applyJ(const MappingParams* mparams, Data<VecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    SOFA_UNUSED(mparams);
    SCOPED_TIMER("AdaptiveBeamMapping_applyJ");

    auto out = sofa::helper::getWriteOnlyAccessor(dOut);
    const InVecDeriv& in= dIn.getValue();

    const Data<InVecCoord>& dataInX = *this->getFromModel()->read(sofa::core::vec_id::read_access::position);
    const InVecCoord& x = dataInX.getValue();

    if (d_useCurvAbs.getValue() && !d_contactDuplicate.getValue())
        computeDistribution();

    // should not be necessary if apply() was called first
    if (!m_isSubMapping)
    {
        if (d_nbPointsPerBeam.getValue() > 0)
        {
            out.resize(m_pointBeamDistribution.size());
        }
        else
        {
            out.resize(d_points.getValue().size());
        }
    }

    // effective computation for each element in pointBeamDistribution
    auto applyJ_impl = [&](const auto& range)
    {
        auto elementID = std::distance(m_pointBeamDistribution.begin(), range.start);
        for (auto it = range.start; it != range.end; ++it, ++elementID)
        {
            const auto& pointBeamDistribution = m_pointBeamDistribution[elementID];

            unsigned int IdxNode0, IdxNode1;
            l_adaptativebeamInterpolation->getNodeIndices(pointBeamDistribution.beamId, IdxNode0, IdxNode1);

            SpatialVector vDOF0{ In::getDRot(in[IdxNode0]), In::getDPos(in[IdxNode0]) };
            SpatialVector vDOF1{ In::getDRot(in[IdxNode1]), In::getDPos(in[IdxNode1]) };

            Deriv vResult(sofa::type::NOINIT);

            applyJonPoint(elementID, vDOF0, vDOF1, vResult, x);

            if (m_isSubMapping)
            {
                if (m_idPointSubMap.size() > 0)
                    out[m_idPointSubMap[elementID]] = vResult;
            }
            else
                out[elementID] = vResult;
        }
    };

    const bool multithreading = d_parallelMapping.getValue();

    const simulation::ForEachExecutionPolicy execution = multithreading ?
        simulation::ForEachExecutionPolicy::PARALLEL :
        simulation::ForEachExecutionPolicy::SEQUENTIAL;

    simulation::TaskScheduler* taskScheduler = simulation::MainTaskSchedulerFactory::createInRegistry();
    assert(taskScheduler);

    simulation::forEachRange(execution, *taskScheduler, m_pointBeamDistribution.begin(), m_pointBeamDistribution.end(), applyJ_impl);
}

// AdaptiveBeamMapping.inl:408
void Mapping::applyJT(const core::ConstraintParams* cparams, Data<InMatrixDeriv>& dOut, const Data<OutMatrixDeriv>& dIn)
{
    SOFA_UNUSED(cparams);

    SCOPED_TIMER("AdaptiveBeamMapping_ApplyJT");
    auto out = sofa::helper::getWriteOnlyAccessor(dOut);
    const OutMatrixDeriv& in = dIn.getValue();
    const Data<InVecCoord>& dataInX = *this->getFromModel()->read(sofa::core::vec_id::read_access::position);
    const InVecCoord& x = dataInX.getValue();

    //////////// What's for ?? it seems not useful//////////
    bool proximity_lever = false;
    if (d_proximity.getValue() > 0.0)
        proximity_lever = true;

    if (proximity_lever && this->f_printLog.getValue() )
        msg_warning() <<"The constraints are contact at the surface of the beam (not at their center)";
    /////////////////////////////////////

    typename Out::MatrixDeriv::RowConstIterator rowItEnd = in.end();
    for (typename Out::MatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        typename Out::MatrixDeriv::ColConstIterator colItEnd = rowIt.end();
        for (typename Out::MatrixDeriv::ColConstIterator colIt = rowIt.begin(); colIt != colItEnd; ++colIt)
        {
            typename In::MatrixDeriv::RowIterator o = out.wref().writeLine(rowIt.index());
            unsigned int indexIn = colIt.index();
            const Deriv data = colIt.val();

            if (m_isSubMapping)
            {
                // look if we get the indexIn in the idPointSubMap:
                unsigned int i=0;
                while( i<m_idPointSubMap.size() && m_idPointSubMap[i]!=indexIn)
                    i++;

                if (i<m_idPointSubMap.size())
                    indexIn = i;
                else
                    continue;
            }

            if(indexIn<m_pointBeamDistribution.size()){
                const PosPointDefinition&  ppd = m_pointBeamDistribution[indexIn];
                unsigned int IdxNode0, IdxNode1;
                l_adaptativebeamInterpolation->getNodeIndices(ppd.beamId,IdxNode0,IdxNode1);

                SpatialVector FNode0, FNode1;
                applyJTonPoint(indexIn, data, FNode0, FNode1, x);

                // Compute the mapped Constraint on the beam nodes
                InDeriv direction0;
                In::setDPos(direction0,FNode0.getForce());
                In::setDRot(direction0,FNode0.getTorque());
                InDeriv direction1;
                In::setDPos(direction1,FNode1.getForce());
                In::setDRot(direction1,FNode1.getTorque());

                o.addCol(IdxNode0, direction0);
                o.addCol(IdxNode1, direction1);
            }
            else
            {
                if ( this->f_printLog.getValue() )
                    msg_warning() <<"Wrong index in VecConst in";
                break;
            }
        }
    }
}
