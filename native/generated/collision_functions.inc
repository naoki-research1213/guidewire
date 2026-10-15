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
#define SOFA_COMPONENT_COLLISION_LOCALMINDISTANCE_CPP

// upstream:343; 7c18e95d5c5f2839079892c69e7d89a313c79603
int LocalMinDistance::computeIntersection(Triangle& e2, Point& e1, OutputVector* contacts, const core::collision::Intersection* currentIntersection)
{
    if(!e1.isActive(e2.getCollisionModel()))
        return 0;

    const SReal alarmDist = currentIntersection->getAlarmDistance() + e1.getContactDistance() + e2.getContactDistance();

    static_assert(std::is_same_v<Triangle::Coord, Point::Coord>, "Data mismatch");
    using Real = Triangle::Coord::value_type;

    const Triangle::Coord AB = e2.p2()-e2.p1();
    const Triangle::Coord AC = e2.p3()-e2.p1();
    const auto AP = e1.p() -e2.p1();
    MatNoInit<2, 2, Real> A;
    VecNoInit<2, Real> b;

    A(0,0) = AB*AB;
    A(1,1) = AC*AC;
    A(0,1) = A(1,0) = AB*AC;
    b[0] = AP*AB;
    b[1] = AP*AC;

    const Real det = type::determinant(A);

    const Real alpha=(b[0]*A(1,1)-b[1]*A(0,1))/det;
    const Real beta=(b[1]*A(0,0)-b[0]*A(1,0))/det;

    if (alpha < 0.000001 ||
            beta  < 0.000001 ||
            alpha + beta  > 0.999999)
        return 0;

    const Point::Coord P = e1.p();
    const Triangle::Coord Q = e2.p1()+AB*alpha+AC*beta;
    const auto PQ = Q-P;

    if (PQ.norm2() >= alarmDist*alarmDist)
        return 0;

    const auto QP = -PQ;

    // filter for LMD

    if (!d_useLMDFilters.getValue())
    {
        if (!testValidity(e1, PQ))
            return 0;

        if (!testValidity(e2, QP))
            return 0;
    }

    //end filter

    contacts->resize(contacts->size()+1);
    DetectionOutput *detection = &*(contacts->end()-1);

#ifdef SOFA_DETECTIONOUTPUT_FREEMOTION
    if (e1.hasFreePosition() && e2.hasFreePosition())
    {
        type::Vec3 Pfree,Qfree,ABfree,ACfree;
        ABfree = e2.p2Free()-e2.p1Free();
        ACfree = e2.p3Free()-e2.p1Free();
        Pfree = e1.pFree();
        Qfree = e2.p1Free() + ABfree * alpha + ACfree * beta;

        detection->freePoint[0] = Qfree;
        detection->freePoint[1] = Pfree;
    }
#endif

    const double contactDist = currentIntersection->getContactDistance() + e1.getContactDistance() + e2.getContactDistance();

    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e2, e1);
    detection->id = e1.getIndex();
    detection->point[0] = Q;
    detection->point[1] = P;
    detection->normal = QP;
    detection->value = detection->normal.norm();
    detection->normal /= detection->value;
    detection->value -= contactDist;
    return 1;
}

// upstream:628; 7c18e95d5c5f2839079892c69e7d89a313c79603
int LocalMinDistance::computeIntersection(Line& e2, Point& e1, OutputVector* contacts, const core::collision::Intersection* currentIntersection)
{
    static_assert(std::is_same_v<Line::Coord, Point::Coord>, "Data mismatch");
    using Real = Triangle::Coord::value_type;

    if(!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return 0;

    const SReal alarmDist = currentIntersection->getAlarmDistance() + e1.getContactDistance() + e2.getContactDistance();
    const Line::Coord AB = e2.p2()-e2.p1();
    const Line::Coord AP = e1.p()-e2.p1();

    if (AB.norm()<0.000000000001*AP.norm())
    {
        return 0;
    }

    const Real A=AB*AB;
    const Real b=AP*AB;

    const Real alpha = b / A;
    if (alpha < 0.000001 || alpha > 0.999999)
        return 0;

    const Point::Coord& P = e1.p();
    const Line::Coord Q = e2.p1() + AB * alpha;
    const auto PQ = Q - P;

    if (PQ.norm2() >= alarmDist*alarmDist)
        return 0;

    const auto QP = -PQ;

    // filter for LMD
    if (!d_useLMDFilters.getValue())
    {
        if (!testValidity(e1, PQ))
            return 0;

        if (!testValidity(e2, QP))
            return 0;
    }

    // end filter

    contacts->resize(contacts->size()+1);
    DetectionOutput *detection = &*(contacts->end()-1);

#ifdef SOFA_DETECTIONOUTPUT_FREEMOTION
    if (e1.hasFreePosition() && e2.hasFreePosition())
    {
        type::Vec3 ABfree = e2.p2Free() - e2.p1Free();
        type::Vec3 Pfree = e1.pFree();
        type::Vec3 Qfree = e2.p1Free() + ABfree * alpha;

        detection->freePoint[0] = Qfree;
        detection->freePoint[1] = Pfree;
    }
#endif

    const double contactDist = currentIntersection->getContactDistance() + e1.getContactDistance() + e2.getContactDistance();

    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e2, e1);
    detection->id = e1.getIndex();
    detection->point[0]=Q;
    detection->point[1]=P;
    detection->normal=QP;
    detection->value = detection->normal.norm();
    detection->normal /= detection->value;
    detection->value -= contactDist;

    return 1;
}

// upstream:161; 7c18e95d5c5f2839079892c69e7d89a313c79603
int LocalMinDistance::computeIntersection(Line& e1, Line& e2, OutputVector* contacts, const core::collision::Intersection* currentIntersection)
{

    if(!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
    {
        dmsg_info_when(EMIT_EXTRA_DEBUG_MESSAGE)
            <<" not activated" ;
        return 0;
    }

    const SReal alarmDist = currentIntersection->getAlarmDistance() + e1.getContactDistance() + e2.getContactDistance();

    // E1 => A-->B
    // E2 => C-->D
    const Line::Coord AB = e1.p2()-e1.p1();
    const Line::Coord CD = e2.p2()-e2.p1();
    const Line::Coord AC = e2.p1()-e1.p1();
    Matrix2 A;
    Vec2 b;

    A(0,0) = AB*AB;
    A(1,1) = CD*CD;
    A(0,1) = A(1,0) = -CD*AB;
    b[0] = AB*AC;
    b[1] = -CD*AC;
    const double det = type::determinant(A);

    double alpha;
    double beta;

    // If lines are not parallel
    if (det < -1.0e-30 || det > 1.0e-30)
    {
        // compute the parametric coordinates along the line
        alpha = (b[0]*A(1,1) - b[1]*A(0,1))/det;
        beta  = (b[1]*A(0,0) - b[0]*A(1,0))/det;

        // if parameters are outside of ]0,1[ then there is no intersection
        // the intersection is outside of the edge supporting point.
        if (alpha < 1e-15 || alpha > (1.0-1e-15) ||
            beta  < 1e-15  || beta  > (1.0-1e-15) )
            return 0;
    }
    else
    {
        // lines are parallel,
        // the alpha/beta parameters are set to 0.5 so the collision
        // output is in the middle of the line
        alpha = 0.5;
        beta = 0.5;
    }


    Vec3 P,Q,PQ;
    P = e1.p1() + AB * alpha;
    Q = e2.p1() + CD * beta;
    PQ = Q-P;

    // If the geometric distance between P and Q is higher than the alarm distance/
    if (PQ.norm2() >= alarmDist*alarmDist)
        return 0;

    // filter for LMD //

    if (!d_useLMDFilters.getValue())
    {
        if (!testValidity(e1, PQ))
        {
            dmsg_info_when(EMIT_EXTRA_DEBUG_MESSAGE)
                 <<" testValidity rejected for the first segment" ;
            return 0;
        }

        const Vec3 QP = -PQ;

        if (!testValidity(e2, QP))
        {
            dmsg_info_when(EMIT_EXTRA_DEBUG_MESSAGE)
                <<" testValidity rejected for the second segment";
            return 0;
        }
    }

    contacts->resize(contacts->size() + 1);
    DetectionOutput *detection = &*(contacts->end() - 1);

#ifdef SOFA_DETECTIONOUTPUT_FREEMOTION

    if (e1.hasFreePosition() && e2.hasFreePosition())
    {
        type::Vec3 Pfree, Qfree, ABfree, CDfree;
        ABfree = e1.p2Free()-e1.p1Free();
        CDfree = e2.p2Free()-e2.p1Free();
        Pfree = e1.p1Free() + ABfree * alpha;
        Qfree = e2.p1Free() + CDfree * beta;
        detection->freePoint[0] = Pfree;
        detection->freePoint[1] = Qfree;
    }

#endif

    const double contactDist = currentIntersection->getContactDistance() + e1.getContactDistance() + e2.getContactDistance();

    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
    detection->id = (e1.getCollisionModel()->getSize() > e2.getCollisionModel()->getSize()) ? e1.getIndex() : e2.getIndex();
    detection->point[0] = P;
    detection->point[1] = Q;
    detection->normal = PQ;
    detection->value = detection->normal.norm();
    detection->normal /= detection->value;
    detection->value -= contactDist;

    return 1;
}

// upstream:1191; 7c18e95d5c5f2839079892c69e7d89a313c79603
bool LocalMinDistance::testValidity(Point &p, const Vec3 &PQ) const
{
    if (!d_filterIntersection.getValue())
        return true;

    const Vec3 pt = p.p();

    const sofa::simulation::Node* node = dynamic_cast<sofa::simulation::Node*>(p.getCollisionModel()->getContext());
    if ( !node || !(node->get< LineCollisionModel<sofa::defaulttype::Vec3Types> >()) )
        return true;

    BaseMeshTopology* topology = p.getCollisionModel()->getCollisionTopology();
    const auto& x =(p.getCollisionModel()->getMechanicalState()->read(core::vec_id::read_access::position)->getValue());

    const auto& trianglesAroundVertex = topology->getTrianglesAroundVertex(p.getIndex());
    const auto& edgesAroundVertex = topology->getEdgesAroundVertex(p.getIndex());
    Vec3 nMean;

    for (unsigned int t : trianglesAroundVertex)
    {
        const auto& ptr = topology->getTriangle(t);
        Vec3 nCur = (x[ptr[1]]-x[ptr[0]]).cross(x[ptr[2]]-x[ptr[0]]);
        nCur.normalize();
        nMean += nCur;
    }

    if (trianglesAroundVertex.empty())
    {
        for (unsigned int e : edgesAroundVertex)
        {
            const auto& ped = topology->getEdge(e);
            Vec3 l = (pt - x[ped[0]]) + (pt - x[ped[1]]);
            l.normalize();
            nMean += l;
        }
    }



    if (nMean.norm()> 0.0000000001)
    {
        /// validity test with nMean, except if bothSide
        const PointCollisionModel<sofa::defaulttype::Vec3Types> *pM = p.getCollisionModel();
        const bool bothSide_computation = pM->d_bothSide.getValue();
        nMean.normalize();
        if (dot(nMean, PQ) < -d_angleCone.getValue() * PQ.norm() && !bothSide_computation)
        {
            return false;
        }
    }

    for (unsigned int e : edgesAroundVertex)
    {
        const auto& ped = topology->getEdge(e);
        Vec3 l = (pt - x[ped[0]]) + (pt - x[ped[1]]);
        l.normalize();
        double computedAngleCone = dot(nMean , l) * d_coneFactor.getValue();
        if (computedAngleCone<0)
            computedAngleCone=0.0;
        computedAngleCone+=d_angleCone.getValue();
        if (dot(l , PQ) < -computedAngleCone*PQ.norm())
        {
            return false;
        }
    }

    return true;
}

// upstream:1260; 7c18e95d5c5f2839079892c69e7d89a313c79603
bool LocalMinDistance::testValidity(Line &l, const Vec3 &PQ) const
{
    if (!d_filterIntersection.getValue())
        return true;

    const LineCollisionModel<sofa::defaulttype::Vec3Types> *lM = l.getCollisionModel();
    const bool bothSide_computation = lM->d_bothSide.getValue();

    Vec3 n1;

    const Line::Coord &pt1 = l.p1();
    const Line::Coord &pt2 = l.p2();

    Line::Coord AB = pt2 - pt1;
    AB.normalize();

    BaseMeshTopology* topology = l.getCollisionModel()->getCollisionTopology();
    const auto& x =(l.getCollisionModel()->getMechanicalState()->read(core::vec_id::read_access::position)->getValue());
    const auto& trianglesAroundEdge = topology->getTrianglesAroundEdge(l.getIndex());

    if ( trianglesAroundEdge.size() == 2)
    {
        // which triangle is left ?
        const BaseMeshTopology::Triangle& triangle0 = topology->getTriangle(trianglesAroundEdge[0]);
        bool triangle0_is_left=false;
        if ( (l.i1()==triangle0[0]&&l.i2()==triangle0[1]) || (l.i1()==triangle0[1]&&l.i2()==triangle0[2]) || (l.i1()==triangle0[2]&&l.i2()==triangle0[0]) )
        {
            triangle0_is_left=true;
        }


        // compute the normal of the triangle situated on the right
        const BaseMeshTopology::Triangle& triangleRight = triangle0_is_left ? topology->getTriangle(trianglesAroundEdge[1]): topology->getTriangle(trianglesAroundEdge[0]);
        n1 = cross(x[triangleRight[1]]-x[triangleRight[0]], x[triangleRight[2]]-x[triangleRight[0]]);
        n1.normalize();
        Vec3 nMean=n1;
        Vec3 t1=cross(n1, AB);
        t1.normalize(); // necessary ?

        // compute the normal of the triangle situated on the left
        const BaseMeshTopology::Triangle& triangleLeft = triangle0_is_left ? topology->getTriangle(trianglesAroundEdge[0]): topology->getTriangle(trianglesAroundEdge[1]);
        Vec3 n2=cross(x[triangleLeft[1]]-x[triangleLeft[0]], x[triangleLeft[2]]-x[triangleLeft[0]]);
        n2.normalize();
        nMean += n2;
        Vec3 t2=cross(AB, n2);
        t2.normalize(); // necessary ?

        nMean.normalize();

        if ((nMean*PQ) < 0  && !bothSide_computation) // test
        {
            msg_info_when(EMIT_EXTRA_DEBUG_MESSAGE)
                    <<" rejected because of nMean: "<<nMean ;
            return false;
        }

        // compute the angle for the cone to filter contacts using the normal of the triangle situated on the right
        double computedAngleCone = (nMean * t1) * d_coneFactor.getValue();
        if (computedAngleCone<0)
            computedAngleCone=0.0;
        computedAngleCone+=d_angleCone.getValue();

        if (t1*PQ < -computedAngleCone*PQ.norm())
        {
            msg_info_when(EMIT_EXTRA_DEBUG_MESSAGE)
                    <<" rejected because of right triangle normal: "<<n1<<" tang "<< t1 ;
            return false;
        }

        // compute the angle for the cone to filter contacts using the normal of the triangle situated on the left
        computedAngleCone = (nMean * t2) * d_coneFactor.getValue();
        if (computedAngleCone<0)
            computedAngleCone=0.0;
        computedAngleCone+=d_angleCone.getValue();

        if (t2*PQ < -computedAngleCone*PQ.norm())
        {
            msg_info_when(EMIT_EXTRA_DEBUG_MESSAGE)
                <<" rejected because of left triangle normal: "<<n2 ;

            return false;
        }

    }
    else
    {
        n1 = PQ;
        n1.normalize();
        if (fabs(dot(AB,n1)) > d_angleCone.getValue() + 0.0001 )		// dot(AB,n1) should be equal to 0
        {
            // means that proximity was detected with a null determinant
            // in function computeIntersection
            msg_info_when(EMIT_EXTRA_DEBUG_MESSAGE)
                <<"bad case detected  -  abs(dot(AB,n1)) ="<<fabs(dot(AB,n1)) ;
            return false;
        }
    }
    return true;
}

// upstream:1360; 7c18e95d5c5f2839079892c69e7d89a313c79603
bool LocalMinDistance::testValidity(Triangle &t, const Vec3 &PQ) const
{
    const TriangleCollisionModel<sofa::defaulttype::Vec3Types> *tM = t.getCollisionModel();
    const bool bothSide_computation = tM->d_bothSide.getValue();

    if (!d_filterIntersection.getValue() || bothSide_computation)
        return true;

    const Vec3& pt1 = t.p1();
    const Vec3& pt2 = t.p2();
    const Vec3& pt3 = t.p3();

    const Vec3 n = cross(pt2-pt1,pt3-pt1);

    return n * PQ >= 0.0;
}
