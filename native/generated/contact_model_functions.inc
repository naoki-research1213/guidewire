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

// BaseContactLagrangianConstraint.inl:57
void ContactConstraint::addContact(const ContactParams& parameters, Deriv norm, Coord P, Coord Q, Real contactDistance, int m1, int m2, long id, PersistentID localid)
{
    addContact(parameters, norm, P, Q, contactDistance, m1, m2,
            this->getMState2()->read(core::vec_id::read_access::freePosition)->getValue()[m2],
            this->getMState1()->read(core::vec_id::read_access::freePosition)->getValue()[m1],
            id, localid);
}

// BaseContactLagrangianConstraint.inl:66
void ContactConstraint::addContact(const ContactParams& parameters, Deriv norm, Real contactDistance, int m1, int m2, long id, PersistentID localid)
{
    addContact(parameters, norm,
            this->getMState2()->read(core::vec_id::read_access::position)->getValue()[m2],
            this->getMState1()->read(core::vec_id::read_access::position)->getValue()[m1],
            contactDistance, m1, m2,
            this->getMState2()->read(core::vec_id::read_access::freePosition)->getValue()[m2],
            this->getMState1()->read(core::vec_id::read_access::freePosition)->getValue()[m1],
            id, localid);
}

// BaseContactLagrangianConstraint.inl:79
void ContactConstraint::addContact(const ContactParams& parameters, Deriv norm, Coord P, Coord Q, Real contactDistance, int m1, int m2, Coord /*Pfree*/, Coord /*Qfree*/, long id, PersistentID localid)
{
    contacts.resize(contacts.size() + 1);
    Contact &c = contacts.back();

    c.P			= P;
    c.Q			= Q;
    c.m1		= m1;
    c.m2		= m2;
    c.norm		= norm;
    c.t			= Deriv(norm.z(), norm.x(), norm.y());
    c.s			= cross(norm, c.t);
    c.s			= c.s / c.s.norm();
    c.t			= cross((-norm), c.s);
    c.parameters= parameters;
    c.contactId = id;
    c.localId	= localid;
    c.contactDistance = contactDistance;
}

// BaseContactLagrangianConstraint.inl:102
void ContactConstraint::buildConstraintMatrix(const core::ConstraintParams *, DataMatrixDeriv &c1_d, DataMatrixDeriv &c2_d, unsigned int &contactId
        , const DataVecCoord &, const DataVecCoord &)
{
    assert(this->mstate1);
    assert(this->mstate2);

    if (this->mstate1 == this->mstate2)
    {
        helper::WriteAccessor<DataMatrixDeriv> c1 = c1_d;


        for (unsigned int i = 0; i < contacts.size(); i++)
        {
            Contact& c = contacts[i];

            c.id = contactId++;

            MatrixDerivRowIterator c1_it = c1->writeLine(c.id);

            c1_it.addCol(c.m1, -c.norm);
            c1_it.addCol(c.m2, c.norm);

            if (c.parameters.hasTangentialComponent())
            {
                c1_it = c1->writeLine(c.id + 1);
                c1_it.setCol(c.m1, -c.t);
                c1_it.setCol(c.m2, c.t);

                c1_it = c1->writeLine(c.id + 2);
                c1_it.setCol(c.m1, -c.s);
                c1_it.setCol(c.m2, c.s);

                contactId += 2;
            }
        }

    }
    else
    {
        helper::WriteAccessor<DataMatrixDeriv> c1 = c1_d;
        helper::WriteAccessor<DataMatrixDeriv> c2 = c2_d;

        for (unsigned int i = 0; i < contacts.size(); i++)
        {
            Contact& c = contacts[i];

            c.id = contactId++;

            MatrixDerivRowIterator c1_it = c1->writeLine(c.id);
            c1_it.addCol(c.m1, -c.norm);

            MatrixDerivRowIterator c2_it = c2->writeLine(c.id);
            c2_it.addCol(c.m2, c.norm);

            if (c.parameters.hasTangentialComponent())
            {
                c1_it = c1->writeLine(c.id + 1);
                c1_it.setCol(c.m1, -c.t);

                c1_it = c1->writeLine(c.id + 2);
                c1_it.setCol(c.m1, -c.s);

                c2_it = c2->writeLine(c.id + 1);
                c2_it.setCol(c.m2, c.t);

                c2_it = c2->writeLine(c.id + 2);
                c2_it.setCol(c.m2, c.s);

                contactId += 2;
            }
        }

    }
}

// BaseContactLagrangianConstraint.inl:179
void ContactConstraint::getPositionViolation(linearalgebra::BaseVector *v)
{
    const VecCoord &PfreeVec = this->getMState2()->read(core::vec_id::read_access::freePosition)->getValue();
    const VecCoord &QfreeVec = this->getMState1()->read(core::vec_id::read_access::freePosition)->getValue();

    Real dfree = (Real)0.0;
    Real dfree_t = (Real)0.0;
    Real dfree_s = (Real)0.0;

    const unsigned int cSize = contacts.size();

    for (unsigned int i = 0; i < cSize; i++)
    {
        const Contact& c = contacts[i];

        // Compute dfree, dfree_t and d_free_s

        const Coord &Pfree =  PfreeVec[c.m2];
        const Coord &Qfree =  QfreeVec[c.m1];

        const Coord PPfree = Pfree - c.P;
        const Coord QQfree = Qfree - c.Q;

        const Real ref_dist = PPfree.norm() + QQfree.norm();

        dfree = dot(Pfree - Qfree, c.norm) - c.contactDistance;
        const Real delta = dot(c.P - c.Q, c.norm) - c.contactDistance;

        if ((helper::rabs(delta) < 0.00001 * ref_dist) && (helper::rabs(dfree) < 0.00001 * ref_dist))
        {
            dfree_t = dot(PPfree, c.t) - dot(QQfree, c.t);
            dfree_s = dot(PPfree, c.s) - dot(QQfree, c.s);
        }
        else if (helper::rabs(delta - dfree) > 0.001 * delta)
        {
            const Real dt = delta / (delta - dfree);

            if (dt > 0.0 && dt < 1.0)
            {
                const Coord Pt		= c.P * (1 - dt) + Pfree * dt;
                const Coord Qt		= c.Q * (1 - dt) + Qfree * dt;
                const Coord PtPfree = Pfree - Pt;
                const Coord QtQfree = Qfree - Qt;

                dfree_t = dot(PtPfree, c.t) - dot(QtQfree, c.t);
                dfree_s = dot(PtPfree, c.s) - dot(QtQfree, c.s);
            }
            else if (dfree < 0.0)
            {
                dfree_t = dot(PPfree, c.t) - dot(QQfree, c.t);
                dfree_s = dot(PPfree, c.s) - dot(QQfree, c.s);
            }
            else
            {
                dfree_t = 0;
                dfree_s = 0;
            }
        }
        else
        {
            dfree_t = dot(PPfree, c.t) - dot(QQfree, c.t);
            dfree_s = dot(PPfree, c.s) - dot(QQfree, c.s);
        }

        // Sets dfree in global violation vector

        v->set(c.id, dfree);

        c.dfree = dfree; // PJ : For isActive() method. Don't know if it's still usefull.

        if (c.parameters.hasTangentialComponent())
        {
            v->set(c.id + 1, dfree_t);
            v->set(c.id + 2, dfree_s);
        }
    }
}

// BaseContactLagrangianConstraint.inl:315
void ContactConstraint::getConstraintInfo(const core::ConstraintParams*, VecConstraintBlockInfo& blocks, VecPersistentID& ids)
{
    if (contacts.empty()) return;
    const bool friction = (contacts[0].parameters.mu > 0.0); /// @todo: can there be both friction-less and friction contacts in the same BaseContactLagrangianConstraint ???
    ConstraintBlockInfo info;
    info.parent = this;
    info.const0 = contacts[0].id;
    info.nbLines = friction ? 3 : 1;
    info.hasId = true;
    info.offsetId = ids.size();
    info.nbGroups = contacts.size();

    for (unsigned int i=0; i<contacts.size(); i++)
    {
        Contact& c = contacts[i];
        ids.push_back( yetIntegrated ? c.contactId : -c.contactId);
    }


    yetIntegrated = true;

    blocks.push_back(info);
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
#pragma once

// BaseUnilateralContactResponse.inl:88
void ContactResponse::setDetectionOutputs(OutputVector* o)
{
    TOutputVector& outputs = *static_cast<TOutputVector*>(o);
    // We need to remove duplicate contacts
    constexpr double minDist2 = 0.00000001f;

    contacts.clear();

    if (model1->getContactStiffness(0) == 0 || model2->getContactStiffness(0) == 0)
    {
        msg_error() << "Disabled BaseUnilateralContactResponse with " << (outputs.size()) << " collision points.";
        return;
    }

    contacts.reserve(outputs.size());

    const int SIZE = outputs.size();

    // the following procedure cancels the duplicated detection outputs
    for (int cpt=0; cpt<SIZE; cpt++)
    {
        sofa::core::collision::DetectionOutput* detectionOutput = &outputs[cpt];

        bool found = false;
        for (unsigned int i=0; i<contacts.size() && !found; i++)
        {
            const sofa::core::collision::DetectionOutput* p = contacts[i];
            if ((detectionOutput->point[0]-p->point[0]).norm2()+(detectionOutput->point[1]-p->point[1]).norm2() < minDist2)
                found = true;
        }

        if (!found)
            contacts.push_back(detectionOutput);
    }

    // DUPLICATED CONTACTS FOUND
    msg_info_when(contacts.size()<outputs.size()) << "Removed " << (outputs.size()-contacts.size()) <<" / " << outputs.size() << " collision points." << msgendl;
}

// BaseUnilateralContactResponse.inl:129
void ContactResponse::activateMappers()
{
    if (!m_constraint)
    {
        // Get the mechanical model from mapper1 to fill the constraint vector
        MechanicalState1* mmodel1 = mapper1.createMapping(getName().c_str());
        // Get the mechanical model from mapper2 to fill the constraints vector
        MechanicalState2* mmodel2;
        if (selfCollision)
        {
            mmodel2 = mmodel1;
        }
        else
        {
            mmodel2 = mapper2.createMapping(getName().c_str());
        }
        setupConstraint(mmodel1,mmodel2);
        m_constraint->setName( getName() );
        setInteractionTags(mmodel1, mmodel2);
        m_constraint->setCustomTolerance(d_tol.getValue() );
    }

    int size = contacts.size();
    m_constraint->clear(size);
    if (selfCollision)
        mapper1.resize(2*size);
    else
    {
        mapper1.resize(size);
        mapper2.resize(size);
    }
    int i = 0;
    const double d0 = intersectionMethod->getContactDistance() + model1->getContactDistance() + model2->getContactDistance(); // - 0.001;

    mappedContacts.resize(contacts.size());
    for (std::vector<sofa::core::collision::DetectionOutput*>::const_iterator it = contacts.begin(); it!=contacts.end(); it++, i++)
    {
        sofa::core::collision::DetectionOutput* o = *it;
        CollisionElement1 elem1(o->elem.first);
        CollisionElement2 elem2(o->elem.second);
        int index1 = elem1.getIndex();
        int index2 = elem2.getIndex();

        typename DataTypes1::Real r1 = 0.;
        typename DataTypes2::Real r2 = 0.;

        // Create mapping for first point
        index1 = mapper1.addPointB(o->point[0], index1, r1);
        // Create mapping for second point
        if (selfCollision)
        {
            index2 = mapper1.addPointB(o->point[1], index2, r2);
        }
        else
        {
            index2 = mapper2.addPointB(o->point[1], index2, r2);
        }
        const double distance = d0 + r1 + r2;

        mappedContacts[i].first.first = index1;
        mappedContacts[i].first.second = index2;
        mappedContacts[i].second = distance;
    }

    // Update mappings
    mapper1.update();
    mapper1.updateXfree();
    if (!selfCollision) mapper2.update();
    if (!selfCollision) mapper2.updateXfree();
}

// BaseUnilateralContactResponse.inl:201
void ContactResponse::createResponse(core::objectmodel::BaseContext* group)
{

    activateMappers();

    if (m_constraint)
    {
        int i=0;
        for (std::vector<sofa::core::collision::DetectionOutput*>::const_iterator it = contacts.begin(); it!=contacts.end(); it++, i++)
        {
            const sofa::core::collision::DetectionOutput* o = *it;
            const int index1 = mappedContacts[i].first.first;
            const int index2 = mappedContacts[i].first.second;
            const double distance = mappedContacts[i].second;

            // Polynome de Cantor de NxN sur N bijectif f(x,y)=((x+y)^2+3x+y)/2
            const long index = cantorPolynomia(o->id /*cantorPolynomia(index1, index2)*/,id);

            const ConstraintParameters params = getParameterFromDatas();

            // Add contact in unilateral constraint
            m_constraint->addContact(params, o->normal, distance, index1, index2, index, o->id);
        }

        if (parent!=nullptr)
        {
            parent->removeObject(this);
            parent->removeObject(m_constraint);
        }

        parent = group;
        if (parent!=nullptr)
        {
            parent->addObject(this);
            parent->addObject(m_constraint);
        }
    }
}
