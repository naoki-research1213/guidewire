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


// upstream:391; 7c18e95d5c5f2839079892c69e7d89a313c79603
int ContactSolver::nlcp_gaussseidel_unbuilt(SReal *dfree, SReal *f, std::vector<SReal>* residuals)
{
    if(!_numConstraints)
        return 0;

    auto _mu = d_mu.getValue();
    if(_mu==0.0)
    {
        msg_error() << "frictionless case with unbuilt nlcp is not implemented";
        return 0;
    }

    if (_numConstraints%3 != 0)
    {
        msg_error() << "dim should be dividable by 3 in nlcp_gaussseidel";
        return 0;
    }

    sofa::helper::advancedtimer::stepBegin("build_constraints");

    int numContacts =  _numConstraints/3;

    int it,c1;

    // data for iterative procedure
    SReal _tol = d_tol.getValue();
    int _maxIt = d_maxIt.getValue();

    /// each constraintCorrection has an internal force vector that is set to "0"

    // indirection of the sequence of contact
    std::list<unsigned int> contact_sequence;

    for (unsigned int c=0; c< _numConstraints; c++)
    {
        contact_sequence.push_back(c);
    }


    for (const auto& cc : l_constraintCorrections)
    {
        cc->resetForUnbuiltResolution(f, contact_sequence);
    }

    sofa::helper::advancedtimer::stepEnd("build_constraints");

    auto linkConstraintTimer = std::make_unique<sofa::helper::ScopedAdvancedTimer>("link_constraints");

    //////// Important component if the LCP is not build :
    // for each contact, the pair of constraint correction that is involved with the contact is memorized
    _cclist_elem1.clear();
    _cclist_elem2.clear();
    std::vector<int> missingConstraintCorrectionContacts;
    for (c1=0; c1<numContacts; c1++)
    {
        bool elem1 = false;
        bool elem2 = false;

        for (const auto& cc : l_constraintCorrections)
        {
            if(cc->hasConstraintNumber(3*c1))
            {
                if(elem1)
                {
                    _cclist_elem2.push_back(cc);
                    elem2=true;
                }
                else
                {
                    _cclist_elem1.push_back(cc);
                    elem1=true;
                }

            }
        }
        if (!elem1)
        {
            _cclist_elem1.push_back(nullptr);
            missingConstraintCorrectionContacts.push_back(c1);
        }
        if(!elem2)
            _cclist_elem2.push_back(nullptr);

    }

    if (!missingConstraintCorrectionContacts.empty())
    {
        std::stringstream ss;
        for (const auto c : missingConstraintCorrectionContacts)
        {
            ss << c << ' ';
        }
        msg_error() << "The following contacts do not have an associated constraint correction component: " << ss.str();
        if (missingConstraintCorrectionContacts.size() == _cclist_elem1.size())
        {
            msg_error() << "None of the contacts has an associated constraint correction component: constraint correction is aborted";
            return 0;
        }
    }

    // memory allocation of vector d
    unbuilt_d.resize(_numConstraints);
    SReal *d = &(unbuilt_d[0]);

    linkConstraintTimer.reset();
    auto buildDiagonalTimer = std::make_unique<sofa::helper::ScopedAdvancedTimer>("build_diagonal");

    //////////////
    // Beginning of iterative computations
    //////////////


    /////////// the 3x3 diagonal block matrix is built:
    /////////// for each contact, the pair of constraintcorrection is called to add the contribution
    for (c1=0; c1<numContacts; c1++)
    {
        dmsg_info() << "contact " << c1 << " cclist_elem1 : " << _cclist_elem1[c1]->getName();

        // compliance of object1
        if (_cclist_elem1[c1] != nullptr)
        {
            _cclist_elem1[c1]->getBlockDiagonalCompliance(&_Wdiag, 3 * c1, 3 * c1 + 2);
        }
        // compliance of object2 (if object2 exists)
        if(_cclist_elem2[c1] != nullptr)
        {
            _cclist_elem2[c1]->getBlockDiagonalCompliance(&_Wdiag, 3*c1, 3*c1+2);


           dmsg_info() <<"  _cclist_elem2 : "<<_cclist_elem2[c1]->getName();
        }
        dmsg_info() <<" "<<msgendl;
    }

    addRegularization(_Wdiag);


    // allocation of the inverted system 3x3
    // TODO: evaluate the cost of this step : it can be avoied by directly feeding W33 in constraint correction
    unbuilt_W33.clear();
    unbuilt_W33.resize(numContacts);
    helper::LocalBlock33 *W33 = &(unbuilt_W33[0]); //new helper::LocalBlock33[numContacts];
    for (c1=0; c1<numContacts; c1++)
    {
        SReal w[6];
        w[0] = _Wdiag.element(3*c1  , 3*c1  );
        w[1] = _Wdiag.element(3*c1  , 3*c1+1);
        w[2] = _Wdiag.element(3*c1  , 3*c1+2);
        w[3] = _Wdiag.element(3*c1+1, 3*c1+1);
        w[4] = _Wdiag.element(3*c1+1, 3*c1+2);
        w[5] = _Wdiag.element(3*c1+2, 3*c1+2);
        W33[c1].compute(w[0], w[1] , w[2], w[3], w[4] , w[5]);
    }


    dmsg_info() <<" Compliance In constraint Space : \n W ="<<(* _W)<<msgendl
                <<"getBlockDiagonalCompliance   \n Wdiag = "<< _Wdiag ;

    buildDiagonalTimer.reset();

    SCOPED_TIMER_VARNAME(gaussSeidelTimer, "GAUSS_SEIDEL");

    SReal error = 0;
    SReal dn, dt, ds, fn, ft, fs, fn0;

    for (it=0; it<_maxIt; it++)
    {
        std::list<unsigned int>::iterator it_c ;
        error =0;

        //constraints are treated 3x3 (friction contact)
        for (it_c = contact_sequence.begin(); it_c != contact_sequence.end() ; std::advance(it_c, 3) )
        {
            int constraint = *it_c;
            c1 = constraint/3;

            // compute the current violation :

            // violation when no contact force
            d[3*c1]=dfree[3*c1]; d[3*c1+1]=dfree[3*c1+1]; d[3*c1+2]=dfree[3*c1+2];


            // set current force in fn, ft, fs
            fn0=fn=f[3*c1]; ft=f[3*c1+1]; fs=f[3*c1+2];

            // displacement of object1 due to contact force
            if(_cclist_elem1[c1] != nullptr)
                _cclist_elem1[c1]->addConstraintDisplacement(d, 3*c1, 3*c1+2);

            // displacement of object2 due to contact force (if object2 exists)
            if(_cclist_elem2[c1] != nullptr)
                _cclist_elem2[c1]->addConstraintDisplacement(d, 3*c1, 3*c1+2);


            // set displacement in dn, dt, ds
            dn=d[3*c1]; dt=d[3*c1+1]; ds=d[3*c1+2];

            // compute a new state for stick/slip
            /// ATTENTION  NOUVEAU GS_STATE : maintenant dn, dt et ds inclue les forces fn, ft, fs
            W33[c1].New_GS_State(_mu,dn,dt,ds,fn,ft,fs);

            // evaluate an error (based on displacement)
            error += helper::absError(dn,dt,ds,d[3*c1],d[3*c1+1],d[3*c1+2]);

            bool update;
            if (fn0 == 0.0 && fn == 0.0)
                update=false;               // the contact is not active and was not active in the previous step
            else
                update=true;

            // set the new force :
            // compute the Delta of contact forces:
            f[3*c1  ] = fn - f[3*c1  ];
            f[3*c1+1] = ft - f[3*c1+1];
            f[3*c1+2] = fs - f[3*c1+2];

            ///////// verifier si Delta force vaut 0 => pas la peine d'ajouter la force

            // set Delta force on object 1 for evaluating the following displacement

            if(update)
            {
                if(_cclist_elem1[c1] != nullptr)
                    _cclist_elem1[c1]->setConstraintDForce(f, 3*c1, 3*c1+2, update);

                // set Delta force on object2 (if object2 exists)
                if(_cclist_elem2[c1] != nullptr)
                    _cclist_elem2[c1]->setConstraintDForce(f, 3*c1, 3*c1+2, update);
            }

            // set force on the contact force vector
            helper::set3Dof(f,c1,fn,ft,fs);
        }

        residuals->push_back(error);


        if (error < _tol*(numContacts+1))
        {
            msg_info() << "convergence after "<<it<<" iterations - error " << error;

            sofa::helper::AdvancedTimer::valSet("GS iterations", it+1);
            return 1;
        }

    }

    sofa::helper::AdvancedTimer::valSet("GS iterations", it);

    dmsg_warning() << "No convergence in unbuilt nlcp gaussseidel function : error ="
                <<error <<" after "<< it<<" iterations";

    return 0;
}
