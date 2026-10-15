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

namespace sofa::helper {

// upstream:549; 7c18e95d5c5f2839079892c69e7d89a313c79603
void LocalBlock33::compute(SReal& w11, SReal& w12, SReal& w13, SReal& w22, SReal& w23, SReal& w33)
{
    w[0]=w11; w[1]=w12; w[2] = w13; w[3]=w22; w[4]=w23; w[5]=w33;
    det = w11*w22*w33-w11*w23*w23-w12*w12*w33+2*w12*w13*w23-w13*w13*w22;
    wInv[0] = (w22*w33-w23*w23)/det;
    wInv[1] = -(w12*w33-w13*w23)/det;
    wInv[2] = (w12*w23-w13*w22)/det;
    wInv[3] = (w11*w33-w13*w13)/det;
    wInv[4] = -(w11*w23-w12*w13)/det;
    wInv[5] = (w11*w22-w12*w12)/det;
    computed=true;
}

// upstream:612; 7c18e95d5c5f2839079892c69e7d89a313c79603
void LocalBlock33::New_GS_State(SReal &mu, SReal &dn, SReal &dt, SReal &ds, SReal &fn, SReal &ft, SReal &fs)
{

    SReal d[3];
    SReal normFt;
    f_1[0]=fn; f_1[1]=ft; f_1[2]=fs;

    // evaluation of the current normal position
    d[0] = dn;
    // evaluation of the new contact force
    fn -= d[0]/w[0];

    if (fn <= 0)
    {
        fn=0; ft=0; fs=0;
        // if the force was previously not null -> update the state
        if (f_1[0]>0)
        {
            SReal df[3];
            df[0] = fn-f_1[0];  df[1] = ft-f_1[1];  df[2] = fs-f_1[2];

            dn += w[0]*df[0] + w[1]*df[1] + w[2]*df[2];
            dt += w[1]*df[0] + w[3]*df[1] + w[4]*df[2];
            ds += w[2]*df[0] + w[4]*df[1] + w[5]*df[2];
        }
        return;
    }


    // evaluation of the current tangent positions
    d[1] = w[1]*(fn-f_1[0]) + dt;
    d[2] = w[2]*(fn-f_1[0]) + ds;

    // envaluation of the new fricton forces
    ft -= 2*d[1]/(w[3]+w[5]);
    fs -= 2*d[2]/(w[3]+w[5]);

    normFt=sqrt(ft*ft+fs*fs);

    if (normFt > mu*fn)
    {
        ft *=mu*fn/normFt;
        fs *=mu*fn/normFt;
    }

    SReal df[3];
    df[0] = fn-f_1[0];  df[1] = ft-f_1[1];  df[2] = fs-f_1[2];

    dn += w[0]*df[0] + w[1]*df[1] + w[2]*df[2];
    dt += w[1]*df[0] + w[3]*df[1] + w[4]*df[2];
    ds += w[2]*df[0] + w[4]*df[1] + w[5]*df[2];



}

}
