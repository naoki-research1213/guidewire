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

// Original RigidCoord<3,real>::operator+=(Deriv), unchanged body.
void integratePose(Pose& pose,const RigidDerivativeView& dg)
{
using real=double;auto& center=pose.p;auto& orientation=pose.q;
{
        // R3 x SO(3) exponential integration
        center += dg.getVCenter();

        const Vec3 omega = dg.getVOrientation() / 2;
        const real theta = omega.norm();

        constexpr const real epsilon = std::numeric_limits<real>::epsilon();

        if (theta < epsilon) {
            // fallback to gnomonic projection
            Quat exp(omega[0], omega[1], omega[2], 1);
            exp.normalize();
            orientation = exp * orientation;
        }
        else {
            // expontential
            const real sinc = std::sin(theta) / theta;
            const Quat exp(sinc * omega[0],
                sinc * omega[1],
                sinc * omega[2],
                std::cos(theta));
            orientation = exp * orientation;
        }

    }
}
