/*
 *   This file is part of the OpenPhase (R) software library.
 *
 *   Copyright (c) 2009-2022 Ruhr-Universitaet Bochum,
 *                 Universitaetsstrasse 150, D-44801 Bochum, Germany
 *             AND 2018-2022 OpenPhase Solutions GmbH,
 *                 Universitaetsstrasse 136, D-44799 Bochum, Germany.
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   File created :   2011
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Philipp Engels;
 *                         Muhammad Adil Ali; Hesham Salama
 *
 */

#include "Tools.h"
#include "Info.h"
#include "PhaseField.h"
#include "Orientations.h"
#include "Crystallography.h"


namespace openphase
{

using namespace std;
/*
void Tools::WriteAbaqusInput(const PhaseField& Phi, const Orientations& OR)
{
    const int Nx = Phi.Nx;
    const int Ny = Phi.Ny;
    const int Nz = Phi.Nz;
    const double dx = Phi.dx;

    Storage3D<dVector3, 0> Coordinates;                                         // Stores gridpoint coordinates.
    Storage3D<int, 0> ElementNumber;                                            // Stores ElementNumber.
    Storage3D<int, 0> ElementSet;                                               // Stores ElementSet.

    Coordinates.Allocate(Nx+1, Ny+1, Nz+1, 0);
    ElementNumber.Allocate(Nx+2, Ny+2, Nz+2, 0);
    ElementSet.Allocate(Nx+2, Ny+2, Nz+2, 0);

    string GeometryFileName = "geometry.inp";
    string MaterialFileName = "materialinput.inp";
    string GrainFileName = "graindata.inp";

    Info::WriteLine();
    Info::WriteStandard("OpenPhase to Abaqus Parser", "Started");
    Info::WriteStandard("Abaqus geometry file ", GeometryFileName);
    Info::WriteStandard("Abaqus material file ", MaterialFileName);

    // Analyze PhaseField
    int indGrains = 0;
    for(unsigned int n = 0; n < Phi.FieldsStatistics.size(); n++)
    {
        if(Phi.FieldsStatistics[n].Exist)
        {
            indGrains += 1;
        }
    }
    Info::WriteStandard("Individual grains found", std::to_string(indGrains));

    // Determine and write node coordinates
    for (int k =    0; k < Nz+1; k++)
    for (int j =    0; j < Ny+1; j++)
    for (int i =    0; i < Nx+1; i++)
    {
        Coordinates(i,j,k)[0]= i*dx;
        Coordinates(i,j,k)[1]= j*dx;
        Coordinates(i,j,k)[2]= k*dx;
    }

    stringstream outbuffer;

    int elementnumber = 1;
    outbuffer << "*Part, name=PART-1" << endl;
    outbuffer << "*Node, nset=nall" << endl;
    for (int k =    0; k < Nz+1; k++)
    for (int j =    0; j < Ny+1; j++)
    for (int i =    0; i < Nx+1; i++)
    {
        outbuffer << elementnumber << ", "
                  << Coordinates(i,j,k)[0] << ", "
                  << Coordinates(i,j,k)[1] << ", "
                  << Coordinates(i,j,k)[2] << endl;
        elementnumber +=1;
    }

    // Determine and write elements

    Info::WriteStandard("Number of elements (C3D8)", std::to_string(Nx*Ny*Nz));
    outbuffer << "*Element, type=C3D8" << endl;

    elementnumber = 1;
    int elementnodes[8];

    for (int k =    1; k < Nz+1; k++)
    for (int j =    1; j < Ny+1; j++)
    for (int i =    1; i < Nx+1; i++)
    {
        int ii = i;
        int jj = j;
        int kk = k;

        elementnodes[0] = ii + (Nx+1)*(jj-1) + (Nx+1)*(Ny+1)*(kk-1);
        elementnodes[1] = ii + (Nx+1)*(jj-1) + (Nx+1)*(Ny+1)*(kk-1) + 1;
        elementnodes[2] = ii + (Nx+1)*jj     + (Nx+1)*(Ny+1)*(kk-1) + 1;
        elementnodes[3] = ii + (Nx+1)*jj     + (Nx+1)*(Ny+1)*(kk-1);
        elementnodes[4] = ii + (Nx+1)*(jj-1) + (Nx+1)*(Ny+1)*kk;
        elementnodes[5] = ii + (Nx+1)*(jj-1) + (Nx+1)*(Ny+1)*kk + 1;
        elementnodes[6] = ii + (Nx+1)*jj     + (Nx+1)*(Ny+1)*kk + 1;
        elementnodes[7] = ii + (Nx+1)*jj     + (Nx+1)*(Ny+1)*kk;

        ElementNumber(i,j,k) = elementnumber;

        outbuffer << elementnumber << ", ";
        for (int n = 0; n < 7; n++)
        {
            outbuffer << elementnodes[n] << ", ";
        }
        outbuffer << elementnodes[7] << endl;

        // Determine element sets
        // In the interface, the majority phase is used.

        if(Phi.Interface(i-1,j-1,k-1))
        {
            int maxindex = 0.0;
            double maxvalue = 0.0;
            for(auto alpha = Phi.Fields(i-1,j-1,k-1).cbegin();
                    alpha < Phi.Fields(i-1,j-1,k-1).cend(); alpha++)
            {
                if (alpha -> value > maxvalue)
                {
                    maxvalue = alpha -> value;
                    maxindex = alpha -> index;
                }
            }
            ElementSet(i,j,k) = maxindex;
        }
        else
        {
            ElementSet(i,j,k) = Phi.FieldIndex(i-1,j-1,k-1);
        }
        elementnumber += 1;
    }
    outbuffer << "*End Part" << endl;
    outbuffer << "*Assembly, name=Assembly" << endl;
    outbuffer << "*Instance, name=PART-1-1, part=PART-1" << endl;

    // Write element sets
    for (unsigned int ph = 0; ph < Phi.FieldsStatistics.size(); ph++)
    if(Phi.FieldsStatistics[ph].Exist)
    {
        int rowcount = 1;
        outbuffer << "*Elset, elset=Grain" << ph + 1 << endl;
        for (int k =    1; k < Nz+1; k++)
        for (int j =    1; j < Ny+1; j++)
        for (int i =    1; i < Nx+1; i++)
        {
            if((int)ph == ElementSet(i,j,k))
            {
                outbuffer << ElementNumber(i,j,k) << ", ";
                if (rowcount == 15) // Size limitation in Abaqus input files
                {
                    outbuffer << endl;
                    rowcount = 1;
                }
                else
                {
                    rowcount +=1;
                }
            }
        }
        outbuffer << endl << "*Solid Section, elset=Grain" << ph + 1 << ", material=Grain" << ph + 1 << "\n";
    }

    // WriteInterface
    if(Phi.Eta > 0)
    {
        int rowcount = 1;
        outbuffer << "*Elset, elset=InterfaceElements" << endl;
        for (int k =    1; k < Nz+1; k++)
        for (int j =    1; j < Ny+1; j++)
        for (int i =    1; i < Nx+1; i++)
        {
            if(Phi.Interface(i-1,j-1,k-1))
            {
                outbuffer << ElementNumber(i,j,k) << ", ";
                if (rowcount == 15) // Size limitation in Abaqus input files
                {
                    outbuffer << endl;
                    rowcount = 1;
                }
                else
                {
                    rowcount +=1;
                }
            }
        }
        outbuffer << endl;
        for (unsigned int ph = 0; ph < Phi.FieldsStatistics.size(); ph++)
        if(Phi.FieldsStatistics[ph].Exist)
        {
            int rowcount = 1;
            outbuffer << "*Elset, elset=Interface_" << ph << endl;
            for (int k =    1; k < Nz+1; k++)
            for (int j =    1; j < Ny+1; j++)
            for (int i =    1; i < Nx+1; i++)
            {
                if(Phi.Interface(i-1,j-1,k-1) and Phi.Fields(i-1,j-1,k-1).get(ph) >= 0.5)
                {
                    outbuffer << ElementNumber(i,j,k) << ", ";
                    if (rowcount == 15) // Size limitation in Abaqus input files
                    {
                        outbuffer << endl;
                        rowcount = 1;
                    }
                    else
                    {
                        rowcount +=1;
                    }
                }
            }
            outbuffer << endl;
        }
    } // end interface existing

    int rowcount = 1;
    outbuffer << "*Elset, elset=MIDDLEELEMENTS" << endl;
    for (int i =    1; i < Nz+1; i++)
    {
        outbuffer << ElementNumber(Nx/2,Ny/2+1, i) << ", ";
        if (rowcount == 15) // Size limitation in Abaqus input files
        {
            outbuffer << endl;
            rowcount = 1;
        }
        else
        {
            rowcount +=1;
        }
    }
    outbuffer << endl;

    outbuffer << "*End Instance" << endl;
    outbuffer << "*End Assembly" << endl;
    outbuffer << "*Include, input=materialinput.inp" << endl;
    //outbuffer << "*Include, input=boundary.inp" << endl;

    ofstream Geometry_file(GeometryFileName.c_str());
    Geometry_file << outbuffer.rdbuf();
    Geometry_file.close();

    outbuffer.str(std::string());                                               // wipe
    for (unsigned int ph = 0; ph < Phi.FieldsStatistics.size(); ph++)
    if(Phi.FieldsStatistics[ph].Exist)
    {
        outbuffer << "*Material, name=Grain" << ph + 1 << endl
                  << "*Depvar" << endl
                  << "     53" << endl
                  << "*User Material, constants = 2" << endl
                  << ph + 1  << " , 2" << endl;
    }

    ofstream Material_file(MaterialFileName.c_str());
    Material_file << outbuffer.rdbuf();
    Material_file.close();

    Info::WriteStandard("Abaqus grain data file ", GrainFileName);

    // Writes seperate file with Euler angles
    ofstream Grain_file(GrainFileName.c_str());
    outbuffer.str(std::string());                                               // wipe

    for (unsigned int ph = 0; ph < Phi.FieldsStatistics.size(); ph++)
    if(Phi.FieldsStatistics[ph].Exist)
    {
        Angles tempOrient = OR.GrainEulerAngles[ph].get_degree();
        outbuffer << "Grain : " << ph + 1 << " : " <<
                           tempOrient.Q[0] <<  " : " <<
                           tempOrient.Q[1] <<  " : " <<
                           tempOrient.Q[2] << endl;
    }

    Grain_file << outbuffer.rdbuf();
    Grain_file.close();

    Info::WriteLine();
}
*/
void Tools::Eigensystem(dMatrix3x3& M, dMatrix3x3& Eigenvectors, dMatrix3x3& Eigenvalues)
{
#ifdef DEBUG
    ///M MUST BE _SYMMETRIC_
    bool symmetric = true;
    for(int i = 0; i < 3; i++)
    for(int j = i; j < 3; j++)
    {
        if(fabs(M(i,j) - M(j,i)) > DBL_EPSILON)
            symmetric = false;
    }
    if(not symmetric)
    {
        Info::WriteExit("M:\n"
                        + M.print() + "is not symmetric!", "Tools", "Eigensystem()");
        exit(13);
    }
#endif
    int max_iterations=1000000;
    double delta = 1.0e-6;
    double Check = 1.0/delta;
    dMatrix3x3 M2,M4,M8,M16,
    //M32,M64,M128, M256, M512,
    MFinal;
    M2=M*M;
    M4=M2*M2;
    M8=M4*M4;
    M16=M8*M8;
    //M32=M16*M16;
    //M64=M32*M32;
    //M128=M64*M64;
    //M256=M128*M128;
    //M512=M256*M256;
    MFinal=M;

    dVector3 u0;
    u0.set_to_unitX();
    dVector3 tmp;
    tmp = M*u0;                                                                 //tmp = M.u0
    tmp.normalize();
    dVector3 u = tmp;

    dVector3 v = u-u0;                                                          //{u[0] - u0[0], u[1] - u0[1], u[2] - u0[2]};
    if(v.length() < delta)
        v.setY(1.0);
    v.normalize();

    int n = 1;
    while(Check >= delta && fabs(2.0 - Check) > delta)
    {
        tmp.set_to_zero();
        tmp -= u;

        u = MFinal*u;                                                           //u_i=M^32 . u_(i-1), high power of M accelerates convergence
        u.normalize();

        tmp += u;

        Check = tmp.length();                                                   //TODO: to avoid sqrt() calculation it's possible to use smth. like a Lengthsqr() fct!

        n++;
        if(n >= max_iterations and MFinal != M16)
        {
            n=0;
            MFinal = M16;
        }

        if(n >= max_iterations and MFinal == M16)
        {
            stringstream msg;
            msg << "Eigensystem() could not converge! M:\n " << M.print()
                << "Last suggested solution for the first eigenvector: "
                << u.print() << "; last correction: " << Check <<" > "
                << delta << endl;
            Info::WriteExit( msg.str(), "Tools", "Eigensystem()");
            exit(13);
        }
    }
//    cout<<n-1<<" iterations, Check: "<<Check<<endl;

    double v_cdot_u = u*v;
    v -= (u*v_cdot_u);
    v.normalize();

    n = 1;
    Check = 1.0/delta;
    MFinal = M;
    while(Check >= delta && fabs(2.0 - Check) > delta)
    {
        tmp.set_to_zero();
        tmp -= v;
        v = MFinal*v;                                                           //v_i=M^32 . v_(i-1) ; it's better not to use to high powers of M here (-;

        v_cdot_u = u*v;
        v -= (u*v_cdot_u);
        v.normalize();

        tmp += v;

        Check = tmp.length();                                                   //TODO: to avoid sqrt() calculation it's possible to use smth. like a Lengthsqr() fct!

        n++;
        if(n >= max_iterations and MFinal != M16)
        {
            n=0;
            MFinal = M16;
        }
        if(n >= max_iterations and MFinal == M16)
        {
            stringstream msg;
            msg << "Eigensystem() could not converge! M:\n " << M.print()
                << "Last suggested solution for the second eigenvector: "
                << v.print() << "; last correction: " << Check <<" > "
                << delta << endl;
            Info::WriteExit(msg.str(), "Tools", "Eigensystem()");
            exit(13);
        }
    }
//    cout<<n-1<<" iterations, Check: "<<Check<<endl;

    dVector3 w;
    w = u.cross(v);

    dMatrix3x3 Eigenvectors_tr;                                                 //Row Eigenvectors
    Eigenvectors_tr(0,0) = u.getX();
    Eigenvectors_tr(0,1) = u.getY();
    Eigenvectors_tr(0,2) = u.getZ();
    Eigenvectors_tr(1,0) = v.getX();
    Eigenvectors_tr(1,1) = v.getY();
    Eigenvectors_tr(1,2) = v.getZ();
    Eigenvectors_tr(2,0) = w.getX();
    Eigenvectors_tr(2,1) = w.getY();
    Eigenvectors_tr(2,2) = w.getZ();

    Eigenvectors = Eigenvectors_tr.transposed();                                //Column Eigenvectors
    Eigenvalues = Eigenvectors_tr * M * Eigenvectors;
    for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
        if(i != j)
        {
            Eigenvalues(i,j) = 0.0;
        }
    }
}

void Tools::Eigensystem2(dMatrix3x3& M, dMatrix3x3& Eigenvectors, dMatrix3x3& Eigenvalues, bool calculateEigenVectors)
{
#ifdef DEBUG
    ///M MUST BE _SYMMETRIC_
    bool symmetric = true;
    for(int i = 0; i < 3; i++)
    for(int j = i; j < 3; j++)
    {
        if(fabs(M(i,j) - M(j,i)) > DBL_EPSILON)
            symmetric = false;
    }
    if(not symmetric)
    {
        Info::WriteExit("M:\n"
                        + M.print() + "is not symmetric!", "Tools", "Eigensystem2()");
        exit(13);
    }
#endif
    double p1 = pow(M(0,1), 2) + pow(M(0,2), 2) + pow(M(1,2), 2);
    if (p1 == 0)
    {
        Eigenvalues(0,0) = M(0,0);
        Eigenvalues(1,1) = M(1,1);
        Eigenvalues(2,2) = M(2,2);
    }
    else
    {
        double q = M.trace() / 3.0;
        double p2 = pow((M(0,0) - q), 2) + pow((M(1,1) - q), 2) + pow((M(2,2) - q), 2) + 2 * p1;
        double p = sqrt(p2 / 6.0);

        if (p != 0.0)
        {
            dMatrix3x3 I = dMatrix3x3::UnitTensor();
            dMatrix3x3 B = (M - I * q) * (1 / p);
            double r = 0.5 * B.determinant();

            double phi = 0;
            if (r <= -1.0)
            {
                phi = Pi / 3.0;
            }
            else if (r < 1.0)
            {
                phi = acos(r) / 3.0;
            }

            //the eigenvalues satisfy eig3 <= eig2 <= eig1
            Eigenvalues(0,0) = q + 2.0 * p * cos(phi);
            Eigenvalues(2,2) = q + 2.0 * p * cos(phi + (2.0*Pi / 3.0));
            Eigenvalues(1,1) = 3.0 * q - Eigenvalues(0,0) - Eigenvalues(2,2);
        }
        else
        {
            Eigenvalues(0,0) = 1;
            Eigenvalues(1,1) = 1;
            Eigenvalues(2,2) = 1;
        }
    }

    Eigenvectors.set_to_unity();
    if (calculateEigenVectors)
    {
        if (Eigenvalues(0,0) != Eigenvalues(1,1) and
            Eigenvalues(0,0) != Eigenvalues(2,2)/* and
            Eigenvalues(1,1) != Eigenvalues(2,2)*/)
        for (int n = 0; n < 3; n++)
        {
            dVector3 row0({ M(0,0) - Eigenvalues(n,n), M(0,1), M(0,2) });
            dVector3 row1({ M(1,0), M(1,1) - Eigenvalues(n,n), M(1,2) });
            dVector3 row2({ M(2,0), M(2,1), M(2,2) - Eigenvalues(n,n) });

            dVector3 r0xr1 = row0.cross(row1);
            dVector3 r0xr2 = row0.cross(row2);
            dVector3 r1xr2 = row1.cross(row2);

            double d0 = r0xr1 * r0xr1;
            double d1 = r0xr2 * r0xr2;
            double d2 = r1xr2 * r1xr2;
            double dMax = d0;
            int iMax = 0;

            if (d1 > dMax)
            {
                dMax = d1;
                iMax = 1;
            }
            if (d2 > dMax)
            {
                iMax = 2;
            }
            if (iMax == 0 and d0 > 0.0)
            {
                dVector3 temp = r0xr1 / sqrt(d0);
                Eigenvectors(n,0) = temp[0];
                Eigenvectors(n,1) = temp[1];
                Eigenvectors(n,2) = temp[2];
            }

            else if (iMax == 1 and d1 > 0.0)
            {
                dVector3 temp = r0xr2 / sqrt(d1);
                Eigenvectors(n,0) = temp[0];
                Eigenvectors(n,1) = temp[1];
                Eigenvectors(n,2) = temp[2];
            }
            else if (d2 > 0.0)
            {
                dVector3 temp = r1xr2 / sqrt(d2);
                Eigenvectors(n,0) = temp[0];
                Eigenvectors(n,1) = temp[1];
                Eigenvectors(n,2) = temp[2];
            }
        }
    }
}

void Tools::Eigensystem3(dMatrix3x3& M, dMatrix3x3& Eigenvectors, dMatrix3x3& Eigenvalues)
{
#ifdef DEBUG
    ///M MUST BE _SYMMETRIC_
    bool symmetric = true;
    for(int i = 0; i < 3; i++)
    for(int j = i; j < 3; j++)
    {
        if(fabs(M(i,j) - M(j,i)) > DBL_EPSILON)
            symmetric = false;
    }
    if(not symmetric)
    {
        Info::WriteExit("M:\n"
                        + M.print() + "is not symmetric!", "Tools", "Eigensystem2()");
        exit(13);
    }
#endif
    double NDvalue = 0.0;
    double theta = 0.0;
    Eigenvectors.set_to_unity();
    do
    {
        bool B12 = false;
        bool B13 = false;
        bool B23 = false;
        if(fabs(M(0,1)) > fabs(M(0,2)) and fabs(M(0,1)) > fabs(M(1,2)))
        {
            NDvalue = M(0,1);
            B12 = true;
        }
        else if (fabs(M(0,2)) > fabs(M(0,1)) and fabs(M(0,2)) > fabs(M(1,2)))
        {
            NDvalue = M(0,2);
            B13 = true;
        }
        else
        {
            NDvalue = M(1,2);
            B23 = true;
        }

        if(M(0,0) == M(1,1) or M(1,1) == M(2,2) or M(0,0) == M(2,2) )
        {
            theta = Pi/4.0;
        }
        else if (M(0,0) > M(1,1))
        {
            theta = 0.5 * atan(2 * NDvalue) / (M(0,0) - M(1,1));

        }
        else if (M(1,1) > M(2,2))
        {
            theta = 0.5 * atan(2 * NDvalue) / (M(1,1) - M(2,2));

        }
        else if (M(1,1) < M(2,2))
        {
            theta =  0.5 * atan(2 * NDvalue) / (M(2,2) - M(1,1));
        }

        if(B12)
        {
            dMatrix3x3 S12;
            S12(0,0) =  cos(theta); S12(0,1) = -sin(theta); S12(0,2) = 0.0;
            S12(1,0) =  sin(theta); S12(1,1) = cos(theta); S12(1,2) = 0.0;
            S12(2,0) =         0.0; S12(2,1) =        0.0; S12(2,2) = 1.0;
            M = S12.transposed() * M * S12;
            Eigenvectors = Eigenvectors * S12;
        }
        if(B13)
        {
            dMatrix3x3 S13;
            S13(0,0) =  cos(theta); S13(0,1) = 0.0; S13(0,2) = -sin(theta);
            S13(1,0) =         0.0; S13(1,1) = 1.0; S13(1,2) = 0.0;
            S13(2,0) = sin(theta); S13(2,1) = 0.0; S13(2,2) = cos(theta);
            M = S13.transposed() * M * S13;
            Eigenvectors = Eigenvectors * S13;
        }
        if(B23)
        {
            dMatrix3x3 S23;
            S23(0,0) = 1.0; S23(0,1) = 0.0;         S23(0,2) =        0.0;
            S23(1,0) = 0.0; S23(1,1) = cos(theta);  S23(1,2) = sin(theta);
            S23(2,0) = 0.0; S23(2,1) = -sin(theta); S23(2,2) = cos(theta);
            M = S23.transposed() * M * S23;
            Eigenvectors = Eigenvectors * S23;
        }
        for(int x = 0; x < 3; x++)
            for(int y = 0; y < 3; y++)
            {
                if(abs(M(x,y)) < (1e-6))
                {
                    M(x,y) = 0.0;
                }
            }
    } while ( abs(M(0,1)) > 0.0 or abs(M(0,2)) > 0.0 or abs(M(1,2)) > 0.0 );
    Eigenvalues = M;
}

dMatrix3x3 Tools::sqrtM3x3(dMatrix3x3 M)
{
    dMatrix3x3 Eigenvectors;
    dMatrix3x3 Eigenvalues;

    //Eigensystem(M,Eigenvectors,Eigenvalues);
    //Eigensystem2(M,Eigenvectors,Eigenvalues,true);
    Eigensystem3(M,Eigenvectors,Eigenvalues);

    for(int i = 0; i < 3; i++)
    {
        Eigenvalues(i,i) = sqrt(Eigenvalues(i,i));
    }
    return Eigenvectors * Eigenvalues * Eigenvectors.transposed();              //sqrtM = Eigenvectors_column . sqrt(Eigenvalues) . Eigenvectors_row;
}

dMatrix3x3 Tools::qurtM3x3(dMatrix3x3 M)
{
    dMatrix3x3 Eigenvectors;
    dMatrix3x3 Eigenvalues;

    //Eigensystem(M,Eigenvectors,Eigenvalues);
    //Eigensystem2(M,Eigenvectors,Eigenvalues,true);
    Eigensystem3(M,Eigenvectors,Eigenvalues);

    for(int i = 0; i < 3; i++)
    {
        Eigenvalues(i,i) = sqrt(sqrt(Eigenvalues(i,i)));
    }
    return Eigenvectors * Eigenvalues * Eigenvectors.transposed();              //sqrtM = Eigenvectors_column . sqrt(Eigenvalues) . Eigenvectors_row;
}

dMatrix3x3 Tools::logM3x3(dMatrix3x3 M)
{
    dMatrix3x3 Eigenvectors;
    dMatrix3x3 Eigenvalues;

    //Eigensystem(M,Eigenvectors,Eigenvalues);
    //Eigensystem2(M,Eigenvectors,Eigenvalues,true);
    Eigensystem3(M,Eigenvectors,Eigenvalues);

    for(int i = 0; i < 3; i++)
    {
        Eigenvalues(i,i) = log(Eigenvalues(i,i));
    }
    return Eigenvectors * Eigenvalues * Eigenvectors.transposed();              //logM = Eigenvectors_column . log(Eigenvalues) . Eigenvectors_row;
}

dMatrix3x3 Tools::expM3x3(dMatrix3x3 M)
{
    dMatrix3x3 Eigenvectors;
    dMatrix3x3 Eigenvalues;
    if(M.norm() < DBL_EPSILON)
    {
        return dMatrix3x3::UnitTensor();
    }
    //Eigensystem(M,Eigenvectors,Eigenvalues);
    //Eigensystem2(M,Eigenvectors,Eigenvalues,true);
    Eigensystem3(M,Eigenvectors,Eigenvalues);

    for(int i = 0; i < 3; i++)
    {
        Eigenvalues(i,i) = exp(Eigenvalues(i,i));
    }
    return Eigenvectors * Eigenvalues * Eigenvectors.transposed();              //expM = Eigenvectors_column . exp(Eigenvalues) . Eigenvectors_row;
}

void Tools::Decompose(dMatrix3x3& Deformation, dMatrix3x3& RotationMatrix, dMatrix3x3& StretchTensor)
{
    dMatrix3x3 M = Deformation.transposed() * Deformation;

    dMatrix3x3 Eigenvectors;
    dMatrix3x3 Eigenvalues;

    Eigensystem3(M,Eigenvectors,Eigenvalues);

    for(int i = 0; i < 3; i++)
    {
        Eigenvalues(i,i) = sqrt(Eigenvalues(i,i));
    }
    StretchTensor = Eigenvectors * Eigenvalues * Eigenvectors.transposed();
    RotationMatrix = Deformation*StretchTensor.inverted();
}

dMatrix3x3 Tools::AlignBaseAxes(dMatrix3x3 M)
{
    /*
     * This function uses Givens rotation algorithm to zero (1,0), (2,0) and
     * (2,1) elements of M.
    */
    double c = 0.0;
    double s = 0.0;
    double den = 0.0;
    dMatrix3x3 G;

    auto denom = [](double x, double y) {return sqrt(x*x + y*y);};

    // Zeroing element (1,0)
    if(M(1,0) > DBL_EPSILON or M(1,0) < -DBL_EPSILON)
    {
        den = denom (M(1,0),M(0,0));
        if(den > DBL_EPSILON)
        {
            c = M(0,0) / den;
            s = M(1,0) / den;
        }
        G(0,0) = c;  G(0,1) = s; G(0,2) = 0;
        G(1,0) = -s; G(1,1) = c; G(1,2) = 0;
        G(2,0) = 0;  G(2,1) = 0; G(2,2) = 1;
        M = G*M;
    }

    // Zeroing element (2,0)
    if(M(2,0) > DBL_EPSILON or M(2,0) < -DBL_EPSILON)
    {
        den = denom (M(2,0),M(0,0));
        if(den > DBL_EPSILON)
        {
            c = M(0,0) / den;
            s = M(2,0) / den;
        }
        G(0,0) = c;  G(0,1) = 0; G(0,2) = s;
        G(1,0) = 0;  G(1,1) = 1; G(1,2) = 0;
        G(2,0) = -s; G(2,1) = 0; G(2,2) = c;
        M = G*M;
    }

    // Zeroing element (2,1)
    if(M(2,1) > DBL_EPSILON or M(2,1) < -DBL_EPSILON)
    {
        den = denom (M(2,1),M(1,1));
        if(den > DBL_EPSILON)
        {
            c = M(1,1) / den;
            s = M(2,1) / den;
        }
        G(0,0) = 1; G(0,1) = 0;  G(0,2) = 0;
        G(1,0) = 0; G(1,1) = c;  G(1,2) = s;
        G(2,0) = 0; G(2,1) = -s; G(2,2) = c;
        M = G*M;
    }

    return M;
}


////////// extract Rotation from total deformation //////////
/// Implemented from https://dl.acm.org/doi/abs/10.1145/1028523.1028541
void Tools::jacobiRotate(dMatrix3x3 &A, dMatrix3x3 &R, int p, int q)
{
    if (A(p, q) == 0.0)
    {
        return;
    }
    double d = (A(p, p) - A(q, q)) / (2.0*A(p, q));
    double t = 1.0 / (fabs(d) + sqrt(d*d + 1.0));
    if (d < 0.0)
    {
        t = -t;
    }
    double c = 1.0 / sqrt(t*t + 1);
    double s = t*c;
    A(p, p) += t*A(p, q);
    A(q, q) -= t*A(p, q);
    A(p, q) = A(q, p) = 0.0;
    int k;
    for (k = 0; k < 3; k++)
    {
        if (k != p && k != q)
        {
            double Akp = c*A(k, p) + s*A(k, q);
            double Akq = -s*A(k, p) + c*A(k, q);
            A(k, p) = A(p, k) = Akp;
            A(k, q) = A(q, k) = Akq;
        }
    }
    for (k = 0; k < 3; k++)
    {
        double Rkp = c*R(k, p) + s*R(k, q);
        double Rkq = -s*R(k, p) + c*R(k, q);
        R(k, p) = Rkp;
        R(k, q) = Rkq;
    }
}
void Tools::eigenDecomposition(const dMatrix3x3 &A, dMatrix3x3 &EigenVectors, dVector3 &EigenValues)
{
    const int numJacobiIterations = 10;
    dMatrix3x3 D = A;
    EigenVectors.set_to_unity();
    int iter = 0;
    while (iter < numJacobiIterations)
    {
        int p, q;
        double a, max;
        max = fabs(D(0, 1));
        p = 0;
        q = 1;
        a = fabs(D(0, 2));
        if (a > max)
        {
            p = 0;
            q = 2;
            max = a;
        }
        a = fabs(D(1, 2));
        if (a > max)
        {
            p = 1;
            q = 2;
            max = a;
        }
        if (max < DBL_EPSILON)
        {
            break;
        }
        jacobiRotate(D, EigenVectors, p, q);
        iter++;
    }
    EigenValues[0] = D(0, 0);
    EigenValues[1] = D(1, 1);
    EigenValues[2] = D(2, 2);
}
void Tools::rotationMatrixIrving(const dMatrix3x3 &A, dMatrix3x3 &R)
{
    dMatrix3x3 AT_A, V;
    AT_A = A.transposed() * A;
    dVector3 S;
    eigenDecomposition(AT_A, V, S);
    const double detV = V.determinant();
    if (detV < 0.0)
    {
        double minLambda = DBL_MAX;
        unsigned char pos = 0;
        for (unsigned char l = 0; l < 3; l++)
        {
            if (S[l] < minLambda)
            {
                pos = l;
                minLambda = S[l];
            }
        }
        V(0, pos) = -V(0, pos);
        V(1, pos) = -V(1, pos);
        V(2, pos) = -V(2, pos);
    }
    if (S[0] < 0.0f)
        S[0] = 0.0f;
    if (S[1] < 0.0f)
        S[1] = 0.0f;
    if (S[2] < 0.0f)
        S[2] = 0.0f;

    dVector3 sigma;
    sigma[0] = sqrt(S[0]);
    sigma[1] = sqrt(S[1]);
    sigma[2] = sqrt(S[2]);
    unsigned char chk = 0;
    unsigned char pos = 0;
    dMatrix3x3 U;
    for (unsigned char l = 0; l < 3; l++)
    {
        if (fabs(sigma[l]) < 1.0e-4)
        {
            pos = l;
            chk++;
        }
    }
    if (chk > 0)
    {
        if (chk > 1)
        {
            U.set_to_unity();
        }
        else
        {
            U = A * V;
            for (unsigned char l = 0; l < 3; l++)
            {
                if (l != pos)
                {
                    for (unsigned char m = 0; m < 3; m++)
                    {
                        U(m, l) *= 1.0f / sigma[l];
                    }
                }
            }
            dVector3 v[2];
            unsigned char index = 0;
            for (unsigned char l = 0; l < 3; l++)
            {
                if (l != pos)
                {
                    v[index++] = dVector3({U(0, l), U(1, l), U(2, l)});
                }
            }
            dVector3 vec = v[0].cross(v[1]);
            vec.normalized();
            U(0, pos) = vec[0];
            U(1, pos) = vec[1];
            U(2, pos) = vec[2];
        }
    }
    else
    {
        dVector3 sigmaInv;
        sigmaInv = {1.0 / sigma[0], 1.0 / sigma[1], 1.0 / sigma[2]};
        U = A * V;
        for (unsigned char l = 0; l < 3; l++)
        {
            for (unsigned char m = 0; m < 3; m++)
            {
                U(m, l) *= sigmaInv[l];
            }
        }
    }
    const double detU = U.determinant();
    if (detU < 0.0)
    {
        double minLambda = DBL_MAX;
        unsigned char pos = 0;
        for (unsigned char l = 0; l < 3; l++)
        {
            if (sigma[l] < minLambda)
            {
                pos = l;
                minLambda = sigma[l];
            }
        }
        sigma[pos] = -sigma[pos];
        U(0, pos) = -U(0, pos);
        U(1, pos) = -U(1, pos);
        U(2, pos) = -U(2, pos);
    }
    R = U * V.transpose();
}
void Tools::GetAxisAngleFromRotationMatrix(const dMatrix3x3 RotMatrix, double& Angle, dVector3& Axis)
{
    ///  Implemented from https://www.euclideanspace.com/maths/geometry/rotations/index.htm
    double epsilon = 0.01; // margin to allow for rounding errors
    double epsilon2 = 0.1; // margin to distinguish between 0 and 180 degrees
    if ((fabs(RotMatrix(0,1)-RotMatrix(1,0)) < epsilon)
     && (fabs(RotMatrix(0,2)-RotMatrix(2,0)) < epsilon)
     && (fabs(RotMatrix(1,2)-RotMatrix(2,1)) < epsilon))
    {
        // singularity found
        // first check for identity matrix which must have +1 for all terms
        // in leading diagonal and zero in other terms
        if ((fabs(RotMatrix(0,1)+RotMatrix(1,0)) < epsilon2)
         && (fabs(RotMatrix(0,2)+RotMatrix(2,0)) < epsilon2)
         && (fabs(RotMatrix(1,2)+RotMatrix(2,1)) < epsilon2)
         && (fabs(RotMatrix(0,0)+RotMatrix(1,1)+RotMatrix(2,2)-3.0) < epsilon2))
        {
            // this singularity is identity matrix so angle = 0
            Angle = 0.0;
            Axis = {1,0,0};
            // zero angle, arbitrary axis
        }
        // otherwise this singularity is angle = 180
        Angle = Pi;
        double xx = (RotMatrix(0,0)+1)/2;
        double yy = (RotMatrix(1,1)+1)/2;
        double zz = (RotMatrix(2,2)+1)/2;
        double xy = (RotMatrix(0,1)+RotMatrix(1,0))/4;
        double xz = (RotMatrix(0,2)+RotMatrix(2,0))/4;
        double yz = (RotMatrix(1,2)+RotMatrix(2,1))/4;
        if ((xx > yy) && (xx > zz))
        {
            // RotMatrix(0,0) is the largest diagonal term
            if (xx< epsilon)
            {
                Axis = {0,0.7071,0.7071};
            }
            else
            {
                double x = sqrt(xx);
                Axis = {x,xy/x,xz/x};
            }
        }
        else if(yy > zz)
        {
            // m[1][1] is the largest diagonal term
            if (yy< epsilon)
            {
                Axis = {0.7071,0,0.7071};
            }
            else
            {
                double y = sqrt(yy);
                Axis = {xy/y,y,yz/y};
            }
        }
        else
        {
            // m[2][2] is the largest diagonal term so base result on this
            if (zz< epsilon)
            {
                Axis = {0.7071,0.7071,0};
            }
            else
            {
                double z = sqrt(zz);
                Axis = {xz/z,yz/z,z};
            }
        }
    }
    // as we have reached here there are no singularities so we can handle normally
    // First, normalise the rotation matrix
    double Norm = pow((RotMatrix(1,2) - RotMatrix(2,1)),2) + pow((RotMatrix(2,0) - RotMatrix(0,2)),2) + pow((RotMatrix(0,1) - RotMatrix(1,0)),2);
    Norm =  sqrt(Norm);

    if (fabs(Norm) < 0.001)
    {
        Norm = 1.0;
        // prevent divide by zero, should not happen if matrix is orthogonal and should be
        // caught by singularity test above, but I've left it in just in case
    }
    Angle = acos(0.5*(RotMatrix.trace() - 1));
    Axis[0] = ( RotMatrix(2,1) - RotMatrix(1,2) ) / Norm;
    Axis[1] = ( RotMatrix(0,2) - RotMatrix(2,0) ) / Norm;
    Axis[2] = ( RotMatrix(1,0) - RotMatrix(0,1) ) / Norm;
}
void Tools::getAxisAngle(const dMatrix3x3& TransformationMatrix, dVector3& Axis, double &Angle )
{
    dMatrix3x3 RotationMatrix;
    rotationMatrixIrving(TransformationMatrix,RotationMatrix);
    GetAxisAngleFromRotationMatrix(RotationMatrix,Angle,Axis);
}

dVector3 Tools::IPFColor(size_t sd, EulerAngles& tempEuler, Crystallography& CR)
{
    size_t index = 0;
    double Red = 0.0;
    double Green = 0.0;
    double Blue = 0.0;
    double MaxRGB = 0.0;

    dVector3 ref_dir;
    dVector3 hkl;
    dMatrix3x3 R;
    dMatrix3x3 tempOM;

    dVector3 RGBint;
    dVector3 RGB;

    double phi1 = tempEuler.Q[0];
    double PHI = tempEuler.Q[1];
    double phi2 = tempEuler.Q[2];

    // Assign reference sample direction
    switch (sd)
    {
        case 1: // 100
        {
            ref_dir = {1,0,0};
            break;
        }

        case 2: // 010
        {
            ref_dir = {0,1,0};
            break;
        }

        case 3: // 001
        {
            ref_dir = {0,0,1};
            break;
        }
    };

    // Start of main routine //
    // Assign black RGB values for bad data points (nsym = 0)
    if (CR.nsym == 0)
    {
        RGB.set_to_zero();
    }

    // Assign black RGB value for Euler angles outside of allowable range
    else if (phi1 > 2.0*M_PI || PHI > M_PI || phi2 > 2.0*M_PI)
    {
        RGB.set_to_zero();
    }

    //  Routine for valid set of Euler angles
    else
    {
        // Construct 3X3 orientation matrix from Euler Angles
        tempOM(0,0) = std::cos(phi1) * std::cos(phi2) - std::sin(phi1) * std::cos(PHI) * std::sin(phi2);
        tempOM(0,1) = std::sin(phi1) * std::cos(phi2) + std::cos(phi1) * std::cos(PHI) * std::sin(phi2);
        tempOM(0,2) = std::sin(phi2) * std::sin(PHI);
        tempOM(1,0) = -std::cos(phi1) * std::sin(phi2) - std::sin(phi1) * std::cos(PHI) * std::cos(phi2);
        tempOM(1,1) = -std::sin(phi1) * std::sin(phi2) + std::cos(phi1) * std::cos(PHI) * std::cos(phi2);
        tempOM(1,2) = std::cos(phi2) * std::sin(PHI);
        tempOM(2,0) = std::sin(phi1) * std::sin(PHI);
        tempOM(2,1) = -std::cos(phi1) * std::sin(PHI);
        tempOM(2,2) = std::cos(PHI);

        //Sorting Euler angles into standard stereographic triangle (SST)
        index = 0;
        while(index < CR.nsym)
        {
            // Form orientation matrix
            for (size_t i = 0; i < 3; ++i)
            {
                for (size_t j = 0; j < 3; ++j)
                {
                    R(i,j) = 0.0;
                    for (size_t k = 0; k < 3; ++k)
                    {
                        R(i,j) += CR.CrystalSymmetries[index](i,k) * tempOM(k,j);
                    }
                }
            }

            // Multiple orientation matrix by reference sample direction
            for (size_t i = 0; i < 3; ++i)
            {
                hkl[i] = 0;
                for (size_t j = 0; j < 3; ++j)
                {
                    hkl[i] += R(i,j) * ref_dir[j];
                }
            }

            // Convert to spherical coordinates (ignore "r" variable since r=1)
            CR.Theta = abs(atan2(hkl[1], hkl[0]));
            CR.Phi = acos(abs(hkl[2]));

            // Continue if Theta and Phi values are within the SST
            if (CR.Theta >= CR.Theta_min && CR.Theta < CR.Theta_max && CR.Phi >= CR.Phi_min && CR.Phi < CR.Phi_max)
            {
                break;
            }

            // Increment to next symmetry operator if not in SST
            else
            {
                index++;
            }
        }

        //  Adjust maximum Phi value to ensure it falls within the SST (cubic materials only)
        if(CR.nsym == 24)
        {
            CR.Phi_max2 = acos(sqrt(1.0 / (2.0 + (pow(tan(CR.Theta),2.0)))));
        }
        else
        {
            CR.Phi_max2 = M_PI / 2;
        }

        // Calculate the RGB color values and make adjustments to maximize colorspace
        Red = abs(1.0 - (CR.Phi / CR.Phi_max2));
        Blue = abs((CR.Theta - CR.Theta_min) / (CR.Theta_max - CR.Theta_min));
        Green = 1.0 - Blue;

        Blue *= (CR.Phi / CR.Phi_max2);
        Green *= (CR.Phi / CR.Phi_max2);

        // Check for negative RGB values before taking square root
        if(Red < 0 || Green < 0 || Blue < 0)
        {
            string msg = "RGB component values must be positive!";
            Info::WriteWarning(msg,"Tools()","Euler2rgb");
        }

        RGB[0] = sqrt(Red);
        RGB[1] = sqrt(Green);
        RGB[2] = sqrt(Blue);

        // Find maximum value of red, green, or blue
        MaxRGB = max({RGB[0], RGB[1], RGB[2]});

        // Normalize position of SST center point
        RGB /= MaxRGB;
    }
    return RGB;
}

/*
 * A Robust Method to Extract the Rotational Part of Deformations
 * https://animation.rwth-aachen.de/media/papers/2016-MIG-StableRotation.pdf
 */

/* The quat can be :
 *
 * 1- The solution of the previous step of an iterative solve.
 * 2- If such a solution is not available, start with q = Quat(Mat)/|Quat(Mat)|.
 * 3- q=(1,0,0,0).
 *
 *
 * No of iterations: tested, maxIter: 20 for best results!
 *
 *
*/
void Tools::ExtractRotation(dMatrix3x3 &Mat,  Quaternion& Quat, const size_t maxIter)
{
    for (size_t iter = 0; iter < maxIter; iter++)
    {
        dMatrix3x3 Rot;
        Quaternion tempQuat;
        Rot = Quat.getRotationMatrix();
        std::vector<dVector3> Rot_col = Col(Rot);
        std::vector<dVector3> Mat_col = Col(Mat);

        double b = (1.0 / fabs(Rot_col[0]*Mat_col[0] + Rot_col[1]*Mat_col[1]
                        + Rot_col[2]*Mat_col[2]) + 1.0e-9);

        dVector3 omega = (Rot_col[0].cross(Mat_col[0]) + Rot_col[1].cross(Mat_col[1])
                + Rot_col[2].cross(Mat_col[2]))*b;

        double Angle = omega.length();

        if (Angle < 1.0e-9) break;

        dVector3 Axis = omega * (1.0 / Angle);
        tempQuat.set(Axis, Angle);
        tempQuat.normalize();
        Quat = tempQuat*Quat;
        Quat.normalize();
    }
}
}// namespace openphase
