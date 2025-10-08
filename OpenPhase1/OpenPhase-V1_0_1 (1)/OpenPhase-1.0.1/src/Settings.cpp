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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */

#include "Base/UserInterface.h"
#include "BoundaryConditions.h"
#include "Info.h"
#include "RunTimeControl.h"
#include "Settings.h"

namespace openphase
{

using namespace std;

void Settings::Initialize(void)
{
    thisclassname = "Settings";
    thisobjectname = thisclassname;

    Nx = 1;
    Ny = 1;
    Nz = 1;

    dNx = 1;
    dNy = 1;
    dNz = 1;

    newNx = Nx;
    newNy = Ny;
    newNz = Nz;

    OffsetX = 0;
    TotalNx = Nx;

    OffsetY = 0;
    TotalNy = Ny;

    OffsetZ = 0;
    TotalNz = Nz;

    maxNx = Nx;
    maxNy = Ny;
    maxNz = Nz;

    dx = 1;
    iWidth = 1;
    Eta = iWidth*dx;

    Nphases = 0;
    Ncomp = 0;

    Resolution = Resolutions::Single;
    ConsiderNucleusVolume = false;

    DiffusionStencil = LaplacianStencils::Isotropic;
    PhaseFieldLaplacianStencil = LaplacianStencils::Isotropic;
    PhaseFieldGradientStencil = GradientStencils::Simple;

    if(VTKDir.size() == 0) VTKDir = DefaultVTKDir;                              //Only overwrite default directory if VTKDir is empty. This way, VTKDir can be defined earlier.
    if(RawDataDir.size() == 0) RawDataDir = DefaultRawDataDir;                  //Only overwrite default directory if RawDataDir is empty. This way, RawDataDir can be defined earlier.
    if(TextDir.size() == 0) TextDir = DefaultTextDir;                           //Only overwrite default directory if TextDir is empty. This way, TextDir can be defined earlier.

    initialized = true;
    Info::WriteStartScreen();
    Info::WriteStandard(thisclassname, "Initialized");
}

void Settings::ReadInput(const string InputFileName)
{
    Info::WriteLineInsert("Settings input");
    Info::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };

    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);

    inp.close();
}

void Settings::ReadInput(std::stringstream& inp)
{
    Info::WriteLine();
    Info::WriteLineInsert("System Dimensions");

    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    Nx = UserInterface::ReadParameterI(inp, moduleLocation, string("Nx"));
    Ny = UserInterface::ReadParameterI(inp, moduleLocation, string("Ny"));
    Nz = UserInterface::ReadParameterI(inp, moduleLocation, string("Nz"));

    Bcells = UserInterface::ReadParameterI(inp, moduleLocation, string("Bcells"), false, 1);

    if(Nx == 0) // reduced X dimension
    {
        Nx  = 1;
        dNx = 0;
    }
    if(Ny == 0) // reduced Y dimension
    {
        Ny  = 1;
        dNy = 0;
    }
    if(Nz == 0) // reduced Z dimension
    {
        Nz  = 1;
        dNz = 0;
    }

    if (dNx + dNy + dNz == 0)
    {
        Info::WriteExit("All dimensions are suppressed: Nx=Ny=Nz=0!\n OpenPhase is not intended for simulations of zero-dimensional systems", "Settings", "ReadInput()");
        exit(1);
    }

    TotalNx = Nx;
    TotalNy = Ny;
    TotalNz = Nz;

    newNx = Nx;
    newNy = Ny;
    newNz = Nz;

    maxNx = Nx;
    maxNy = Ny;
    maxNz = Nz;

    DimensionsHistory.push_back(std::array<int, 4>({0, TotalNx, TotalNy, TotalNz}));

    dx     = UserInterface::ReadParameterD(inp, moduleLocation, string("dx"));
    iWidth = UserInterface::ReadParameterD(inp, moduleLocation, string("IWidth"));

    ConsiderNucleusVolume = UserInterface::ReadParameterB(inp, moduleLocation, string("ConsiderNucleusVolume"), false, false);

    string locResolution = UserInterface::ReadParameterK(inp, moduleLocation, string("Resolution"), false, "SINGLE");
    bool resolution_set = false;
    if(locResolution == "DOUBLE")
    {
        Resolution = Resolutions::Double;
        dx_2 = 0.5*dx;
        Eta = iWidth*dx_2;
        resolution_set = true;
    }
    if(locResolution == "SINGLE")
    {
        Resolution = Resolutions::Single;
        dx_2 = dx;
        Eta = iWidth*dx;
        resolution_set = true;
    }
    if(!resolution_set)
    {
        Info::WriteExit("Wrong resolution selected -> " + locResolution, "Settings", "ReadInput()");
        exit(1);
    }

    string tmp1 = UserInterface::ReadParameterK(inp, moduleLocation, string("DiffusionStencil"), false, string("ISOTROPIC"));
    if(tmp1 == "SIMPLE")
    {
        DiffusionStencil = LaplacianStencils::Simple;
    }
    if(tmp1 == "ISOTROPIC")
    {
        DiffusionStencil = LaplacianStencils::Isotropic;
    }
    if(tmp1 == "LB")
    {
        DiffusionStencil = LaplacianStencils::LB;
    }

    string tmp2 = UserInterface::ReadParameterK(inp, moduleLocation, string("PhaseFieldLaplacianStencil"), false, string("ISOTROPIC"));
    if(tmp2 == "SIMPLE")
    {
        PhaseFieldLaplacianStencil = LaplacianStencils::Simple;
    }
    if(tmp2 == "ISOTROPIC")
    {
        PhaseFieldLaplacianStencil = LaplacianStencils::Isotropic;
    }
    if(tmp2 == "LB")
    {
        PhaseFieldLaplacianStencil = LaplacianStencils::LB;
    }
    string tmp3 = UserInterface::ReadParameterK(inp, moduleLocation, string("PhaseFieldGradientStencil"), false, string("SIMPLE"));
    if(tmp3 == "SIMPLE")
    {
        PhaseFieldGradientStencil = GradientStencils::Simple;
    }
    if(tmp3 == "ISOTROPIC")
    {
        PhaseFieldGradientStencil = GradientStencils::Isotropic;
    }
    if(tmp3 == "LB")
    {
        PhaseFieldGradientStencil = GradientStencils::LB;
    }

    move_frame_phase     = UserInterface::ReadParameterI(inp, moduleLocation, string("MoveFramePhase"), false, -1);
    move_frame_pos       = UserInterface::ReadParameterI(inp, moduleLocation, string("MoveFramePosition"), false, -1);
    move_frame_direction = UserInterface::ReadParameterS(inp, moduleLocation, string("MoveFrameDirection"), false, "No");

    Info::WriteLineInsert("Active phases");

    bool endofnames = false;
    size_t n = 0;

    while(!(endofnames))
    {
        stringstream converter;
        converter << string("Phase_") << n;
        if(UserInterface::FindParameter(inp, moduleLocation, converter.str()) != -1)
        {
            string tmp = UserInterface::ReadParameterK(inp, moduleLocation, converter.str());
            PhaseNames.push_back(tmp);
            n++;
        }
        else
        {
            endofnames = true;
        }
    }
    Nphases = PhaseNames.size();

    // Reading states of matter for all phases.
    PhaseAggregateStates.resize(Nphases);
    for(size_t m = 0; m < Nphases; m++)
    {
        stringstream converter;
        converter << string("State_") << m;
        string state_of_matter = UserInterface::ReadParameterK(inp, moduleLocation, converter.str(), false, "SOLID");
        bool state_set = false;
        if(state_of_matter == "SOLID")
        {
            PhaseAggregateStates[m] = AggregateStates::Solid;
            state_set = true;
        }
        if(state_of_matter == "LIQUID" )
        {
            PhaseAggregateStates[m] = AggregateStates::Liquid;
            state_set = true;
        }
        if(state_of_matter == "GAS")
        {
            PhaseAggregateStates[m] = AggregateStates::Gas;
            state_set = true;
        }
        if(!state_set)
        {
            Info::WriteExit("Wrong state of matter is selected for phase " + to_string(m) + " -> " + state_of_matter, "Settings", "ReadInput()");
            exit(1);
        }
    }

    // Reading number of crystallographic variants of each phase
    Nvariants.resize(Nphases);
    for(size_t pIndex = 0; pIndex < Nphases; pIndex++)
    {
        stringstream converter;
        converter << "Nvariants_" << pIndex;
        Nvariants[pIndex] = UserInterface::ReadParameterI(inp, moduleLocation, converter.str(),false,0);
    }

    endofnames = false;
    n = 0;

    while(!(endofnames))
    {
        stringstream converter;
        converter << string("Comp_") << n;
        if(UserInterface::FindParameter(inp, moduleLocation, converter.str()) != -1)
        {
            string tmp = UserInterface::ReadParameterK(inp, moduleLocation, converter.str());
            if(tmp == "VA")
            {
                cerr << "\"VA\" as a component is not allowed in "
                     << "ChemicalProperties.opi" << endl;
                exit(0);
            }
            ElementNames.push_back(tmp);
            n++;
        }
        else
        {
            endofnames = true;
        }
    }
    Ncomp = ElementNames.size();

    VTKDir     = UserInterface::ReadParameterF(inp, moduleLocation,"VTKDir",false,VTKDir);
    RawDataDir = UserInterface::ReadParameterF(inp, moduleLocation,"RAWDir",false,RawDataDir);
    TextDir    = UserInterface::ReadParameterF(inp, moduleLocation,"DATADir",false,TextDir);

    string from = "\\";
    string to = "/";

#ifdef _WIN32
    from = "/";
    to = "\\";
#endif

    size_t start_pos = 0;
    while((start_pos = VTKDir.find(from, start_pos)) != std::string::npos)
    {
        VTKDir.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }

    start_pos = 0;
    while((start_pos = RawDataDir.find(from, start_pos)) != std::string::npos)
    {
        RawDataDir.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }

    start_pos = 0;
    while((start_pos = TextDir.find(from, start_pos)) != std::string::npos)
    {
        TextDir.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }

    if (!VTKDir.empty())
    {
        stringstream ss; ss << *VTKDir.rbegin();
        string s; ss >> s;
        if(s != dirSeparator) VTKDir += dirSeparator;
    }

    if (!RawDataDir.empty())
    {
        stringstream ss; ss << *RawDataDir.rbegin();
        string s; ss >> s;
        if(s != dirSeparator) RawDataDir += dirSeparator;
    }

    if (!TextDir.empty())
    {
        stringstream ss; ss << *TextDir.rbegin();
        string s; ss >> s;
        if(s != dirSeparator) TextDir += dirSeparator;
    }

    Info::WriteLine();
    Info::WriteBlankLine();

#ifdef MPI_PARALLEL
    MPI_3D_DECOMPOSITION = UserInterface::ReadParameterB(inp, moduleLocation, std::string("MPI3D"), false, false);

    if(MPI_3D_DECOMPOSITION)
    {
        MPI_CART_SIZE[0] = UserInterface::ReadParameterI(inp, moduleLocation, string("Ncx"), false, 1);
        MPI_CART_SIZE[1] = UserInterface::ReadParameterI(inp, moduleLocation, string("Ncy"), false, 1);
        MPI_CART_SIZE[2] = UserInterface::ReadParameterI(inp, moduleLocation, string("Ncz"), false, 1);

        if(MPI_CART_SIZE[0] * MPI_CART_SIZE[1] * MPI_CART_SIZE[2] != MPI_SIZE)
        {
            string message = "The requested number of MPI parallel blocks can not be decomposed by requested number of 3D blocks.";
            Info::WriteExit(message, thisclassname, "ReadInput()");
            exit(1);
        }
        Setup_MPI3D();
    }
    else
    {
        Setup_MPI();
    }
#endif
}
int Settings::ActiveDimensions(void) const
{
    return dNx + dNy + dNz;
}
void Settings::Resize(int size_x, int size_y, int size_z, int tStep, const BoundaryConditions& BC)
{
    if(dNx) newNx = size_x;
    if(dNy) newNy = size_y;
    if(dNz) newNz = size_z;

    if(RemeshingAllowed and (Nx != newNx or Ny != newNy or Nz != newNz))
    {
        RemeshAll(BC);

        if(dNx) Nx = newNx;
        if(dNy) Ny = newNy;
        if(dNz) Nz = newNz;

        /// TODO: add MPI update of system dimensions
        TotalNx = Nx;
        TotalNy = Ny;
        TotalNz = Nz;

        maxNx = max(Nx, maxNx);
        maxNy = max(Ny, maxNy);
        maxNz = max(Nz, maxNz);

        DimensionsHistory.push_back(std::array<int, 4>({tStep, TotalNx, TotalNy, TotalNz}));
    }
}

void Settings::ReadDimensions(int tStep)
{
    string FileName = UserInterface::MakeFileName(RawDataDir,"SystemDimensions_", tStep, ".dat");

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be opened", thisclassname, "Read");
        exit(1);
    };
    string tmp;

    getline(inp, tmp); // discards the headline

    inp >> dNx >> dNy >> dNz; //reads active dimensions

    getline(inp, tmp); // discards another headline

    while(!inp.eof())
    {
        int loctStep;

        inp >> loctStep >> Nx >> Ny >> Nz;

        /// TODO: add MPI update of system dimensions
        TotalNx = Nx;
        TotalNy = Ny;
        TotalNz = Nz;

        DimensionsHistory.push_back(std::array<int, 4>({loctStep, TotalNx, TotalNy, TotalNz}));

        maxNx = max(Nx, maxNx);
        maxNy = max(Ny, maxNy);
        maxNz = max(Nz, maxNz);
    }
}

void Settings::WriteDimensions(int tStep)
{
    string FileName = UserInterface::MakeFileName(RawDataDir,"SystemDimensions_", tStep, ".dat");

    ofstream out(FileName.c_str(), ios::out);

    if (!out)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be created", thisclassname, "Write");
        exit(1);
    };
    out << "Active dimensions:" << endl;
    out << dNx << "\t" << dNy << "\t" << dNz << endl;

    out << "tStep" << "\t" << "Nx" << "\t" << "Ny" << "\t" << "Nz" << endl;
    for(size_t n = 0; n < DimensionsHistory.size(); n++)
    {
        out << DimensionsHistory[n][0] << "\t" << DimensionsHistory[n][1] << "\t" << DimensionsHistory[n][2] << "\t" << DimensionsHistory[n][3] << endl;
    }
}

void Settings::AddForRemesh(OPObject& Obj)
{
    ObjectsToRemesh.push_back(&Obj);
    NamesOfObjectsToRemesh.push_back(Obj.thisclassname);
}

void Settings::RemeshAll(const BoundaryConditions& BC)
{
    for(size_t n = 0; n < ObjectsToRemesh.size(); n++)
    {
        ObjectsToRemesh[n]->Remesh(newNx, newNy, newNz, BC);
    }
}

void Settings::Remesh(std::string ObjNameBase, int size_x, int size_y, int size_z, const BoundaryConditions& BC)
{
    if(dNx) newNx = size_x;
    if(dNy) newNy = size_y;
    if(dNz) newNz = size_z;

    for(size_t n = 0; n < ObjectsToRemesh.size(); n++)
    if(NamesOfObjectsToRemesh[n] == ObjNameBase)
    {
        ObjectsToRemesh[n]->Remesh(newNx, newNy, newNz, BC);
    }
}

void Settings::AddForAdvection(OPObject& Obj)
{
    ObjectsToAdvect.push_back(&Obj);
    NamesOfObjectsToAdvect.push_back(Obj.thisclassname);
}

void Settings::AdvectAll(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi, const BoundaryConditions& BC, RunTimeControl& RTC)
{
    for(size_t n = 0; n < ObjectsToAdvect.size(); n++)
    {
        ObjectsToAdvect[n]->Advect(Adv, Vel, Phi, BC, RTC.dt, RTC.tStep);
    }
}

void Settings::Advect(std::string ObjNameBase, AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi, const BoundaryConditions& BC, RunTimeControl& RTC)
{
    for(size_t n = 0; n < ObjectsToAdvect.size(); n++)
    if(NamesOfObjectsToAdvect[n] == ObjNameBase)
    {
        ObjectsToAdvect[n]->Advect(Adv, Vel, Phi, BC, RTC.dt, RTC.tStep);
    }
}

Settings& Settings::operator= (const Settings& rhs)
{
    if (this != &rhs) // protect against self-assignment
    {
        thisclassname = rhs.thisclassname;
        thisobjectname = rhs.thisobjectname;

        TotalNx = rhs.TotalNx;
        OffsetX = rhs.OffsetX;
        TotalNy = rhs.TotalNy;
        OffsetY = rhs.OffsetY;
        TotalNz = rhs.TotalNz;
        OffsetZ = rhs.OffsetZ;

        Nx = rhs.Nx;
        Ny = rhs.Ny;
        Nz = rhs.Nz;

        dNx = rhs.dNx;
        dNy = rhs.dNy;
        dNz = rhs.dNz;

        newNx = rhs.newNx;
        newNy = rhs.newNy;
        newNz = rhs.newNz;

        maxNx = rhs.maxNx;
        maxNy = rhs.maxNy;
        maxNz = rhs.maxNz;

        dx     = rhs.dx;
        iWidth = rhs.iWidth;
        Eta    = rhs.Eta;
        Bcells = rhs.Bcells;

        ConsiderNucleusVolume = rhs.ConsiderNucleusVolume;

        Nphases = rhs.Nphases;
        Ncomp   = rhs.Ncomp;

        Resolution = rhs.Resolution;
        DiffusionStencil = rhs.DiffusionStencil;
        PhaseFieldLaplacianStencil = rhs.PhaseFieldLaplacianStencil;
        PhaseFieldGradientStencil = rhs.PhaseFieldGradientStencil;

        RemeshingAllowed = rhs.RemeshingAllowed;
        initialized = rhs.initialized;

        DimensionsHistory = rhs.DimensionsHistory;
        ObjectsToRemesh = rhs.ObjectsToRemesh;
        NamesOfObjectsToRemesh = rhs.NamesOfObjectsToRemesh;

        PhaseNames = rhs.PhaseNames;
        ElementNames = rhs.ElementNames;

        VTKDir = rhs.VTKDir;
        RawDataDir = rhs.RawDataDir;
        TextDir = rhs.TextDir;
    }
    return *this;
}

#ifdef MPI_PARALLEL
void Settings::Setup_MPI()
{
    /// 1D MPI domain decomposition along outer (X) dimension

    int blocksize = Nx / MPI_SIZE;
    int blocksize_rest = Nx % MPI_SIZE;
    int position = blocksize * MPI_RANK;
    if (blocksize_rest != 0)
    {
        if (MPI_RANK < MPI_SIZE - 1)
        {
            blocksize += 1;
            position += MPI_RANK;
        }
        else if (MPI_RANK == MPI_SIZE - 1)
        {
            blocksize -= MPI_SIZE - 1 - blocksize_rest;
            position  += MPI_RANK;
        }
    }

    Nx = blocksize;
    OffsetX = position;
    if (Nx < iWidth)
    {
        Info::WriteExit("Dimension X is too small for the requested number of MPI processes", thisclassname, "Setup_MPI()");
        exit(1);
    }
    Info::WriteStandard(thisclassname + "(RANK " + std::to_string(MPI_RANK) + ")", "MPI environment is initialized");
}

void Settings::Setup_MPI3D()
{
    MPI_Comm cart_comm;
    int reorder = 0;
    int periodic[3] = {0,0,0};

    MPI_Cart_create(MPI_COMM_WORLD, 3, MPI_CART_SIZE, periodic, reorder, &cart_comm);

    MPI_Cart_coords(cart_comm, MPI_RANK, 3, MPI_CART_RANK);

    TotalNx = Nx;
    TotalNy = Ny;
    TotalNz = Nz;

    int blocksizeX = Nx / MPI_CART_SIZE[0];
    int blocksizeY = Ny / MPI_CART_SIZE[1];
    int blocksizeZ = Nz / MPI_CART_SIZE[2];

    int blocksize_restX = Nx % MPI_CART_SIZE[0];
    int blocksize_restY = Ny % MPI_CART_SIZE[1];
    int blocksize_restZ = Nz % MPI_CART_SIZE[2];

    int positionX = MPI_CART_RANK[0] * blocksizeX;
    int positionY = MPI_CART_RANK[1] * blocksizeY;
    int positionZ = MPI_CART_RANK[2] * blocksizeZ;

    if (blocksize_restX != 0 and MPI_CART_RANK[0] < MPI_CART_SIZE[0] - 1)
    {
        blocksizeX += 1;
        positionX += MPI_CART_RANK[0];
    }
    if(blocksize_restX != 0 and MPI_CART_RANK[0] == MPI_CART_SIZE[0] - 1)
    {
        blocksizeX -= MPI_CART_SIZE[0] - 1 - blocksize_restX;
        positionX  += MPI_CART_RANK[0];
    }

    if (blocksize_restY != 0 and MPI_CART_RANK[1] < MPI_CART_SIZE[1] - 1)
    {
        blocksizeY += 1;
        positionY += MPI_CART_RANK[1];
    }
    if(blocksize_restY != 0 and MPI_CART_RANK[1] == MPI_CART_SIZE[1] - 1)
    {
        blocksizeY -= MPI_CART_SIZE[1] - 1 - blocksize_restY;
        positionY  += MPI_CART_RANK[1];
    }

    if (blocksize_restZ != 0 and MPI_CART_RANK[2] < MPI_CART_SIZE[2] - 1)
    {
        blocksizeZ += 1;
        positionZ += MPI_CART_RANK[2];
    }
    if(blocksize_restZ != 0 and MPI_CART_RANK[2] == MPI_CART_SIZE[2] - 1)
    {
        blocksizeZ -= MPI_CART_SIZE[2] - 1 - blocksize_restZ;
        positionZ  += MPI_CART_RANK[2];
    }

    Nx = blocksizeX;
    OffsetX = positionX;

    Ny = blocksizeY;
    OffsetY = positionY;

    Nz = blocksizeZ;
    OffsetZ = positionZ;

    if (dNx > 0)
    if (Nx < iWidth)
    {
        Info::WriteExit("Dimension X is too small for the number of MPI processes", thisclassname, "Setup_MPI3D()");
        exit(1);
    }
    if (dNy > 0)
    if (Ny < iWidth)
    {
        Info::WriteExit("Dimension Y is too small for the number of MPI processes", thisclassname, "Setup_MPI3D()");
        exit(1);
    }
    if (dNz > 0)
    if (Nz < iWidth)
    {
        Info::WriteExit("Dimension Z is too small for the number of MPI processes", thisclassname, "Setup_MPI3D()");
        exit(1);
    }
    Info::WriteStandard(thisclassname, "MPI3D environment is initialized");
}
#endif

} //namespace openphase
