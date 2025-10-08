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
 *   File created :   2014
 *   Main contributors :   Efim Borukhovich; Oleg Shchyglo
 *
 */

/*
 * OPObject is a parent virtual class for most of the stand-alone classes in
 * the OpenPhase. It defines the names of common methods and parameters which
 * unifies their syntax across the library. The common base also allows creating
 * lists of heterogeneous objects and allows batch execution of common methods,
 * e.g. Remesh() and Advect() on a list of different objects.
 */

#ifndef OPOBJECT_H
#define OPOBJECT_H

#include "Base/Includes.h"
#include "Info.h"
namespace openphase
{

class Settings;
class BoundaryConditions;
class PhaseField;
class Velocities;
class AdvectionHR;

class OP_EXPORTS OPObject
{
 public:
    virtual ~OPObject(void);

    std::string thisclassname;                                                  ///< Object's implementation class name
    std::string thisobjectname;                                                 ///< Object's name
    bool initialized = false;                                                   ///< Object's initialization status flag
    bool remeshable = false;                                                    ///< True if the object has non-empty Remesh() method
    bool advectable = false;                                                    ///< True if the object has non-empty Advect() method

    virtual void Initialize(Settings& locSettings){(void)locSettings;};         ///< Initializes internal variables and storages
    virtual void ReadInput(const std::string InputFileName){(void)InputFileName;};///< Reads input data from the user specified input file
    virtual void ReadInput(std::stringstream& inp){(void)inp;};///< Reads input data from the user specified input file

    virtual void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)///< Remeshes the storages
    {
        (void) newNx; //unused
        (void) newNy; //unused
        (void) newNz; //unused
        (void) BC;    //unused
    };

    virtual void Advect(AdvectionHR& Adv, const Velocities& Vel,
                        PhaseField& Phi, const BoundaryConditions& BC,
                        const double dt, const double tStep)
    {
        (void) Adv;   //unused
        (void) Vel;   //unused
        (void) Phi;   //unused
        (void) BC;    //unused
        (void) dt;    //unused
        (void) tStep; //unused
    };

    static OPObject* findOPObject(std::vector<OPObject*> ObjectList,
                           std::string ObjectName, std::string thisclassname,
                           std::string thisfunctionname, bool necessary, bool verbose = true)
    /* Returns a reference to the first OPObject in the list whose class name
     * starts with ObjectName. I.e.: if ObjectName is "Elasticity", the first
     * object in the list which has either the name ElasticitySteinbach or
     * ElasticityKhachaturyan will be returned.
     */
    {
        for(unsigned int i = 0; i < ObjectList.size(); i++)
        {
            if(ObjectList[i]->thisclassname.size() >= ObjectName.size())
            {
                std::string locObjectName = ObjectList[i]->thisclassname;
                locObjectName.resize(ObjectName.size());

                if(locObjectName == ObjectName)
                    return ObjectList[i];
            }
        }
        if (verbose)
        {
            std::cout << "No " << ObjectName << " object found! " << thisclassname
                      << "."<<thisfunctionname <<"()"<< std::endl;
        }

        if(necessary)
        {
            exit(13);
        }

        return nullptr;
    }

    static std::vector<OPObject*> findOPObjects(std::vector<OPObject*> ObjectList,
                           std::string ObjectName, std::string thisclassname,
                           std::string thisfunctionname, bool necessary, bool verbose = true)
    /* Returns a list of OPObject references found in the passed ObjectList
     * whose class name start with ObjectName. I.e.: if ObjectName is
     * "Elasticity", all objects in the ObjectList which have a name starting
     * with "Elasticity", e.g. ElasticitySteinbach or ElasticityKhachaturyan
     * will be returned in the result list.
     */
    {
        std::vector<OPObject*> result;

        for(unsigned int i = 0; i < ObjectList.size(); i++)
        {
            if(ObjectList[i]->thisclassname.size() >= ObjectName.size())
            {
                std::string locObjectName = ObjectList[i]->thisclassname;
                locObjectName.resize(ObjectName.size());

                if(locObjectName == ObjectName)
                    result.push_back(ObjectList[i]);
            }
        }

        if(result.size() == 0)
        {
            if (verbose)
            {
                std::cout << "No " << ObjectName << " objects found! " << thisclassname
                          << thisfunctionname << std::endl;
            }
            if(necessary)
            {
                exit(13);
            }
        }
        return result;
    }

 protected:
 private:
};

}// namespace openphase
#endif
