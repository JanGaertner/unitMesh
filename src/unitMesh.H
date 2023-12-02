/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
 
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
 
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
 
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
 
Class
    Foam::unitMesh
 
Description
    Create a cubic fvMesh object with length 1m and nCells in each dimension
 
SourceFiles
    unitMesh.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2023
 
\*---------------------------------------------------------------------------*/


#ifndef unitMesh_H
#define unitMesh_H


#include "fvMesh.H"

namespace Foam
{

class unitMesh
{
    private:

        //- Pointer to store the unit mesh
        autoPtr<fvMesh> mesh_;

        //- Default system dictionary to read for runTime
        dictionary systemDict_;

        //- Store a time object as a new objectRegistry
        autoPtr<Time> runTime_;

        //- Create the Time object and object registry
        void createObjectRegistry(const fileName& rootPath);

        //- Add boundaries
        void addBoundaries
        (
            const label nCells,
            const label startIndex
        );

    public:

        //- Construct a unit mesh with nCells in each dimension 
        unitMesh
        (
            const label& nCells,
            const fileName& rootPath="./"
        );

    // Access

        //- Access to the unit mesh
        const fvMesh& mesh();

        //- Access to the objectRegistry over Time
        const Time& runTime();

};
}

#endif
