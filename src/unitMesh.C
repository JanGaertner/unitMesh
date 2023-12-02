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

\*---------------------------------------------------------------------------*/

#include "unitMesh.H"

Foam::unitMesh::unitMesh
(
    const label& nCells,
    const fileName& rootPath    
)
{

    // Create object registry
    createObjectRegistry(rootPath);

    // Create cells and faces
    pointField points(std::pow(nCells+1,3));
    faceList   faces(3*(std::pow(nCells,2)*(nCells+1)));
    labelList   owner(faces.size());
    labelList   neighbor(faces.size()-6*std::pow(nCells,2));

    // To transform from the 1D representation to the 3D access for points
    const label yOffset = nCells+1;
    const label zOffset = std::pow(nCells+1,2);


    // Construction of the point field
    for (label k=0; k < nCells+1; k++)
    {
        for (label j=0; j < nCells+1; j++)
        {
            for (label i=0; i < nCells+1; i++)
            {
                points[i+j*yOffset+k*zOffset] = point
                (
                    i/double(nCells),
                    j/double(nCells),
                    k/double(nCells)
                );
            }
        }
    }


    // ========================================================================
    // Construction of the faces
    // ========================================================================

    // First all internal faces
    label ind = 0;

    // Faces along the x direction
    for (label xInd=1; xInd < nCells; xInd++)
    {
        for (label zInd=0; zInd < nCells; zInd++)
        {
            for (label yInd=0; yInd < nCells; yInd++)
            {
                face& f = faces[ind];
                f.resize(4);
                f[0] = xInd +     yInd*yOffset +     zInd*zOffset;
                f[1] = xInd + (yInd+1)*yOffset +     zInd*zOffset;
                f[2] = xInd + (yInd+1)*yOffset + (zInd+1)*zOffset;
                f[3] = xInd +     yInd*yOffset + (zInd+1)*zOffset;
                owner[ind]    = (xInd-1) + yInd*nCells + zInd*std::pow(nCells,2);
                neighbor[ind] = xInd     + yInd*nCells + zInd*std::pow(nCells,2);
                ind++;
            }
        }
    }


    // Faces along the y direction

    for (label yInd=1; yInd < nCells; yInd++)
    {
        for (label zInd=0; zInd < nCells; zInd++)
        {
            for (label xInd=0; xInd < nCells; xInd++)
            {
                face& f = faces[ind];
                f.resize(4);
                f[0] = xInd     + yInd*yOffset +     zInd*zOffset;
                f[1] = xInd     + yInd*yOffset + (zInd+1)*zOffset;
                f[2] = (xInd+1) + yInd*yOffset + (zInd+1)*zOffset;
                f[3] = (xInd+1) + yInd*yOffset +     zInd*zOffset;
                owner[ind]    = xInd + (yInd-1)*nCells + zInd*std::pow(nCells,2);
                neighbor[ind] = xInd +     yInd*nCells + zInd*std::pow(nCells,2);
                ind++;
            }
        }
    }

    // Faces along the z direction
    for (label zInd=1; zInd < nCells; zInd++)
    {
        for (label yInd=0; yInd < nCells; yInd++)
        {
            for (label xInd=0; xInd < nCells; xInd++)
            {
                face& f = faces[ind];
                f.resize(4);
                f[0] = xInd     +     yInd*yOffset + zInd*zOffset;
                f[1] = (xInd+1) +     yInd*yOffset + zInd*zOffset;
                f[2] = (xInd+1) + (yInd+1)*yOffset + zInd*zOffset;
                f[3] = xInd     + (yInd+1)*yOffset + zInd*zOffset;
                owner[ind]    = xInd +     yInd*nCells + (zInd-1)*std::pow(nCells,2);
                neighbor[ind] = xInd +     yInd*nCells + (zInd-1)*std::pow(nCells,2);
                ind++;
            }
        }
    }

    // ========================================================================
    // Construction of the boundary faces
    // ========================================================================

    const label startIndexBoundaryFaces = ind;

    // Faces along the x direction
    label xInd = 0;
    for (label zInd=0; zInd < nCells; zInd++)
    {
        for (label yInd=0; yInd < nCells; yInd++)
        {
            face& f = faces[ind];
            f.resize(4);
            f[0] = xInd +     yInd*yOffset     +     zInd*zOffset;
            f[1] = xInd +     yInd*yOffset     + (zInd+1)*zOffset;
            f[2] = xInd +     (yInd+1)*yOffset + (zInd+1)*zOffset;
            f[3] = xInd +     (yInd+1)*yOffset +     zInd*zOffset;
            owner[ind]    = xInd + yInd*nCells + zInd*std::pow(nCells,2);
            ind++;
        }
    }

    xInd = nCells;
    for (label zInd=0; zInd < nCells; zInd++)
    {
        for (label yInd=0; yInd < nCells; yInd++)
        {
            face& f = faces[ind];
            f.resize(4);
            f[0] = xInd +     yInd*yOffset +     zInd*zOffset;
            f[1] = xInd + (yInd+1)*yOffset +     zInd*zOffset;
            f[2] = xInd + (yInd+1)*yOffset + (zInd+1)*zOffset;
            f[3] = xInd +     yInd*yOffset + (zInd+1)*zOffset;
            owner[ind]    = (xInd-1) + yInd*nCells + zInd*std::pow(nCells,2);
            ind++;
        }
    }


    label yInd = 0;
    for (label zInd=0; zInd < nCells; zInd++)
    {
        for (label xInd=0; xInd < nCells; xInd++)
        {
            face& f = faces[ind];
            f.resize(4);
            f[0] = xInd     + yInd*yOffset +     zInd*zOffset;
            f[1] = (xInd+1) + yInd*yOffset +     zInd*zOffset;
            f[2] = (xInd+1) + yInd*yOffset + (zInd+1)*zOffset;
            f[3] = xInd     + yInd*yOffset + (zInd+1)*zOffset;
            owner[ind]    = xInd + yInd*nCells + zInd*std::pow(nCells,2);
            ind++;
        }
    }

    yInd = nCells;
    for (label zInd=0; zInd < nCells; zInd++)
    {
        for (label xInd=0; xInd < nCells; xInd++)
        {
            face& f = faces[ind];
            f.resize(4);
            f[0] = xInd     + yInd*yOffset +     zInd*zOffset;
            f[1] = xInd     + yInd*yOffset + (zInd+1)*zOffset;
            f[2] = (xInd+1) + yInd*yOffset + (zInd+1)*zOffset;
            f[3] = (xInd+1) + yInd*yOffset +     zInd*zOffset;
            owner[ind]    = xInd + (yInd-1)*nCells + zInd*std::pow(nCells,2);
            ind++;
        }
    }

    label zInd = 0;
    for (label yInd=0; yInd < nCells; yInd++)
    {
        for (label xInd=0; xInd < nCells; xInd++)
        {
            face& f = faces[ind];
            f.resize(4);
            f[0] = xInd     +     yInd*yOffset + zInd*zOffset;
            f[1] = xInd     + (yInd+1)*yOffset + zInd*zOffset;
            f[2] = (xInd+1) + (yInd+1)*yOffset + zInd*zOffset;
            f[3] = (xInd+1) +     yInd*yOffset + zInd*zOffset;
            owner[ind]    = xInd +     yInd*nCells + zInd*std::pow(nCells,2);
            ind++;
        }
    }

    zInd = nCells;
    for (label yInd=0; yInd < nCells; yInd++)
    {
        for (label xInd=0; xInd < nCells; xInd++)
        {
            face& f = faces[ind];
            f.resize(4);
            f[0] = xInd     +     yInd*yOffset + zInd*zOffset;
            f[1] = (xInd+1) +     yInd*yOffset + zInd*zOffset;
            f[2] = (xInd+1) + (yInd+1)*yOffset + zInd*zOffset;
            f[3] = xInd     + (yInd+1)*yOffset + zInd*zOffset;
            owner[ind]    = xInd +     yInd*nCells + (zInd-1)*std::pow(nCells,2);
            ind++;
        }
    }

    mesh_.reset
    (
        new fvMesh
        (
            IOobject
            (
                Foam::fvMesh::defaultRegion,
                runTime_().timeName(),
                runTime_(),
                IOobject::NO_READ
            ),
            std::move(points),
            std::move(faces),
            std::move(owner),
            std::move(neighbor),
            false       // Do not sync
        )
    );

    addBoundaries(nCells,startIndexBoundaryFaces);
}


void Foam::unitMesh::createObjectRegistry(const fileName& rootPath)
{
    IStringStream is
    (
        "deltaT          1;"
        "writeFrequency  1;"
    );
    systemDict_.read(is);
    runTime_.reset(
        new Time
        (
            systemDict_,
            rootPath,
            "case/",
            "system",
            "constant"   
        )
    );
}


void Foam::unitMesh::addBoundaries
(
    const label nCells,
    const label startIndex
)
{
    List<polyPatch*> pList(6,nullptr);

    wordList patchNames = {"left","front","bottom","right","back","top"};

    const label patchSize = nCells*nCells;

    forAll(patchNames,i)
    {
        // Create a list of 6 patches
        pList[i] = new polyPatch
        (
            patchNames[i],
            patchSize,
            startIndex + i*patchSize,
            i,
            mesh().boundaryMesh(),
            "calculated"
        );
    }
    mesh_.ref().addFvPatches(pList);
}


const Foam::fvMesh&
Foam::unitMesh::mesh()
{
    return mesh_();
}


const Foam::Time&
Foam::unitMesh::runTime()
{
    return runTime_();
}
