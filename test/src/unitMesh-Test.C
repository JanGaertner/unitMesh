/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------

License
    This file is part of mmcFoam.

    mmcFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    mmcFoam is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with  mmcFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    Test the unitMesh class

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2023
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared

#include "globalFoamArgs.H"

#include "fvCFD.H"
#include "unitMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


TEST_CASE("unitMesh Test")
{

    SECTION("One Cell Mesh")
    {
        // Call the constructor
        unitMesh uMesh(1);

        const fvMesh& mesh = uMesh.mesh();

        // Check that the number of cells match
        REQUIRE(mesh.C().size()==1);
    }

    SECTION("Multi cell test")
    {
        // Call the constructor
        unitMesh uMesh(10);

        const fvMesh& mesh = uMesh.mesh();

        // Check that the number of cells match
        REQUIRE(mesh.C().size()==std::pow(10,3));
    }

    SECTION("Create a volScalarField on new mesh")
    {
        // Call the constructor
        unitMesh uMesh(10);

        const fvMesh& mesh = uMesh.mesh();

        const Time& runTime = uMesh.runTime();

        volScalarField test
        (
            IOobject
            (
                "test",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("f",dimPressure,1E+5)
        );
        forAll(test,i)
        {
            REQUIRE(test[i]==1E+5);
        }
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

