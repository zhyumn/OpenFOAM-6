/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "ISAT.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//const dataType Foam::ISAT::staticData();


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//Foam::ISAT::ISAT()
//:
    //baseClassName(),
    //data_()
//{}

Foam::ISAT::ISAT(SUPstream::mpi_manager& manager_in,label nmem)
:nmem_leaf_(nmem),nmem_node_(nmem),leafmem(manager_in,nmem_leaf_),
emptylist_leaf(manager_in,nmem_leaf_),nodemem(manager_in,nmem_node_),
emptylist_node(manager_in,nmem_node_),pv_(manager_in),
size_leaf(pv_().size_leaf),head_leaf(pv_().head_leaf),
tail_leaf(pv_().tail_leaf),temp_tail_leaf(pv_().temp_tail_leaf),
size_node(pv_().size_node),heaf_node(pv_().heaf_node),
tail_node(pv_().tail_node),temp_tail_node(pv_().temp_tail_node)
{
    for(label i=0;i<nmem_leaf_;i++)
    {
        emptylist_leaf[i]=i;
    }
    for(label i=0;i<nmem_node_;i++)
    {
        emptylist_node[i]=i;
    }
    size_leaf=0;
    head_leaf=0;
    tail_leaf=nmem_leaf_-1;
    temp_tail_leaf=nmem_leaf_-1;

    size_node=0;
    heaf_node=0;
    tail_node=nmem_node_-1;
    temp_tail_node=nmem_node_-1;


    tleaf.node_=&tnode;
    tleaf.v=1;
    tnode.leaf_=&tleaf;
    tnode.v=2; 

}
/*
Foam::ISAT::ISAT(const dataType& data)
:
    baseClassName(),
    data_(data)
{}


Foam::ISAT::ISAT(const ISAT&)
:
    baseClassName(),
    data_()
{}
*/

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
/*
Foam::autoPtr<Foam::ISAT>
Foam::ISAT::New()
{
    return autoPtr<ISAT>(new ISAT);
}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ISAT::~ISAT()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //
/*
void Foam::ISAT::operator=(const ISAT& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}
*/
// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
