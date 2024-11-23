/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           |
     \\/     M anipulation  | 2020, Aalto University, Finland
-------------------------------------------------------------------------------
License

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
    Foam::LoadBalanced_pyJacChemistryModel

Description
    Extends LoadBalancedChemistryModel model by adding analytical Jacobian.

SourceFiles
    LoadBalanced_pyJacChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef LoadBalanced_pyJacChemistryModel_H
#define LoadBalanced_pyJacChemistryModel_H

#include "clockTime.H"
#include "LoadBalancedChemistryModel.H"

// pyJac dependencies
extern "C" {
    #include "pyJacInclude/chem_utils.h"
    #include "pyJacInclude/dydt.h"
    #include "pyJacInclude/jacob.h"
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class LoadBalanced_pyJacChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class LoadBalanced_pyJacChemistryModel
: 
	public LoadBalancedChemistryModel<ReactionThermo, ThermoType> 
{
    template <class DeltaTType> scalar solve(const DeltaTType& deltaT);

public:
    //- Runtime type information
    TypeName("LB_pyJac");

    // Constructors

    //- Construct from thermo
//    LoadBalanced_pyJacChemistryModel(const ThermoType& thermo);
	LoadBalanced_pyJacChemistryModel(ReactionThermo& thermo);

    //- Destructor
    virtual ~LoadBalanced_pyJacChemistryModel();

    //- Disallow default bitwise copy construction
    LoadBalanced_pyJacChemistryModel(const LoadBalanced_pyJacChemistryModel&);

    inline virtual label nEqns() const;

    //- Override jacobian and derivatives function from standardChemistryModel
    virtual void jacobian(const scalar        t,
                          const scalarField&  c,
//                          const label         li,
                          scalarField&        dcdt,
                          scalarSquareMatrix& J) const override;

    virtual void derivatives(const scalar       t,
                             const scalarField& c,
//                             const label        li,
                             scalarField&       dcdt) const override;

    //- Override Qdot function from standardChemistryModel for using pyJac's enthalpy of formation
    virtual tmp<volScalarField> Qdot() const override;

private:
    ///
    ///@brief List of enthalpy of formation from pyJac
    ///
    ///
    scalarList sp_enth_form;

    //- Update the reaction rate of cell i
    void
    updateReactionRate(const ChemistrySolution& solution, const label& i);

    //- Return the main chemistry variable
    //  (concentration for standard, mass fraction for pyJac)
    scalarField getVariable(const scalarField& concentration, const scalarField& massFraction);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "LoadBalanced_pyJacChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "LoadBalanced_pyJacChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
