namespace Foam
{

    //- Interpolate field vf according to direction dir
    template<class Type>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
    (
        const GeometricField<Type, fvPatchField, volMesh>& vf,
        const surfaceScalarField& dir,
        const word& reconFieldName = word::null
    )
    {
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf
        (
            fvc::interpolate
            (
                vf,
                dir,
                "reconstruct("
                + (reconFieldName != word::null ? reconFieldName : vf.name())
                + ')'
            )
        );

        GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();

        sf.rename(vf.name() + '_' + dir.name());

        return tsf;
    }



    namespace fvc {

        template<class Type>
        void surfaceIntegrate_doubleflux
        (
            Field<Type>& ivf,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf_pos,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf_neg
        )
        {
            const fvMesh& mesh = ssf_pos.mesh();

            const labelUList& owner = mesh.owner();
            const labelUList& neighbour = mesh.neighbour();

            const Field<Type>& issf_pos = ssf_pos;
            const Field<Type>& issf_neg = ssf_neg;

            forAll(owner, facei)
            {
                ivf[owner[facei]] += issf_pos[facei];
                ivf[neighbour[facei]] -= issf_neg[facei];
            }

            forAll(mesh.boundary(), patchi)
            {
                const labelUList& pFaceCells =
                    mesh.boundary()[patchi].faceCells();

                const fvsPatchField<Type>& pssf_pos = ssf_pos.boundaryField()[patchi];
                //const fvsPatchField<Type>& pssf_neg = ssf_neg.boundaryField()[patchi];

                forAll(mesh.boundary()[patchi], facei)
                {
                    ivf[pFaceCells[facei]] += pssf_pos[facei];
                    //ivf[pFaceCells[facei]] += pssf_neg[facei] / 2;
                }
            }

            ivf /= mesh.Vsc();
        }

        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
            surfaceIntegrate_doubleflux
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf_pos,
                const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf_neg
            )
        {
            const fvMesh& mesh = ssf_pos.mesh();

            tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
            (
                new GeometricField<Type, fvPatchField, volMesh>
                (
                    IOobject
                    (
                        "surfaceIntegrate(" + ssf_pos.name() + ',' + ssf_neg.name() + ')',
                        ssf_pos.instance(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensioned<Type>
                    (
                        "0",
                        ssf_pos.dimensions() / dimVol,
                        Zero
                        ),
                    extrapolatedCalculatedFvPatchField<Type>::typeName
                    )
            );
            GeometricField<Type, fvPatchField, volMesh>& vf = tvf.ref();

            surfaceIntegrate_doubleflux(vf.primitiveFieldRef(), ssf_pos, ssf_neg);
            vf.correctBoundaryConditions();

            return tvf;
        }

        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
            div_doubleflux
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf_pos,
                const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf_neg
            )
        {
            return tmp<GeometricField<Type, fvPatchField, volMesh>>
                (
                    new GeometricField<Type, fvPatchField, volMesh>
                    (
                        "div_doubleflux(" + ssf_pos.name() + ',' + ssf_neg.name() + ')',
                        surfaceIntegrate_doubleflux(ssf_pos, ssf_neg)
                        )
                    );
        }
    }

}
