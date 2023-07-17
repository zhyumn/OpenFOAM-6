namespace Foam
{
    dimensionSet minMod(const dimensionSet& ds1, const dimensionSet& ds2);
}
Foam::dimensionSet Foam::minMod(const dimensionSet& ds1, const dimensionSet& ds2)
{
    if (dimensionSet::debug && ds1 != ds2)
    {
        FatalErrorInFunction
            << "Arguments of minMod have different dimensions" << endl
            << "     dimensions : " << ds1 << " and " << ds2 << endl
            << abort(FatalError);
    }

    return ds1;
}

#define TEMPLATE template<class Type>
#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    Field<ReturnType>& res,                                                    \
    const UList<Type1>& f1,                                                    \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    TFOR_ALL_F_OP_FUNC_F_F                                                     \
    (                                                                          \
        ReturnType, res, =, ::Foam::Func, Type1, f1, Type2, f2                 \
    )                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func                                                    \
(                                                                              \
    const UList<Type1>& f1,                                                    \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    tmp<Field<ReturnType>> tRes(new Field<ReturnType>(f1.size()));             \
    Func(tRes.ref(), f1, f2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func                                                    \
(                                                                              \
    const UList<Type1>& f1,                                                    \
    const tmp<Field<Type2>>& tf2                                               \
)                                                                              \
{                                                                              \
    tmp<Field<ReturnType>> tRes = reuseTmp<ReturnType, Type2>::New(tf2);       \
    Func(tRes.ref(), f1, tf2());                                               \
    tf2.clear();                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func                                                    \
(                                                                              \
    const tmp<Field<Type1>>& tf1,                                              \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    tmp<Field<ReturnType>> tRes = reuseTmp<ReturnType, Type1>::New(tf1);       \
    Func(tRes.ref(), tf1(), f2);                                               \
    tf1.clear();                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func                                                    \
(                                                                              \
    const tmp<Field<Type1>>& tf1,                                              \
    const tmp<Field<Type2>>& tf2                                               \
)                                                                              \
{                                                                              \
    tmp<Field<ReturnType>> tRes =                                              \
        reuseTmpTmp<ReturnType, Type1, Type1, Type2>::New(tf1, tf2);           \
    Func(tRes.ref(), tf1(), tf2());                                            \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    Field<ReturnType>& res,                                                    \
    const Type1& s1,                                                           \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    TFOR_ALL_F_OP_FUNC_S_F                                                     \
    (                                                                          \
        ReturnType, res, =, ::Foam::Func, Type1, s1, Type2, f2                 \
    )                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func                                                    \
(                                                                              \
    const Type1& s1,                                                           \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    tmp<Field<ReturnType>> tRes(new Field<ReturnType>(f2.size()));             \
    Func(tRes.ref(), s1, f2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func                                                    \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<Field<Type2>>& tf2                                               \
)                                                                              \
{                                                                              \
    tmp<Field<ReturnType>> tRes = reuseTmp<ReturnType, Type2>::New(tf2);       \
    Func(tRes.ref(), s1, tf2());                                               \
    tf2.clear();                                                               \
    return tRes;                                                               \
}


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    Field<ReturnType>& res,                                                    \
    const UList<Type1>& f1,                                                    \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    TFOR_ALL_F_OP_FUNC_F_S                                                     \
    (                                                                          \
        ReturnType, res, =, ::Foam::Func, Type1, f1, Type2, s2                 \
    )                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func                                                    \
(                                                                              \
    const UList<Type1>& f1,                                                    \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<Field<ReturnType>> tRes(new Field<ReturnType>(f1.size()));             \
    Func(tRes.ref(), f1, s2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func                                                    \
(                                                                              \
    const tmp<Field<Type1>>& tf1,                                              \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<Field<ReturnType>> tRes = reuseTmp<ReturnType, Type1>::New(tf1);       \
    Func(tRes.ref(), tf1(), s2);                                               \
    tf1.clear();                                                               \
    return tRes;                                                               \
}


#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)

namespace Foam
{
    BINARY_FUNCTION(Type, Type, Type, minMod)
    BINARY_TYPE_FUNCTION(Type, Type, Type, minMod)

}
#undef TEMPLATE
#undef BINARY_FUNCTION
#undef BINARY_TYPE_FUNCTION_SF
#undef BINARY_TYPE_FUNCTION_FS
#undef BINARY_TYPE_FUNCTION

#define TEMPLATE template<template<class> class Field, class Type>

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    FieldField<Field, ReturnType>& f,                                          \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    forAll(f, i)                                                               \
    {                                                                          \
        Func(f[i], f1[i], f2[i]);                                              \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        FieldField<Field, Type1>::NewCalculatedType(f1)                        \
    );                                                                         \
    Func(tRes.ref(), f1, f2);                                                  \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        reuseTmpFieldField<Field, ReturnType, Type2>::New(tf2)                 \
    );                                                                         \
    Func(tRes.ref(), f1, tf2());                                               \
    tf2.clear();                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1)                 \
    );                                                                         \
    Func(tRes.ref(), tf1(), f2);                                               \
    tf1.clear();                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        reuseTmpTmpFieldField<Field, ReturnType, Type1, Type1, Type2>::        \
            New(tf1, tf2)                                                      \
    );                                                                         \
    Func(tRes.ref(), tf1(), tf2());                                            \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    FieldField<Field, ReturnType>& f,                                          \
    const Type1& s,                                                            \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    forAll(f, i)                                                               \
    {                                                                          \
        Func(f[i], s, f2[i]);                                                  \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const Type1& s,                                                            \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        FieldField<Field, Type2>::NewCalculatedType(f2)                        \
    );                                                                         \
    Func(tRes.ref(), s, f2);                                                   \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const Type1& s,                                                            \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        reuseTmpFieldField<Field, ReturnType, Type2>::New(tf2)                 \
    );                                                                         \
    Func(tRes.ref(), s, tf2());                                                \
    tf2.clear();                                                               \
    return tRes;                                                               \
}


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    FieldField<Field, ReturnType>& f,                                          \
    const FieldField<Field, Type1>& f1,                                        \
    const Type2& s                                                             \
)                                                                              \
{                                                                              \
    forAll(f, i)                                                               \
    {                                                                          \
        Func(f[i], f1[i], s);                                                  \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const Type2& s                                                             \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        FieldField<Field, Type1>::NewCalculatedType(f1)                        \
    );                                                                         \
    Func(tRes.ref(), f1, s);                                                   \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const Type2& s                                                             \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1)                 \
    );                                                                         \
    Func(tRes.ref(), tf1(), s);                                                \
    tf1.clear();                                                               \
    return tRes;                                                               \
}


#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)


namespace Foam
{
    BINARY_FUNCTION(Type, Type, Type, minMod)
    BINARY_TYPE_FUNCTION(Type, Type, Type, minMod)

}
#undef TEMPLATE
#undef BINARY_FUNCTION
#undef BINARY_TYPE_FUNCTION_SF
#undef BINARY_TYPE_FUNCTION_FS
#undef BINARY_TYPE_FUNCTION


#define TEMPLATE template<class Type, class GeoMesh>
#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1,                               \
    const DimensionedField<Type2, GeoMesh>& df2                                \
)                                                                              \
{                                                                              \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        new DimensionedField<ReturnType, GeoMesh>                              \
        (                                                                      \
            IOobject                                                           \
            (                                                                  \
                #Func "(" + df1.name() + ',' + df2.name() + ')',               \
                df1.instance(),                                                \
                df1.db()                                                       \
            ),                                                                 \
            df1.mesh(),                                                        \
            Func(df1.dimensions(), df2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), df1.field(), df2.field());                        \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1,                               \
    const tmp<DimensionedField<Type2, GeoMesh>>& tdf2                          \
)                                                                              \
{                                                                              \
    const DimensionedField<Type2, GeoMesh>& df2 = tdf2();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpDimensionedField<ReturnType, Type2, GeoMesh>::New              \
        (                                                                      \
            tdf2,                                                              \
            #Func "(" + df1.name() + ',' + df2.name() + ')',                   \
            Func(df1.dimensions(), df2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), df1.field(), df2.field());                        \
                                                                               \
    tdf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,                         \
    const DimensionedField<Type2, GeoMesh>& df2                                \
)                                                                              \
{                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1 = tdf1();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpDimensionedField<ReturnType, Type1, GeoMesh>::New              \
        (                                                                      \
            tdf1,                                                              \
            #Func "(" + df1.name() + ',' + df2.name() + ')',                   \
            Func(df1.dimensions(), df2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), df1.field(), df2.field());                        \
                                                                               \
    tdf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,                         \
    const tmp<DimensionedField<Type2, GeoMesh>>& tdf2                          \
)                                                                              \
{                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1 = tdf1();                      \
    const DimensionedField<Type2, GeoMesh>& df2 = tdf2();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpTmpDimensionedField                                            \
            <ReturnType, Type1, Type1, Type2, GeoMesh>::New                    \
        (                                                                      \
            tdf1,                                                              \
            tdf2,                                                              \
            #Func "(" + df1.name() + ',' + df2.name() + ')',                   \
            Func(df1.dimensions(), df2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), df1.field(), df2.field());                        \
                                                                               \
    tdf1.clear();                                                              \
    tdf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const DimensionedField<Type2, GeoMesh>& df2                                \
)                                                                              \
{                                                                              \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        new DimensionedField<ReturnType, GeoMesh>                              \
        (                                                                      \
            IOobject                                                           \
            (                                                                  \
                #Func "(" + dt1.name() + ',' + df2.name() + ')',               \
                df2.instance(),                                                \
                df2.db()                                                       \
            ),                                                                 \
            df2.mesh(),                                                        \
            Func(dt1.dimensions(), df2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), dt1.value(), df2.field());                        \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const Type1& t1,                                                           \
    const DimensionedField<Type2, GeoMesh>& df2                                \
)                                                                              \
{                                                                              \
    return Func(dimensioned<Type1>(t1), df2);                                  \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const tmp<DimensionedField<Type2, GeoMesh>>& tdf2                          \
)                                                                              \
{                                                                              \
    const DimensionedField<Type2, GeoMesh>& df2 = tdf2();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpDimensionedField<ReturnType, Type2, GeoMesh>::New              \
        (                                                                      \
            tdf2,                                                              \
            #Func "(" + dt1.name() + ',' + df2.name() + ')',                   \
            Func(dt1.dimensions(), df2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), dt1.value(), df2.field());                        \
                                                                               \
    tdf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const Type1& t1,                                                           \
    const tmp<DimensionedField<Type2, GeoMesh>>& tdf2                          \
)                                                                              \
{                                                                              \
    return Func(dimensioned<Type2>(t1), tdf2);                                 \
}


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1,                               \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        new DimensionedField<ReturnType, GeoMesh>                              \
        (                                                                      \
            IOobject                                                           \
            (                                                                  \
                #Func "(" + df1.name() + ',' + dt2.name() + ')',               \
                df1.instance(),                                                \
                df1.db()                                                       \
            ),                                                                 \
            df1.mesh(),                                                        \
            Func(df1.dimensions(), dt2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), df1.field(), dt2.value());                        \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1,                               \
    const Type2& t2                                                            \
)                                                                              \
{                                                                              \
    return Func(df1, dimensioned<Type2>(t2));                                  \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,                         \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1 = tdf1();                      \
                                                                               \
    tmp<DimensionedField<ReturnType, GeoMesh>> tRes                            \
    (                                                                          \
        reuseTmpDimensionedField<ReturnType, Type1, GeoMesh>::New              \
        (                                                                      \
            tdf1,                                                              \
            #Func "(" + df1.name() + ',' + dt2.name() + ')',                   \
            Func(df1.dimensions(), dt2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Func(tRes.ref().field(), df1.field(), dt2.value());                        \
                                                                               \
    tdf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<DimensionedField<ReturnType, GeoMesh>> Func                                \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,                         \
    const Type2& t2                                                            \
)                                                                              \
{                                                                              \
    return Func(tdf1, dimensioned<Type2>(t2));                                 \
}


#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)

namespace Foam
{

    BINARY_FUNCTION(Type, Type, Type, minMod)
    BINARY_TYPE_FUNCTION(Type, Type, Type, minMod)

}
#undef TEMPLATE
#undef BINARY_FUNCTION
#undef BINARY_TYPE_FUNCTION_SF
#undef BINARY_TYPE_FUNCTION_FS
#undef BINARY_TYPE_FUNCTION



#define TEMPLATE \
    template<class Type, template<class> class PatchField, class GeoMesh>


#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& res,                      \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    Foam::Func                                                                 \
    (                                                                          \
        res.primitiveFieldRef(),                                               \
        gf1.primitiveField(),                                                  \
        gf2.primitiveField()                                                   \
    );                                                                         \
    Foam::Func                                                                 \
    (                                                                          \
        res.boundaryFieldRef(),                                                \
        gf1.boundaryField(),                                                   \
        gf2.boundaryField()                                                    \
    );                                                                         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        new GeometricField<ReturnType, PatchField, GeoMesh>                    \
        (                                                                      \
            IOobject                                                           \
            (                                                                  \
                #Func "(" + gf1.name() + ',' + gf2.name() + ')',               \
                gf1.instance(),                                                \
                gf1.db(),                                                      \
                IOobject::NO_READ,                                             \
                IOobject::NO_WRITE                                             \
            ),                                                                 \
            gf1.mesh(),                                                        \
            Func(gf1.dimensions(), gf2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), gf1, gf2);                                          \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tgf2                \
)                                                                              \
{                                                                              \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2 = tgf2();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpGeometricField<ReturnType, Type2, PatchField, GeoMesh>::New    \
        (                                                                      \
            tgf2,                                                              \
            #Func "(" + gf1.name() + ',' + gf2.name() + ')',                   \
            Func(gf1.dimensions(), gf2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), gf1, gf2);                                          \
                                                                               \
    tgf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,               \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1 = tgf1();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            tgf1,                                                              \
            #Func "(" + gf1.name() + ',' + gf2.name() + ')',                   \
            Func(gf1.dimensions(), gf2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), gf1, gf2);                                          \
                                                                               \
    tgf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,               \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tgf2                \
)                                                                              \
{                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1 = tgf1();            \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2 = tgf2();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpTmpGeometricField                                              \
            <ReturnType, Type1, Type1, Type2, PatchField, GeoMesh>             \
        ::New                                                                  \
        (                                                                      \
            tgf1,                                                              \
            tgf2,                                                              \
            #Func "(" + gf1.name() + ',' + gf2.name() + ')',                   \
            Func(gf1.dimensions(), gf2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), gf1, gf2);                                          \
                                                                               \
    tgf1.clear();                                                              \
    tgf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& res,                      \
    const dimensioned<Type1>& dt1,                                             \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    Foam::Func(res.primitiveFieldRef(), dt1.value(), gf2.primitiveField());    \
    Foam::Func(res.boundaryFieldRef(), dt1.value(), gf2.boundaryField());      \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        new GeometricField<ReturnType, PatchField, GeoMesh>                    \
        (                                                                      \
            IOobject                                                           \
            (                                                                  \
                #Func "(" + dt1.name() + ',' + gf2.name() + ')',               \
                gf2.instance(),                                                \
                gf2.db(),                                                      \
                IOobject::NO_READ,                                             \
                IOobject::NO_WRITE                                             \
            ),                                                                 \
            gf2.mesh(),                                                        \
            Func(dt1.dimensions(), gf2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), dt1, gf2);                                          \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const Type1& t1,                                                           \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                      \
)                                                                              \
{                                                                              \
    return Func(dimensioned<Type1>(t1), gf2);                                  \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tgf2                \
)                                                                              \
{                                                                              \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2 = tgf2();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpGeometricField<ReturnType, Type2, PatchField, GeoMesh>::New    \
        (                                                                      \
            tgf2,                                                              \
            #Func "(" + dt1.name() + gf2.name() + ',' + ')',                   \
            Func(dt1.dimensions(), gf2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), dt1, gf2);                                          \
                                                                               \
    tgf2.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const Type1& t1,                                                           \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tgf2                \
)                                                                              \
{                                                                              \
    return Func(dimensioned<Type1>(t1), tgf2);                                 \
}


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& res,                      \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    Foam::Func(res.primitiveFieldRef(), gf1.primitiveField(), dt2.value());    \
    Foam::Func(res.boundaryFieldRef(), gf1.boundaryField(), dt2.value());      \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        new GeometricField<ReturnType, PatchField, GeoMesh>                    \
        (                                                                      \
            IOobject                                                           \
            (                                                                  \
                #Func "(" + gf1.name() + ',' + dt2.name() + ')',               \
                gf1.instance(),                                                \
                gf1.db(),                                                      \
                IOobject::NO_READ,                                             \
                IOobject::NO_WRITE                                             \
            ),                                                                 \
            gf1.mesh(),                                                        \
            Func(gf1.dimensions(), dt2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), gf1, dt2);                                          \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                     \
    const Type2& t2                                                            \
)                                                                              \
{                                                                              \
    return Func(gf1, dimensioned<Type2>(t2));                                  \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,               \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1 = tgf1();            \
                                                                               \
    tmp<GeometricField<ReturnType, PatchField, GeoMesh>> tRes                  \
    (                                                                          \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            tgf1,                                                              \
            #Func "(" + gf1.name() + ',' + dt2.name() + ')',                   \
            Func(gf1.dimensions(), dt2.dimensions())                           \
        )                                                                      \
    );                                                                         \
                                                                               \
    Foam::Func(tRes.ref(), gf1, dt2);                                          \
                                                                               \
    tgf1.clear();                                                              \
                                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,               \
    const Type2& t2                                                            \
)                                                                              \
{                                                                              \
    return Func(tgf1, dimensioned<Type2>(t2));                                 \
}


#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)


namespace Foam
{
BINARY_FUNCTION(Type, Type, Type, minMod)
BINARY_TYPE_FUNCTION(Type, Type, Type, minMod)
}



#undef TEMPLATE
#undef BINARY_FUNCTION
#undef BINARY_TYPE_FUNCTION_SF
#undef BINARY_TYPE_FUNCTION_FS
#undef BINARY_TYPE_FUNCTION