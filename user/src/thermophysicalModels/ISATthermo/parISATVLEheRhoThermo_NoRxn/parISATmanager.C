/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "parISATmanager.H"
#include "LUscalarMatrix.H"

#define REUSELIST

#define LIKELY(exp) __builtin_expect(exp, 1)
#define UNLIKELY(exp) __builtin_expect(exp, 0)

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template <class FuncType>
Foam::ISATmanager<FuncType>::ISATmanager(label in_n, label out_n, FuncType &func, const word &name_in, const dictionary &ISATDict)
    : Ninput(in_n), Noutput(out_n), ISATDict_(ISATDict.subDict(name_in)),
      tableTree_(in_n, out_n, ISATDict_), // readLabel(ISATDict_.lookup("maxNLeafs")), readLabel(ISATDict_.lookup("NtimeTag"))),
      pfunc(&func),
      epsilon_(1e-3),
      relepsilon_(0.0),
      scaleFactor_(out_n, out_n),
      scaleIn_(in_n, in_n),
      toleranceOut2_(out_n, out_n),
      initToleranceIn2_(in_n, in_n),
      init_elp_(in_n, in_n),
      modified_(false),
      timeSteps_(0),
      noISAT_(ISATDict_.lookupOrDefault<bool>("noISAT", false)),
      checkInterval_(readLabel(ISATDict_.lookup("checkInterval"))),
      maxDepthFactor_(readScalar(ISATDict_.lookup("maxDepthFactor"))),
      maxLeafsize_(in_n),
      nRetrieved_(0),
      nGrowth_(0),
      nAdd_(0),
      nCall_(0),
      nRRetrieved_(0),
      nRGrowth_(0),
      nRAdd_(0),
      nRCall_(0),
      treename_(name_in),
      muted_(ISATDict_.lookupOrDefault<bool>("muted", false))
{
    for (int i = 0; i < out_n; i++)
        for (int j = 0; j < out_n; j++)
        {
            scaleFactor_[i][j] = 0;
            toleranceOut2_[i][j] = 0;
        }
    for (int i = 0; i < in_n; i++)
        for (int j = 0; j < in_n; j++)
        {
            init_elp_[i][j] = 0;
            scaleIn_[i][j] = 0;
            initToleranceIn2_[i][j] = 0;
        }
    for (int i = 0; i < out_n; i++)
    {
        scaleFactor_[i][i] = 1.0;
    }
    for (int i = 0; i < in_n; i++)
        scaleIn_[i][i] = 1.0;

    scalarList toleranceOut_temp(out_n);
    scalarList initToleranceIn_temp(in_n);
    scalarList scaleIn_temp(in_n);
    FuncType::read(ISATDict_, maxLeafsize_, toleranceOut_temp, initToleranceIn_temp, scaleIn_temp);
    for (int i = 0; i < out_n; i++)
    {
        toleranceOut2_[i][i] = 1.0 / sqr(toleranceOut_temp[i]);
    }
    for (int i = 0; i < in_n; i++)
    {
        initToleranceIn2_[i][i] = 1.0 / sqr(initToleranceIn_temp[i]);
        scaleIn_[i][i] = scaleIn_temp[i];
    }
}

template <class FuncType>
Foam::ISATmanager<FuncType>::~ISATmanager()
{
    if (tableTree_.size() > 0)
        showPerformance();
}
template <class FuncType>
template <class... Args>
void Foam::ISATmanager<FuncType>::add(const scalarList &value, scalarList &out, bool growflag, Args &... arg)
{
    parISATbinaryTree &T = tableTree_;
    std::atomic<size_t> *root_atomic = reinterpret_cast<std::atomic<long unsigned int> *>(&T.root_.offset);

    static scalarRectangularMatrix A_tmp(Ninput, Noutput);
    static scalarRectangularMatrix EOA_tmp(Ninput, Ninput);

    SharedPointer<parISATleaf> pleaf;
    /*     while (T.size() >= T.maxNLeafs())
    {
        T.deleteLeaf(T.timeTagList().pop());
    } */
    if (T.size() >= T.maxNLeafs())
    {
        return;
    }
    //pleaf = T.insertNewLeaf(value, out);
    if (T.size() == 0)
    {
        T.write_lock.lock();
        if (T.size() == 0)
        {

#ifdef REUSELIST
            SharedPointer<parISATNode> pnode = T.node_manager.New_reuse();
            SharedPointer<parISATleaf> pleaf = T.leaf_manager.New_reuse();
            if (UNLIKELY(pnode.isNULL() || pleaf.isNULL()))
            {
                T.mem_lock.lock();
                if (pnode.isNULL())
                    pnode = T.node_manager.New();
                if (pleaf.isNULL())
                    pleaf = T.leaf_manager.New();
                T.mem_lock.unlock();
            }
#else

            SharedPointer<parISATNode> pnode = T.node_manager.New();
            SharedPointer<parISATleaf> pleaf = T.leaf_manager.New();
#endif

            pnode->parent_.setNULL();
            pnode->nodeLeft_.setNULL();
            pnode->nodeRight_.setNULL();
            pnode->leafLeft_ = pleaf;
            pnode->leafRight_.setNULL();
            //pnode->set();

            pleaf->node_ = pnode;
            pleaf->set(value, pnode, out);
            pleaf->lastUsed = timeSteps_;
            pfunc->derive(value, out, pleaf->A_, arg...);

            for (int i = 0; i < Ninput; i++)
                for (int j = 0; j < Noutput; j++)
                {
                    A_tmp[i][j] = pleaf->A_(i, j);
                }
            EOA_tmp = scaleIn_ * (initToleranceIn2_ + (A_tmp)*toleranceOut2_ * (A_tmp.T())) * scaleIn_;
            //Pout << "!!!!!!~~~" << A_tmp << "," << Ninput << "," << Noutput << endl;
            //auto xx = A_tmp * toleranceOut2_;
            for (int i = 0; i < Ninput; i++)
                for (int j = 0; j < Ninput; j++)
                {
                    pleaf->EOA_(i, j) = EOA_tmp[i][j];
                }
            root_atomic->store(pnode.offset, std::memory_order_release);
            T.size_leaf_++;
            T.size_node++;
            modified_ = true;
            nRAdd_++;
            T.NAdd++;

            T.write_lock.unlock();
            return;
        }
        else
        {
            //T.NF1++;
            T.write_lock.unlock();
        }
    }
    SharedPointer<parISATNode> parentNode;
    int ret_side; // 0 left, 1 right
    SharedPointer<parISATleaf> pleaf2;

    //Pout << "!!!!!!herexxxzzz0" << endl;
    T.binaryTreeSearch(value, SharedPointer<parISATNode>(root_atomic->load(std::memory_order_acquire)), pleaf2, parentNode, ret_side);

    if (pleaf2.isNULL() || pleaf2->inEOA(value, scaleIn_)) //pleaf2.isNULL() is necessary
    {
        //T.NF2++;
        return;
    }
    //Pout << "!!!!!!herexxxzzz1" << endl;
    if (growflag == true && (*parentNode).mutex.try_lock())
    {
        //Pout << "!!!!!!herexxxzzz1==" << T.size() << endl;
        bool tmp_flag = false;
        if (ret_side == 0)
        {
            tmp_flag = (parentNode->leafLeft_ == pleaf2 && pleaf2->node_ == parentNode);
        }
        else
        {
            tmp_flag = (parentNode->leafRight_ == pleaf2 && pleaf2->node_ == parentNode);
        }
        /*         FatalErrorInFunction
            << " test!!!!!1\n"
            << ret_side << "\n"
            << (parentNode->leafLeft_ == pleaf2 && pleaf2->node_ == parentNode) << "\n"
            << (parentNode->leafLeft_ == pleaf2) << "\n"
            << (pleaf2->node_ == parentNode) << "\n"

            << T.root_.offset << "\n"
            << parentNode.offset << "\n"
            << parentNode->leafLeft_.offset << "\n"
            << pleaf2.offset << "\n"
            << pleaf2->node_.offset << "\n"

            << (parentNode->leafRight_ == pleaf2 && pleaf2->node_ == parentNode) << "\n"
            << "size=" << T.size()
            << exit(FatalError); */
        //Pout << "!!!!!!herexxxzzz2" << endl;
        if (tmp_flag)
        {

            static scalarList dvalue(value.size(), Zero);
            for (int i = 0; i < tableTree_.n_in_; i++)
            {
                dvalue[i] = value[i] - pleaf2->value_[i];
            }
            //Pout << "!!!!!!herexxxzzz3" << endl;
            //FatalErrorInFunction<<" test!!!!!xz\n"<< exit(FatalError);
            if (grow2(pleaf2, dvalue, out))
            {
                pleaf2->lastUsed = timeSteps_;
                /*                 FatalErrorInFunction << " test!!!!!xzz\n"
                                     << exit(FatalError); */
                //Pout << "!!!!!!herexxxzzz4" << endl;
                (*parentNode).mutex.unlock();
                return;
            }
            else
            {
                //T.NF5++;
            }
        }
        else
        {
            //T.NF4++;
        }
        (*parentNode).mutex.unlock();
    }
    else
    {
        //T.NF3++;
    }
    //Pout << "!!!!!!herexxxzzz1" << endl;
#ifdef REUSELIST
    SharedPointer<parISATNode> pnode_new = T.node_manager.New_reuse();
    SharedPointer<parISATleaf> pleaf_new = T.leaf_manager.New_reuse();
    if (UNLIKELY(pnode_new.isNULL() || pleaf_new.isNULL()))
    {
        T.mem_lock.lock();
        if (pnode_new.isNULL())
            pnode_new = T.node_manager.New();
        if (pleaf_new.isNULL())
            pleaf_new = T.leaf_manager.New();
        T.mem_lock.unlock();
    }
#else

    T.mem_lock.lock();
    SharedPointer<parISATNode> pnode_new = T.node_manager.New();
    SharedPointer<parISATleaf> pleaf_new = T.leaf_manager.New();
    T.mem_lock.unlock();
#endif

    pleaf_new->set(value, pnode_new, out);
    pleaf_new->lastUsed = timeSteps_;
    //pleaf_new->node_2 = pnode_new;
    pfunc->derive(value, out, pleaf_new->A_, arg...);
    for (int i = 0; i < Ninput; i++)
        for (int j = 0; j < Noutput; j++)
        {
            A_tmp[i][j] = pleaf_new->A_(i, j);
        }
    EOA_tmp = scaleIn_ * (initToleranceIn2_ + (A_tmp)*toleranceOut2_ * (A_tmp.T())) * scaleIn_;
    for (int i = 0; i < Ninput; i++)
        for (int j = 0; j < Ninput; j++)
        {
            pleaf_new->EOA_(i, j) = EOA_tmp[i][j];
        }

    pnode_new->parent_.setNULL();
    pnode_new->nodeLeft_.setNULL();
    pnode_new->nodeRight_.setNULL();
    pnode_new->set(pleaf_new, pleaf2); //v = (leafmem[new_leaf_index].v + leafmem[leaf_index].v) / 2;
    //Pout << "!!!!!!herexxxzzz2" << endl;

    if (!(*parentNode).mutex.try_lock())
    {
        //T.NF6++;
        // T.mem_lock.lock();
        T.node_manager.Delete(pnode_new);
        T.leaf_manager.Delete(pleaf_new);
        //T.mem_lock.unlock();
        return;
    }
    /*     if (T.size() >= 4)
    {
          FatalErrorInFunction << " test!!!!!xzzAAAAAAAA\n"
                             << exit(FatalError); 

        FatalErrorInFunction
            << " test!!!!!1\n"
            << "ret_side=" << ret_side << "\n"
            << (parentNode->leafLeft_ == pleaf2 && pleaf2->node_ == parentNode) << "\n"
            << (parentNode->leafLeft_ == pleaf2) << "\n"
            << (pleaf2->node_ == parentNode) << "\n"

            << T.root_.offset << "\n"
            << parentNode.offset << "\n"
            << parentNode->leafLeft_.offset << "\n"
            << pleaf2.offset << "\n"
            << pleaf2->node_.offset << "\n"

            << (parentNode->leafRight_ == pleaf2 && pleaf2->node_ == parentNode) << " !!\n"
            << parentNode->leafRight_.offset << "\n"
            << "size=" << T.size()
            << exit(FatalError);
    } */
    bool tmp_flag = false;
    if (ret_side == 0)
    {
        tmp_flag = (parentNode->leafLeft_ == pleaf2 && pleaf2->node_ == parentNode);
    }
    else
    {
        tmp_flag = (parentNode->leafRight_ == pleaf2 && pleaf2->node_ == parentNode);
    }
    if (!tmp_flag)
    {
        //T.NF7++;
        //T.mem_lock.lock();
        T.node_manager.Delete(pnode_new);
        T.leaf_manager.Delete(pleaf_new);
        //T.mem_lock.unlock();
        (*parentNode).mutex.unlock();
        return;
    }

    /*     FatalErrorInFunction
        << " test!!!!!2"
        << exit(FatalError); */
    //Pout << "!!!!!!herexxxzzz3" << endl;
    pleaf2->node_ = pnode_new;

    if (pnode_new->goLeft(value))
    {
        pnode_new->leafLeft_ = pleaf_new;
        pnode_new->leafRight_ = pleaf2;
    }
    else
    {
        pnode_new->leafLeft_ = pleaf2;
        pnode_new->leafRight_ = pleaf_new;
    }

    if (T.size_leaf_ > 1)
    {
        pnode_new->parent_ = parentNode;
        //SharedPointer<parISATNode> pnode_new = node_manager.New();
        //SharedPointer<parISATleaf> pleaf_new = leaf_manager.New();
        SharedPointer<parISATleaf> *leaf;
        SharedPointer<parISATNode> *node;
        if (parentNode->leafLeft_ == pleaf2)
        {
            leaf = &(parentNode->leafLeft_);
            node = &(parentNode->nodeLeft_);
        }
        else
        {
            leaf = &(parentNode->leafRight_);
            node = &(parentNode->nodeRight_);
        }

        //std::atomic<size_t> *node_atomic = (std::atomic<long unsigned int> *)(&node->offset);
        //std::atomic<size_t> *leaf_atomic = (std::atomic<long unsigned int> *)(&leaf->offset);

        std::atomic<size_t> *node_atomic = reinterpret_cast<std::atomic<long unsigned int> *>(&node->offset);
        std::atomic<size_t> *leaf_atomic = reinterpret_cast<std::atomic<long unsigned int> *>(&leaf->offset);

        node_atomic->store(pnode_new.offset, std::memory_order_release);
        leaf_atomic->store(sptr_NULL, std::memory_order_release);

        T.size_leaf_++;
        T.size_node++;
        nRAdd_++;
        T.NAdd++;
        (*parentNode).mutex.unlock();
        return;
    }
    else
    {
        //std::atomic<size_t> *root_atomic = (std::atomic<long unsigned int> *)(&T.root_.offset);
        std::atomic<size_t> *root_atomic = reinterpret_cast<std::atomic<long unsigned int> *>(&T.root_.offset);

        root_atomic->store(pnode_new.offset, std::memory_order_release);

        //T.mem_lock.lock();
        T.node_manager.Delete(parentNode);
        //T.mem_lock.unlock();
        T.size_leaf_++;
        nRAdd_++;
        T.NAdd++;

        return;
    }
}
template <class FuncType>
Foam::SharedPointer<parISATleaf> Foam::ISATmanager<FuncType>::search(const scalarList &value)
{
    SharedPointer<parISATleaf> p;
    tableTree_.binaryTreeSearch(value, tableTree_.root_, p);
    return p;
}
template <class FuncType>
template <class... Args>
bool Foam::ISATmanager<FuncType>::call(
    const Foam::scalarList &value, scalarList &out, Args &... arg)
{
    SharedPointer<parISATleaf> pleaf;
    if (noISAT_)
    {
        pfunc->value(value, out, arg...);
        return true;
    }
    if (nloop_ > 0)
    {
        return retrieve(value, out);
    }

    if (!retrieve(value, out))
    {
        pleaf = search(value);

        if (pleaf.isNULL())
        {
            pfunc->value(value, out, arg...);
            add(value, out, false, arg...);
        }
        else
        {

            //pfunc->value(leafvalue - dvalue, out2, arg...);
            static scalarList dvalue(value.size(), Zero);
            for (int i = 0; i < tableTree_.n_in_; i++)
            {
                dvalue[i] = value[i] - pleaf->value_[i];
            }

            bool flag = true;
            for (int i = 0; i < tableTree_.n_in_; i++)
            {
                flag = flag && mag(dvalue[i]) <= maxLeafsize_[i];
            }
            pfunc->value(value, out, arg...);

            //if (!flag || !grow2(pleaf, dvalue, out))
            //    add(value, out, arg...);
            //if (!flag)
            add(value, out, flag, arg...);
        }
    }
    else
    {
        tableTree_.NRetrieved++;
    }
    tableTree_.NCall++;
    nRCall_++;
    return true;
}

template <class FuncType>
void Foam::ISATmanager<FuncType>::tablevalue(
    const Foam::scalarList &value, scalarList &out)
{
    SharedPointer<parISATleaf> pleaf;
    pleaf = search(value);
    pleaf->eval(value, out);
}
template <class FuncType>
bool Foam::ISATmanager<FuncType>::retrieve(
    const Foam::scalarList &value, scalarList &out)
{
    bool retrieved(false);
    SharedPointer<parISATleaf> plf;

    // If the tree is not empty
    if (LIKELY(tableTree_.notempty()))
    {
        plf = search(value);
        if (UNLIKELY(plf.isNULL()))
        {
            return false;
        }

        // lastSearch keeps track of the chemPoint we obtain by the regular
        // binary tree search
        //lastSearch_ = phi0;
        if (plf->inEOA(value, scaleIn_))
        {
            //retrieved = true;
            nRRetrieved_++;
            plf->eval(value, out);
            //plf->lastUsed = timeSteps_;
            plf->increaseNumRetrieve();
            return true;
        }
    }

    return false;
}

/* template <class FuncType>
bool Foam::ISATmanager<FuncType>::grow(
    SharedPointer<parISATleaf> plf,
    const scalarList &dvalue,
    scalarList &data1,
    const scalarList &data2)
{
    // If the pointer to the chemPoint is nullptr, the function stops
    if (!plf)
    {
        return false;
    }

    // Raise a flag when the chemPoint used has been grown more than the
    // allowed number of time

    // If the solution RphiQ is still within the tolerance we try to grow it
    // in some cases this might result in a failure and the grow function of
    // the chemPoint returns false
    scalarList ret1(data1.size()), ret2(data2.size()); //, ret2(data.size());
    scalarList dvalue_m(plf->value_.size());
    for (int i = 0; i < dvalue_m.size(); i++)
    {
        dvalue_m[i] = dvalue[i];
    }

    plf->eval(plf->value_ + dvalue_m, ret1);
    plf->eval(plf->value_ - dvalue_m, ret2);

    if (normalized_distance(ret1, data1) <= 1.0 && normalized_distance(ret2, data2) <= 1.0)
    {
        data1 = ret1;
        plf->grow(plf->value_ + dvalue_m * (1 + 1e-3), scaleIn_);
        tableTree_.timeTagList().renew(plf->pTimeTagList());
        nRGrowth_++;
        modified_ = true;
        return true;
    }
    // The actual solution and the approximation given by ISAT are too different
    else
    {
        return false;
    }
} */

template <class FuncType>
bool Foam::ISATmanager<FuncType>::grow2(
    SharedPointer<parISATleaf> plf,
    const scalarList &dvalue,
    scalarList &data1)
{
    //bool success = false;
    // If the pointer to the chemPoint is nullptr, the function stops
    if (plf.isNULL())
    {
        return false;
    }

    // If the solution RphiQ is still within the tolerance we try to grow it
    // in some cases this might result in a failure and the grow function of
    // the chemPoint returns false
    static scalarList ret1(data1.size()); //, ret2(data.size());
    static scalarList dvalue_m(plf->value_.size());
    static scalarList value_tmp(plf->value_.size());
    for (int i = 0; i < dvalue_m.size(); i++)
    {
        dvalue_m[i] = dvalue[i];
        value_tmp[i] = plf->value_[i];
    }

    //plf->eval(value2, ret2);

    plf->eval(value_tmp + dvalue_m, ret1);
    //plf->eval(plf->value_ - dvalue_m, ret2);

    //Pout << "!!!!!!!!here0" << endl;
    /*     FatalErrorInFunction << " ret1=" << ret1 << "\ndata1=" << data1
                         << normalized_distance(ret1, data1)
                         << exit(FatalError); */
    if (normalized_distance(ret1, data1) <= 1.0)
    {
        /*         FatalErrorInFunction << " test!!!!!xzzA\n"
                             << exit(FatalError); */
        data1 = ret1;
#ifdef REUSELIST
        SharedPointer<parISATleaf> pleaf_new = tableTree_.leaf_manager.New_reuse();
        if (UNLIKELY(pleaf_new.isNULL()))
        {
            tableTree_.mem_lock.lock();
            pleaf_new = tableTree_.leaf_manager.New();
            tableTree_.mem_lock.unlock();
        }
#else
        tableTree_.mem_lock.lock();
        SharedPointer<parISATleaf> pleaf_new = tableTree_.leaf_manager.New();
        tableTree_.mem_lock.unlock();
#endif
        if (pleaf_new.notNULL())
        {
            pleaf_new->node_ = plf->node_;
            pleaf_new->set(*plf);
            pleaf_new->grow(value_tmp + dvalue_m * (1 + 1e-3), scaleIn_);
            SharedPointer<parISATleaf> *leaf;

            if (plf->node_->leafLeft_ == plf)
            {
                leaf = &plf->node_->leafLeft_;
            }
            else
            {
                leaf = &plf->node_->leafRight_;
            }

            //std::atomic<size_t> *leaf_atomic = (std::atomic<long unsigned int> *)(&leaf->offset);
            std::atomic<size_t> *leaf_atomic = reinterpret_cast<std::atomic<long unsigned int> *>(&leaf->offset);
            leaf_atomic->store(pleaf_new.offset, std::memory_order_release);
            //tableTree_.mem_lock.lock();
            tableTree_.leaf_manager.Delete(plf);
            //tableTree_.mem_lock.unlock();
            nRGrowth_++;
            tableTree_.NGrowth++;
            modified_ = true;
            return true;
        }

        /*         plf->grow(value_tmp + dvalue_m * (1 + 1e-3), scaleIn_);
        nRGrowth_++;
        modified_ = true;
        return true; */

        //tableTree_.timeTagList().renew(plf->pTimeTagList());

        return false;
    }
    // The actual solution and the approximation given by ISAT are too different
    else
    {
        return false;
    }
}
template <class FuncType>
double Foam::ISATmanager<FuncType>::distance(const scalarList &l, const scalarList &r)
{
    double sum = 0;
    for (int i = 0; i < l.size(); i++)
        sum += sqr((l[i] - r[i]) * scaleFactor_[i][i]);
    return sqrt(sum);
}
template <class FuncType>
double Foam::ISATmanager<FuncType>::normalized_distance(const scalarList &l, const scalarList &r)
{
    double sum = 0;
    for (int i = 0; i < l.size(); i++)
        sum += sqr((l[i] - r[i])) * toleranceOut2_[i][i];
    return sqrt(sum);
}
template <class FuncType>
double Foam::ISATmanager<FuncType>::norm(const scalarList &l)
{
    double sum = 0;
    for (int i = 0; i < l.size(); i++)
        sum += sqr(l[i] * scaleFactor_[i][i]);
    return sqrt(sum);
}
template <class FuncType>
void Foam::ISATmanager<FuncType>::showPerformance() const
{

    if (muted_)
        return;
    nCall_ += nRCall_;

    nRetrieved_ += nRRetrieved_;
    nGrowth_ += nRGrowth_;
    nAdd_ += nRAdd_;
    if (tableTree_.size() > 0)
    {
        std::cout << treename_ << ", ISAT performance: nCall = " << nCall_ << ", notCall = " << notCall << ", nRetrieved = " << nRetrieved_ << ", nGrowth = " << nGrowth_ << ", nAdd = " << nAdd_ << std::endl;
        std::cout << treename_ << ", recent info: nCall = " << nRCall_ << ", nRetrieved = " << nRRetrieved_ << ", nGrowth = " << nRGrowth_ << ", nAdd = " << nRAdd_ << std::endl;

        std::cout << "tree size = " << tableTree_.size() << std::endl;
        std::cout << "NtimeSteps: " << timeSteps_ << ", Treedepth: " << tableTree_.depth() << ", Mindepth: " << ceil(log2(tableTree_.size() + 1)) << std::endl;
        std::cout << "maxNLeafs: " << tableTree_.maxNLeafs() << std::endl;
        Info << treename_ << ", shared ISAT performance: nCall = " << tableTree_.NCall << ", nRetrieved = " << tableTree_.NRetrieved << ", nGrowth = " << tableTree_.NGrowth << ", nAdd = " << tableTree_.NAdd << endl;
/*         Info << "NF1 = " << tableTree_.NF1 << ", NF2 = " << tableTree_.NF2 << ", NF3 = " << tableTree_.NF3 << endl;
        Info << "NF4 = " << tableTree_.NF4 << ", NF5 = " << tableTree_.NF5 << ", NF6 = " << tableTree_.NF6 << endl;
        Info << "NF7 = " << tableTree_.NF7 << endl; */
    }
    nRCall_ = 0;
    nRRetrieved_ = 0;
    nRGrowth_ = 0;
    nRAdd_ = 0;
    //if (tableTree_.size() > 0)
    //    tableTree_.timeTagList().print_simple();
}

// ************************************************************************* //
