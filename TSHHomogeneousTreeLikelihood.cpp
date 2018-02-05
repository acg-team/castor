/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2018 by Lorenzo Gatti & Massimo Maiolo
 *******************************************************************************
 *
 * This file is part of miniJATI
 *
 * miniJATI is a free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * miniJATI is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with miniJATI. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file TSHLIBHomogeneousTreeLikelihood.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 05 02 2018
 * @version 1.0
 * @maintainer Lorenzo Gatti
 * @email lg@lorenzogatti.me
 * @maintainer Massimo Maiolo
 * @email massimo.maiolo@zhaw.ch
 * @status Development
 *
 * @brief
 * @details
 * @pre
 * @bug
 * @warning
 *
 * @see For more information visit: 
 */
#include "TSHHomogeneousTreeLikelihood.hpp"
#include <Bpp/Numeric/AutoParameter.h>

using namespace bpp;


TSHHomogeneousTreeLikelihood::TSHHomogeneousTreeLikelihood(
        const Tree &tree,
        TransitionModel *model,
        DiscreteDistribution *rDist,
        bool checkRooted,
        bool verbose)
throw(Exception) :
        DRHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
        //brLikFunction_(0),
        brentOptimizer_(0),
        brLenTSHValues_(),
        brLenTSHParams_() {
    brentOptimizer_ = new BrentOneDimension();
    brentOptimizer_->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    brentOptimizer_->setProfiler(0);
    brentOptimizer_->setMessageHandler(0);
    brentOptimizer_->setVerbose(0);
}


TSHHomogeneousTreeLikelihood::TSHHomogeneousTreeLikelihood(
        const Tree &tree,
        const SiteContainer &data,
        TransitionModel *model,
        DiscreteDistribution *rDist,
        bool checkRooted,
        bool verbose)
throw(Exception) :
        DRHomogeneousTreeLikelihood(tree, data, model, rDist, checkRooted, verbose),
        //brLikFunction_(0),
        brentOptimizer_(0),
        brLenTSHValues_(),
        brLenTSHParams_() {
    brentOptimizer_ = new BrentOneDimension();
    brentOptimizer_->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    brentOptimizer_->setProfiler(0);
    brentOptimizer_->setMessageHandler(0);
    brentOptimizer_->setVerbose(0);
    // We have to do this since the DRHomogeneousTreeLikelihood constructor will not call the overloaded setData method:
    //brLikFunction_ = new BranchLikelihood(getLikelihoodData()->getWeights());
}


TSHHomogeneousTreeLikelihood::TSHHomogeneousTreeLikelihood(const TSHHomogeneousTreeLikelihood &lik) :
        DRHomogeneousTreeLikelihood(lik),
        //brLikFunction_(0),
        brentOptimizer_(0),
        brLenTSHValues_(),
        brLenTSHParams_() {
    //brLikFunction_  = dynamic_cast<BranchLikelihood*>(lik.brLikFunction_->clone());
    brentOptimizer_ = dynamic_cast<BrentOneDimension *>(lik.brentOptimizer_->clone());
    brLenTSHValues_ = lik.brLenTSHValues_;
    brLenTSHParams_ = lik.brLenTSHParams_;
}

TSHHomogeneousTreeLikelihood &TSHHomogeneousTreeLikelihood::operator=(const TSHHomogeneousTreeLikelihood &lik) {
    DRHomogeneousTreeLikelihood::operator=(lik);
    //if (brLikFunction_) delete brLikFunction_;
    //brLikFunction_  = dynamic_cast<BranchLikelihood*>(lik.brLikFunction_->clone());
    if (brentOptimizer_) delete brentOptimizer_;
    brentOptimizer_ = dynamic_cast<BrentOneDimension *>(lik.brentOptimizer_->clone());
    brLenTSHValues_ = lik.brLenTSHValues_;
    brLenTSHParams_ = lik.brLenTSHParams_;
    return *this;
}

TSHHomogeneousTreeLikelihood::~TSHHomogeneousTreeLikelihood() {
    //if (brLikFunction_) delete brLikFunction_;
    delete brentOptimizer_;
}


double TSHHomogeneousTreeLikelihood::testTSHRearrangement(int nodeId) const throw(NodeException) {
    return 0;
}

void TSHHomogeneousTreeLikelihood::doTSHRearrangement(int nodeId) throw(NodeException) {


}

