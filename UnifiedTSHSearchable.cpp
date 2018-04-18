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
 * @file UnifiedTSHSearchable.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 18 04 2018
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
#include <Bpp/Numeric/Function/SimpleMultiDimensions.h>
#include <Bpp/Phyl/Likelihood/PseudoNewtonOptimizer.h>
#include "UnifiedTSHSearchable.hpp"


using namespace bpp;

void UnifiedTSHSearchable::setOptimiser(AbstractHomogeneousTreeLikelihood *lk,
                                                      bool optNumericalDerivatives,
                                                      std::map<std::string, std::string> &params,
                                                      const string &suffix,
                                                      bool suffixIsOptional,
                                                      bool verbose,
                                                      int warn) {
    // -------------------------------------------------------------------------
    // Optimisation algorithm
    std::string optBrLenMethod = ApplicationTools::getStringParameter("optimization.topology.brlen_optimization", params, "Brent", suffix, suffixIsOptional, warn);
    optMethodModel_ = optBrLenMethod;

    // -------------------------------------------------------------------------
    // Message handler
    std::string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional, "none", warn + 1);
    auto *messageHandler = static_cast<OutputStream *>((mhPath == "none") ? 0 : (mhPath == "std") ? ApplicationTools::message : new StlOutputStream(
            new std::ofstream(mhPath.c_str(), std::ios::app)));

    // -------------------------------------------------------------------------
    // Profiler
    std::string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional, "none", warn + 1);
    auto *profiler = static_cast<OutputStream *>((prPath == "none") ? nullptr : (prPath == "std") ? ApplicationTools::message : new StlOutputStream(
            new std::ofstream(prPath.c_str(), std::ios::app)));
    if (profiler) profiler->setPrecision(20);


    auto nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", params, 1000000, suffix, suffixIsOptional, warn + 1);

    if (optNumericalDerivatives) {

        AbstractNumericalDerivative *dn_3points = new ThreePointsNumericalDerivative(lk);
        AbstractNumericalDerivative *dn_5points = new FivePointsNumericalDerivative(lk);

        // Initialise branch optimiser
        if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_BRENT) {
            optimiser_ = new SimpleMultiDimensions(dn_5points);
        } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_BFGS) {
            optimiser_ = new BfgsMultiDimensions(dn_5points);
        } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_NEWTON) {
            optimiser_ = new PseudoNewtonOptimizer(dn_3points);
        } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_GRADIENT) {
            optimiser_ = new ConjugateGradientMultiDimensions(dn_5points);
        }
    } else {
        // Initialise branch optimiser
        if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_BRENT) {
            optimiser_ = new SimpleMultiDimensions(lk);
        } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_BFGS) {
            optimiser_ = new BfgsMultiDimensions(lk);
        } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_NEWTON) {
            optimiser_ = new PseudoNewtonOptimizer(lk);
        } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_GRADIENT) {
            optimiser_ = new ConjugateGradientMultiDimensions(lk);
        }
    }
    optimiser_->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimiser_->setProfiler(profiler);
    optimiser_->setMessageHandler(messageHandler);
    optimiser_->setMaximumNumberOfEvaluations(nbEvalMax);
    optimiser_->setVerbose(1);

}


void UnifiedTSHSearchable::fireBranchOptimisation(AbstractHomogeneousTreeLikelihood *lk, std::vector<bpp::Node *> extractionNodes) {

    ParameterList parameters;
    // For each node involved in the move, get the corrisponding branch parameter (no root)
    for (auto &bnode:extractionNodes) {
        if (bnode->hasFather()) {
            Parameter brLen = lk->getParameter("BrLen" + TextTools::toString(bnode->getId()));
            brLen.setName("BrLen" + TextTools::toString(bnode->getId()));
            parameters.addParameter(brLen);
        }
    }

    // set parameters on the likelihood function (inherited)
    lk->setParameters(parameters);

    // Re-estimate branch length:
    if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_BRENT) {
        auto optimiserInstance = dynamic_cast<SimpleMultiDimensions *>(optimiser_);
        optimiserInstance->init(parameters);
        optimiserInstance->doInit(parameters);
        optimiserInstance->getStopCondition()->setTolerance(0.001);
        optimiserInstance->optimize();

        // set parameters on the likelihood function (inherited)
        lk->setParameters(optimiserInstance->getParameters());

    } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_BFGS) {
        auto optimiserInstance = dynamic_cast<BfgsMultiDimensions *>(optimiser_);
        optimiserInstance->init(parameters);
        optimiserInstance->doInit(parameters);
        optimiserInstance->getStopCondition()->setTolerance(0.001);
        optimiserInstance->optimize();

        // set parameters on the likelihood function (inherited)
        lk->setParameters(optimiserInstance->getParameters());

    } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_NEWTON) {
        auto optimiserInstance = dynamic_cast<PseudoNewtonOptimizer *>(optimiser_);
        optimiserInstance->init(parameters);
        optimiserInstance->doInit(parameters);
        optimiserInstance->getStopCondition()->setTolerance(0.001);
        optimiserInstance->optimize();

        // set parameters on the likelihood function (inherited)
        lk->setParameters(optimiserInstance->getParameters());

    } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_GRADIENT) {
        auto optimiserInstance = dynamic_cast<ConjugateGradientMultiDimensions *>(optimiser_);
        optimiserInstance->init(parameters);
        optimiserInstance->doInit(parameters);
        optimiserInstance->getStopCondition()->setTolerance(0.001);
        optimiserInstance->optimize();

        // set parameters on the likelihood function (inherited)
        lk->setParameters(optimiserInstance->getParameters());
    }


}

