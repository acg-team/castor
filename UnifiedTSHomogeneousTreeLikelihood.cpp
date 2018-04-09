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
 * @file UnifiedTSHomogeneousTreeLikelihood.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 09 04 2018
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
#include "UnifiedTSHomogeneousTreeLikelihood.hpp"

/*
 * Implementation for the interface  likelihood under tree search engines (all canonical models)
 */

UnifiedTSHomogeneousTreeLikelihood::UnifiedTSHomogeneousTreeLikelihood(const Tree &tree,
                                                                       const SiteContainer &data,
                                                                       TransitionModel *model,
                                                                       DiscreteDistribution *rDist,
                                                                       bool checkRooted,
                                                                       bool verbose,
                                                                       bool usePatterns) :
        RHomogeneousTreeLikelihood(tree, data, model, rDist, checkRooted, verbose, usePatterns) {


}

UnifiedTSHomogeneousTreeLikelihood::UnifiedTSHomogeneousTreeLikelihood(const Tree &tree,
                                                                       TransitionModel *model,
                                                                       DiscreteDistribution *rDist,
                                                                       bool checkRooted,
                                                                       bool verbose,
                                                                       bool usePatterns)
        : RHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose, usePatterns) {


}

UnifiedTSHomogeneousTreeLikelihood::UnifiedTSHomogeneousTreeLikelihood(const RHomogeneousTreeLikelihood &lik) : RHomogeneousTreeLikelihood(lik) {}

UnifiedTSHomogeneousTreeLikelihood::~UnifiedTSHomogeneousTreeLikelihood() {

}


/*
 * Implementation for the interface  likelihood under tree search engines (all mixed models)
 */


UnifiedTSHomogeneousTreeLikelihood_PIP::UnifiedTSHomogeneousTreeLikelihood_PIP(const Tree &tree,
                                                                               TransitionModel *model,
                                                                               DiscreteDistribution *rDist,
                                                                               UtreeBppUtils::treemap *tm,
                                                                               bool checkRooted,
                                                                               bool verbose,
                                                                               bool usePatterns) :
        RHomogeneousTreeLikelihood_PIP(tree, model, rDist, tm, checkRooted, verbose, usePatterns) {


}

UnifiedTSHomogeneousTreeLikelihood_PIP::UnifiedTSHomogeneousTreeLikelihood_PIP(const Tree &tree,
                                                                               const SiteContainer &data,
                                                                               TransitionModel *model,
                                                                               DiscreteDistribution *rDist,
                                                                               UtreeBppUtils::treemap *tm,
                                                                               bool checkRooted,
                                                                               bool verbose,
                                                                               bool usePatterns) :
        RHomogeneousTreeLikelihood_PIP(tree, data, model, rDist, tm, checkRooted, verbose, usePatterns) {


}

UnifiedTSHomogeneousTreeLikelihood_PIP::UnifiedTSHomogeneousTreeLikelihood_PIP(const RHomogeneousTreeLikelihood_PIP &lik) : RHomogeneousTreeLikelihood_PIP(lik) {

}

UnifiedTSHomogeneousTreeLikelihood_PIP::~UnifiedTSHomogeneousTreeLikelihood_PIP() {}
