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
 * @file UnifiedTSHSearchable.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 18 04 2018
 * @version 1.0.7
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
 * @see For more information visit: https://bitbucket.org/acg-team/minijati/wiki/Home
 */
#ifndef MINIJATI_UNIFIEDTSHSEARCHABLE_HPP
#define MINIJATI_UNIFIEDTSHSEARCHABLE_HPP

#include <Bpp/Numeric/Function/AbstractNumericalDerivative.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/BfgsMultiDimensions.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/FivePointsNumericalDerivative.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/Likelihood/AbstractHomogeneousTreeLikelihood.h>

namespace bpp {

    class UnifiedTSHSearchable {

    protected:

        mutable AbstractHomogeneousTreeLikelihood *lk_;
        mutable std::vector<Parameter> testBrLens_;
        mutable AbstractOptimizer *optimiser_;
        std::string optMethodModel_;
        bool optNumericalDerivatives_;


    public:

        std::string OPTIMIZATION_NEWTON = "newton";
        std::string OPTIMIZATION_GRADIENT = "gradient";
        std::string OPTIMIZATION_BRENT = "Brent";
        std::string OPTIMIZATION_BFGS = "BFGS";

        UnifiedTSHSearchable() : lk_(nullptr), optimiser_(nullptr), optMethodModel_(""), optNumericalDerivatives_(false) {}

        virtual ~UnifiedTSHSearchable() = default;

        void setOptimiser(AbstractHomogeneousTreeLikelihood *lk,
                                        bool optNumericalDerivatives,
                                        std::map<std::string, std::string> &params,
                                        const string &suffix,
                                        bool suffixIsOptional,
                                        bool verbose,
                                        int warn);

        void fireBranchOptimisation(std::vector<bpp::Node *> extractionNodes);


    };
}


#endif //MINIJATI_UNIFIEDTSHSEARCHABLE_HPP
