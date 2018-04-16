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
 * @file UnifiedTSHomogeneousTreeLikelihood.hpp
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
#ifndef MINIJATI_UNIFIEDTSHOMOGENEOUSTREELIKELIHOOD_HPP
#define MINIJATI_UNIFIEDTSHOMOGENEOUSTREELIKELIHOOD_HPP

#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Numeric/Function/AbstractNumericalDerivative.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/BfgsMultiDimensions.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/FivePointsNumericalDerivative.h>
#include "RHomogeneousTreeLikelihood_PIP.hpp"
#include "TSHSearchable.hpp"

namespace bpp {


    class UnifiedTSHSearchable : public virtual Clonable {

    protected:

        mutable std::vector<Parameter> testBrLens_;
        mutable AbstractOptimizer *optimiser_;
        std::string optMethodModel_;
        bool optNumericalDerivatives_;


    public:
        static std::string OPTIMIZATION_GRADIENT;
        static std::string OPTIMIZATION_NEWTON;
        static std::string OPTIMIZATION_BRENT;
        static std::string OPTIMIZATION_BFGS;


        UnifiedTSHSearchable(bool optNumericalDerivatives,
                             std::map<std::string, std::string> &params,
                             const std::string &suffix,
                             bool suffixIsOptional,
                             bool verbose,
                             int warn) {}


        virtual ~UnifiedTSHSearchable() = default;


        virtual UnifiedTSHSearchable *clone() const = 0;


        void setOptimiser(AbstractHomogeneousTreeLikelihood *lk,
                          bool optNumericalDerivatives,
                          std::map<std::string, std::string> &params,
                          const std::string &suffix,
                          bool suffixIsOptional,
                          bool verbose,
                          int warn);


        void fireBranchOptimisation(AbstractHomogeneousTreeLikelihood *lk, std::vector<bpp::Node *> extractionNodes);


        void commitChangesOnTopology(tshlib::Utree *inUTree) {}


        void topologyChangeTested() {


        }


        void topologyChangeSuccessful(const std::vector<tshlib::VirtualNode *> listNodes, tshlib::Utree *inUTree) {

            // Update BPP tree using the structure in Utree
            //fixTopologyChanges(inUTree);

            // Add virtual root to compute the likelihood
            //inUTree->addVirtualRootNode();

            // remap the virtual nodes to the bpp nodes
            //std::vector<Node *> extractionNodes = UtreeBppUtils::remapNodeLists(listNodes, tree_, getTreeMap());

            // Optimise branches involved in the tree rearrangement
            //fireBranchOptimisation(listNodes);

            // Remove the virtual root to allow for further tree topology improvements
            //inUTree->removeVirtualRootNode();

        }


    };


    class UnifiedTSHomogeneousTreeLikelihood : public RHomogeneousTreeLikelihood, public virtual UnifiedTSHSearchable {

    private:
        //AbstractOptimizer *optimiser_;
        //std::string optMethodModel_;
        UtreeBppUtils::treemap *treemap_;
        mutable tshlib::Utree *utree_;
        mutable DRASRTreeLikelihoodData *likelihoodData_;
        mutable DRASRTreeLikelihoodData *likelihoodDataTest_;

    public:
        UnifiedTSHomogeneousTreeLikelihood(const Tree &tree,
                                           TransitionModel *model,
                                           DiscreteDistribution *rDist,
                                           tshlib::Utree *utree_,
                                           UtreeBppUtils::treemap *treemap_,
                                           bool optNumericalDerivatives,
                                           std::map<std::string, std::string> &params,
                                           const std::string &suffix,
                                           bool checkRooted,
                                           bool verbose,
                                           bool usePatterns);


        UnifiedTSHomogeneousTreeLikelihood(const Tree &tree,
                                           const SiteContainer &data,
                                           TransitionModel *model,
                                           DiscreteDistribution *rDist,
                                           tshlib::Utree *utree_,
                                           UtreeBppUtils::treemap *treemap_,
                                           bool optNumericalDerivatives,
                                           std::map<std::string, std::string> &params,
                                           const std::string &suffix,
                                           bool checkRooted,
                                           bool verbose,
                                           bool usePatterns);

        ~UnifiedTSHomogeneousTreeLikelihood() override;

        UnifiedTSHomogeneousTreeLikelihood *clone() const override { return new UnifiedTSHomogeneousTreeLikelihood(*this); }

        const Tree &getTopology() const { return getTree(); }

        double getTopologyValue() const throw(Exception) { return getValue(); }


        void init_(bool usePatterns);

    };

    class UnifiedTSHomogeneousTreeLikelihood_PIP : public RHomogeneousTreeLikelihood_PIP, public virtual UnifiedTSHSearchable {

    private:

        //AbstractOptimizer *optimiser_;
        //std::string optMethodModel_;
        UtreeBppUtils::treemap *treemap_;
        mutable tshlib::Utree *utree_;
        mutable DRASRTreeLikelihoodData *likelihoodData_;
        mutable DRASRTreeLikelihoodData *likelihoodEmptyData_;
        mutable DRASRTreeLikelihoodData *likelihoodDataTest_;
        mutable DRASRTreeLikelihoodData *likelihoodEmptyDataTest_;


    public:

        UnifiedTSHomogeneousTreeLikelihood_PIP(const Tree &tree,
                                               TransitionModel *model,
                                               DiscreteDistribution *rDist,
                                               tshlib::Utree *utree_,
                                               UtreeBppUtils::treemap *treemap_,
                                               bool optNumericalDerivatives,
                                               std::map<std::string, std::string> &params,
                                               const std::string &suffix,
                                               bool checkRooted,
                                               bool verbose,
                                               bool usePatterns);

        UnifiedTSHomogeneousTreeLikelihood_PIP(const Tree &tree,
                                               const SiteContainer &data,
                                               TransitionModel *model,
                                               DiscreteDistribution *rDist,
                                               tshlib::Utree *utree_,
                                               UtreeBppUtils::treemap *treemap_,
                                               bool optNumericalDerivatives,
                                               std::map<std::string, std::string> &params,
                                               const std::string &suffix,
                                               bool checkRooted,
                                               bool verbose,
                                               bool usePatterns);

        ~UnifiedTSHomogeneousTreeLikelihood_PIP() override;

        UnifiedTSHomogeneousTreeLikelihood_PIP *clone() const override { return new UnifiedTSHomogeneousTreeLikelihood_PIP(*this); }

        const Tree &getTopology() const { return getTree(); }

        double getTopologyValue() const throw(Exception) { return getValue(); }

        void init_(bool usePatterns);

    };
}

#endif //MINIJATI_UNIFIEDTSHOMOGENEOUSTREELIKELIHOOD_HPP
