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
 * @file TSHTopologySearch.hpp
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
#ifndef MINIJATI_TSHTOPOLOGYSEARCH_HPP
#define MINIJATI_TSHTOPOLOGYSEARCH_HPP

#include <chrono>
#include <Bpp/Phyl/TopologySearch.h>
#include <Bpp/Phyl/Likelihood/AbstractHomogeneousTreeLikelihood.h>
#include <Utree.hpp>
#include <TreeRearrangment.hpp>

#include "Utilities.hpp"
#include "UnifiedTSHomogeneousTreeLikelihood_Generic.hpp"
#include "UnifiedTSHomogeneousTreeLikelihood_PIP.hpp"

namespace tshlib {

    class TreeSearch {

    private:
        bpp::AbstractHomogeneousTreeLikelihood *likelihoodFunc;
        double initialLikelihoodValue;
        TreeSearchHeuristics tshStrategy;
        TreeRearrangmentOperations tshOperations;
        TreeSearchStopCondition stopConditionMethod;
        double stopConditionValue;
        std::string scoringMethod;
        int performed_moves;
        int search_startingnodes;
        mutable tshlib::Utree *utree_;


    public:


        TreeSearch() {
            likelihoodFunc = nullptr;
            initialLikelihoodValue = -std::numeric_limits<double>::infinity();
            tshStrategy = TreeSearchHeuristics::nosearch;
            tshOperations = TreeRearrangmentOperations::classic_Mixed;
            stopConditionMethod = TreeSearchStopCondition::convergence;
            stopConditionValue = 0;
            performed_moves = 0;
            search_startingnodes = 0;
            scoringMethod = "";
            utree_ = nullptr;
        };

        ~TreeSearch() = default;

        void setUtree(tshlib::Utree *inUTree) {
            utree_ = inUTree;
        }

        void setLikelihoodFunc(bpp::AbstractHomogeneousTreeLikelihood *in_likelihoodFunc) {
            likelihoodFunc = in_likelihoodFunc;

            if (dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)) {
                utree_ = dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)->getUtreeTopology();
            } else {
                utree_ = dynamic_cast<UnifiedTSHomogeneousTreeLikelihood *>(likelihoodFunc)->getUtreeTopology();
            }

        }

        void setInitialLikelihoodValue(double in_initialLikelihoodValue) {
            initialLikelihoodValue = in_initialLikelihoodValue;
        }

        void setTreeSearchStrategy(TreeSearchHeuristics in_tshStrategy, TreeRearrangmentOperations in_tshOperations) {
            tshStrategy = in_tshStrategy;
            tshOperations = in_tshOperations;
        }

        void setStopCondition(TreeSearchStopCondition in_stopConditionMethod, double in_stopConditionValue) {
            stopConditionMethod = in_stopConditionMethod;
            stopConditionValue = in_stopConditionValue;
        }

        bpp::AbstractHomogeneousTreeLikelihood *getLikelihoodFunc() const {
            return likelihoodFunc;
        }

        double getInitialLikelihoodValue() const {
            return initialLikelihoodValue;
        }

        TreeSearchHeuristics getTshStrategy() const {
            return tshStrategy;
        }

        TreeRearrangmentOperations getTshOperations() const {
            return tshOperations;
        }

        TreeSearchStopCondition getStopConditionMethod() const {
            return stopConditionMethod;
        }

        double getStopConditionValue() const {
            return stopConditionValue;
        }

        const std::string getScoringMethod() const {
            return scoringMethod;
        }

        void setScoringMethod(const std::string &inScoringMethod) {
            scoringMethod = inScoringMethod;
        }

        int getStartingNodes() const {
            return search_startingnodes;
        }

        void setStartingNodes(int in_search_startingnodes) {
            search_startingnodes = in_search_startingnodes;
        }


        double performTreeSearch();

    protected:

        tshlib::TreeRearrangment *defineCandidateMoves();

        void testCandidateMoves(tshlib::TreeRearrangment *candidateMoves);

        double greedy();

        double hillclimbing();

        double particleswarming();
    };
}


#endif //MINIJATI_TSHTOPOLOGYSEARCH_HPP
