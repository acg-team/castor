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
 * @file UnifiedTSHomogeneousTreeLikelihood_PIP.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 18 04 2018
 * @version 1.0.10
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
#ifndef MINIJATI_UNIFIEDTSHOMOGENEOUSTREELIKELIHOOD_PIP_HPP
#define MINIJATI_UNIFIEDTSHOMOGENEOUSTREELIKELIHOOD_PIP_HPP

#include "RHomogeneousTreeLikelihood_PIP.hpp"
#include "Utilities.hpp"
#include "UnifiedTSHSearchable.hpp"

#ifdef INTELTBB

#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#endif


namespace bpp {

    class UnifiedTSHomogeneousTreeLikelihood_PIP : public RHomogeneousTreeLikelihood_PIP, public UnifiedTSHSearchable {

    protected:

        mutable tshlib::Utree *utree_;


    public:

        UnifiedTSHomogeneousTreeLikelihood_PIP(const Tree &tree,
                                               TransitionModel *model,
                                               DiscreteDistribution *rDist,
                                               tshlib::Utree *utree_,
                                               UtreeBppUtils::treemap *tm,
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
                                               UtreeBppUtils::treemap *tm,
                                               bool optNumericalDerivatives,
                                               std::map<std::string, std::string> &params,
                                               const std::string &suffix,
                                               bool checkRooted,
                                               bool verbose,
                                               bool usePatterns);


        ~UnifiedTSHomogeneousTreeLikelihood_PIP();

        UnifiedTSHomogeneousTreeLikelihood_PIP *clone() const override { return new UnifiedTSHomogeneousTreeLikelihood_PIP(*this); }

        const Tree &getTopology() const { return getTree(); }

        double getTopologyValue() const throw(Exception) { return getValue(); }

        tshlib::Utree *getUtreeTopology() { return utree_; }

        void init_(bool usePatterns);

        void fireTopologyChange(std::vector<int> nodeList);

        double updateLikelihoodOnTreeRearrangement(std::vector<tshlib::VirtualNode *> &nodeList);

        double getLogLikelihoodOnTreeRearrangement() const;

        void topologyChangeSuccessful(std::vector<tshlib::VirtualNode *> listNodes);

        void topologyCommitTree();

        void addTestLikelihoodData(int idxThread) override;

        void removeTestLikelihoodData(int idxThread) override;


#ifdef INTELTBB
        void recomputeSiteLikelihoodUsingPartitions(const tbb::blocked_range<size_t> &range, std::vector<double> *lk_sites) const;
#endif
    };


#ifdef INTELTBB

    class PartitionedSiteLikelihood {

        friend class UnifiedTSHomogeneousTreeLikelihood_PIP;

    private:

        std::vector<double> *siteLogLK_;
        UnifiedTSHomogeneousTreeLikelihood_PIP *lkFunc_;

    public:
        PartitionedSiteLikelihood(std::vector<double> *inSiteLogLK, UnifiedTSHomogeneousTreeLikelihood_PIP *inLkFunc) :
                siteLogLK_(inSiteLogLK), lkFunc_(inLkFunc) {};

        void operator()(const tbb::blocked_range<size_t> &range) const {

            for (size_t i = range.begin(); i != range.end(); ++i) {

                std::vector<int> tempExtendedNodeList;
                const std::vector<unsigned int> *rootWeights = &lkFunc_->likelihoodData_->getWeights();

                // Extend it
                lkFunc_->_extendNodeListOnSetA(lkFunc_->likelihoodNodes_.back(), tempExtendedNodeList, i);

                // call to function which retrieves the lk value for each site
                (*siteLogLK_)[i] = = log(lkFunc_->computeLikelihoodForASite(tempExtendedNodeList, i)) * rootWeights->at(i);
            }

        };


    };

#endif

}


#endif //MINIJATI_UNIFIEDTSHOMOGENEOUSTREELIKELIHOOD_PIP_HPP
