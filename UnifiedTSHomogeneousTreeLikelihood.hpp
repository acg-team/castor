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
#include "RHomogeneousTreeLikelihood_PIP.hpp"
#include "TSHSearchable.hpp"

namespace bpp {


    class UnifiedTSHResources : public Clonable {


    public:

        UnifiedTSHResources(tshlib::Utree *utree_, UtreeBppUtils::treemap *treemap_);

        virtual ~UnifiedTSHResources();

        UnifiedTSHResources *clone() const override { return new UnifiedTSHResources(*this); }

    protected:

        AbstractOptimizer *optimiser_;
        UtreeBppUtils::treemap *treemap_;
        mutable tshlib::Utree *utree_;
        mutable DRASRTreeLikelihoodData *likelihoodData_;
        std::string optMethodModel_;


    };


    class UnifiedTSHomogeneousTreeLikelihood : public RHomogeneousTreeLikelihood, public virtual TSHSearchable, public UnifiedTSHResources {

    public:
        UnifiedTSHomogeneousTreeLikelihood(const Tree &tree,
                                           TransitionModel *model,
                                           DiscreteDistribution *rDist,
                                           tshlib::Utree *utree_,
                                           UtreeBppUtils::treemap *treemap_,
                                           bool checkRooted,
                                           bool verbose,
                                           bool usePatterns);


        UnifiedTSHomogeneousTreeLikelihood(const Tree &tree,
                                           const SiteContainer &data,
                                           TransitionModel *model,
                                           DiscreteDistribution *rDist,
                                           tshlib::Utree *utree_,
                                           UtreeBppUtils::treemap *treemap_,
                                           bool checkRooted,
                                           bool verbose,
                                           bool usePatterns);

        UnifiedTSHomogeneousTreeLikelihood(const RHomogeneousTreeLikelihood &lik, tshlib::Utree *utree_, UtreeBppUtils::treemap *treemap_);

        ~UnifiedTSHomogeneousTreeLikelihood() override;

        UnifiedTSHomogeneousTreeLikelihood *clone() const override { return new UnifiedTSHomogeneousTreeLikelihood(*this); }

        const Tree &getTopology() const override { return getTree(); }

        double getTopologyValue() const throw(Exception) override { return getValue(); }

        void topologyChangeTested(const TopologyChangeEvent &event) override;

        void topologyChangeSuccessful(const TopologyChangeEvent &event) override;


    };

    class UnifiedTSHomogeneousTreeLikelihood_PIP : public RHomogeneousTreeLikelihood_PIP, public virtual TSHSearchable, public UnifiedTSHResources {

    public:

        UnifiedTSHomogeneousTreeLikelihood_PIP(const Tree &tree,
                                               TransitionModel *model,
                                               DiscreteDistribution *rDist,
                                               tshlib::Utree *utree_,
                                               UtreeBppUtils::treemap *treemap_,
                                               bool checkRooted,
                                               bool verbose,
                                               bool usePatterns);

        UnifiedTSHomogeneousTreeLikelihood_PIP(const Tree &tree,
                                               const SiteContainer &data,
                                               TransitionModel *model,
                                               DiscreteDistribution *rDist,
                                               tshlib::Utree *utree_,
                                               UtreeBppUtils::treemap *treemap_,
                                               bool checkRooted, bool verbose,
                                               bool usePatterns);

        UnifiedTSHomogeneousTreeLikelihood_PIP(const RHomogeneousTreeLikelihood_PIP &lik, tshlib::Utree *utree_, UtreeBppUtils::treemap *treemap_);

        ~UnifiedTSHomogeneousTreeLikelihood_PIP() override;

        UnifiedTSHomogeneousTreeLikelihood_PIP *clone() const override { return new UnifiedTSHomogeneousTreeLikelihood_PIP(*this); }

        const Tree &getTopology() const override { return getTree(); }

        double getTopologyValue() const throw(Exception) override { return getValue(); }

        void topologyChangeTested(const TopologyChangeEvent &event) override;

        void topologyChangeSuccessful(const TopologyChangeEvent &event) override;

    };
}

#endif //MINIJATI_UNIFIEDTSHOMOGENEOUSTREELIKELIHOOD_HPP
