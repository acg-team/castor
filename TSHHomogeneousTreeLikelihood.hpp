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
 * @file TSHLIBHomogeneousTreeLikelihood.hpp
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
#ifndef MINIJATI_TSHLIBHOMOGENEOUSTREELIKELIHOOD_HPP
#define MINIJATI_TSHLIBHOMOGENEOUSTREELIKELIHOOD_HPP


#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>
#include <Utree.hpp>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>
#include <Bpp/Numeric/Function/BfgsMultiDimensions.h>
#include "TSHSearchable.hpp"
#include "Utilities.hpp"
#include "RHomogeneousTreeLikelihood_PIP.hpp"

namespace bpp {

    /**
     * @brief This class adds support for topology search to the DRHomogeneousTreeLikelihood class.
     */
    class TSHHomogeneousTreeLikelihood : public RHomogeneousTreeLikelihood, public virtual TSHSearchable {
    protected:
        //BranchLikelihood* brLikFunction_;


        /**
         * @brief Optimizer used for testing NNI.
         */
        BrentOneDimension *brentOptimizer_;
        //PowellMultiDimensions *optimiser_;
        BfgsMultiDimensions *optimiser_;
        AbstractHomogeneousTreeLikelihood *likelihoodFunc_;

        /**
         * @brief Hash used for backing up branch lengths when testing NNIs.
         */
        mutable std::map<int, double> brLenTSHValues_;

        ParameterList brLenTSHParams_;

    public:
        /**
         * @brief Build a new NNIHomogeneousTreeLikelihood object.
         *
         * @param tree The tree to use.
         * @param model The substitution model to use.
         * @param rDist The rate across sites distribution to use.
         * @param checkRooted Tell if we have to check for the tree to be unrooted.
         * If true, any rooted tree will be unrooted before likelihood computation.
         * @param verbose Should I display some info?
         * @throw Exception in an error occured.
         */
        TSHHomogeneousTreeLikelihood(AbstractHomogeneousTreeLikelihood *lk,
                                     const SiteContainer &data,
                                     TransitionModel *model,
                                     DiscreteDistribution *rDist)
        throw(Exception);

        /**
         * @brief Copy constructor.
         */
        TSHHomogeneousTreeLikelihood(const TSHHomogeneousTreeLikelihood &lik);

        TSHHomogeneousTreeLikelihood &operator=(const TSHHomogeneousTreeLikelihood &lik);

        virtual ~TSHHomogeneousTreeLikelihood();

        TSHHomogeneousTreeLikelihood *clone() const { return new TSHHomogeneousTreeLikelihood(*this); }

    public:
        void setData(const SiteContainer &sites) throw(Exception) {
            RHomogeneousTreeLikelihood::setData(sites);

            // The following calls are made for interfaces to node likelihood
            //if (brLikFunction_) delete brLikFunction_;
            //brLikFunction_ = new BranchLikelihood(getLikelihoodData()->getWeights());
        }

        AbstractHomogeneousTreeLikelihood *getLikelihoodFunction() const;

        UtreeBppUtils::treemap &getTreeMap() { return dynamic_cast<RHomogeneousTreeLikelihood_PIP *>(likelihoodFunc_)->getTreemap(); };

        /**
         * @name The NNISearchable interface.
         *
         * Current implementation:
         * When testing a particular NNI, only the branch length of the parent node is optimized (and roughly).
         * All other parameters (substitution model, rate distribution and other branch length are kept at there current value.
         * When performing a NNI, only the topology change is performed.
         * This is up to the user to re-initialize the underlying likelihood data to match the new topology.
         * Usually, this is achieved by calling the topologyChangePerformed() method, which call the reInit() method of the LikelihoodData object.
         * @{
         */
        const Tree &getTopology() const { return getTree(); }


        double getTopologyValue() const throw(Exception) { return getValue(); }


        void fixTopologyChanges(tshlib::Utree *inUTree);

        void optimiseBranches(std::vector<tshlib::VirtualNode *> listNodes);


        void topologyChange(std::vector<tshlib::VirtualNode *> listNodes, tshlib::Utree *inUTree) {

            fixTopologyChanges(inUTree);
            inUTree->addVirtualRootNode();
            optimiseBranches(listNodes);
            inUTree->removeVirtualRootNode();
        }

        void topologyChangeTested(const TopologyChangeEvent &event) {
            // getLikelihoodData()->reInit();

            // if(brLenNNIParams_.size() > 0)
            fireParameterChanged(brLenTSHParams_);
            brLenTSHParams_.reset();

        }

        void topologyChangeSuccessful(const TopologyChangeEvent &event) {
            brLenTSHValues_.clear();
        }

        /** @} */
    };

} // end of namespace bpp.


#endif //MINIJATI_TSHLIBHOMOGENEOUSTREELIKELIHOOD_HPP
