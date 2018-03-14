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

        /**
         * @brief Optimizer used for testing TSH moves.
         */

        BfgsMultiDimensions *optimiser_;
        AbstractHomogeneousTreeLikelihood *likelihoodFunc_;
        const UtreeBppUtils::treemap &treemap_;
        mutable tshlib::Utree *utree_;
        mutable std::vector<DRASRTreeLikelihoodData *> referenceLikelihoodComponents_;
        mutable std::vector<Node *> likelihoodNodes_;
        mutable DRASRTreeLikelihoodData *likelihoodData_;
        VVVdouble *pxy_;


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
                                     DiscreteDistribution *rDist,
                                     tshlib::Utree *inUtree,
                                     UtreeBppUtils::treemap &inTreeMap)
        throw(Exception);

        /**
         * @brief Copy constructor.
         */
        TSHHomogeneousTreeLikelihood(const TSHHomogeneousTreeLikelihood &lik);

        TSHHomogeneousTreeLikelihood &operator=(const TSHHomogeneousTreeLikelihood &lik);

        virtual ~TSHHomogeneousTreeLikelihood();

        TSHHomogeneousTreeLikelihood *clone() const { return new TSHHomogeneousTreeLikelihood(*this); }

        void setData(const SiteContainer &sites) throw(Exception) { RHomogeneousTreeLikelihood::setData(sites); }

        AbstractHomogeneousTreeLikelihood *getLikelihoodFunction() const;

        double getLikelihoodValue() { return likelihoodFunc_->getValue(); };

        double updateLikelihood(std::vector<tshlib::VirtualNode *> &nodeList);

        const UtreeBppUtils::treemap &getTreeMap() const;

        //UtreeBppUtils::treemap &getTreeMap() const { return treemap_; };

        tshlib::Utree *getUtree() const { return utree_; }
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

        const ParameterList &getParameters() const { return likelihoodFunc_->getParameters(); }

        void fixTopologyChanges(tshlib::Utree *inUTree);

        void optimiseBranches(std::vector<tshlib::VirtualNode *> listNodes);


        void topologyChange(std::vector<tshlib::VirtualNode *> listNodes, tshlib::Utree *inUTree) {
            // Update BPP tree using the structure in Utree
            fixTopologyChanges(inUTree);
            // Add virtual root to compute the likelihood
            inUTree->addVirtualRootNode();
            // Optimise branches involved in the tree rearrangement
            optimiseBranches(listNodes);
            // Remove the virtual root to allow for further tree topology improvements
            inUTree->removeVirtualRootNode();
        }

        /*!
         * @brief This method switch the DRASRT Likelihood data arrays with the original ones
         */
        void topology() {


        }

        /*!
        * @brief This method switch the DRASRT Likelihood data arrays with the original ones
        */
        void topologyChangeSuccess() {


        }



        void topologyChangeTested(const TopologyChangeEvent &event) {
            // getLikelihoodData()->reInit();

            // if(brLenNNIParams_.size() > 0)
            //fireParameterChanged(brLenTSHParams_);
            //brLenTSHParams_.reset();

        }

        void topologyChangeSuccessful(const TopologyChangeEvent &event) {
            //brLenTSHValues_.clear();
        }


        std::vector<Node *> remapVirtualNodeLists(std::vector<tshlib::VirtualNode *> &inputList) const {

            std::vector<Node *> newList;

            for (auto &vnode:inputList) {

                newList.push_back(tree_->getNode(treemap_.right.at(vnode)));
            }

            return newList;
        }

    protected:

        void updateLikelihoodArrays(const Node *node);


        /** @} */
    };

} // end of namespace bpp.


#endif //MINIJATI_TSHLIBHOMOGENEOUSTREELIKELIHOOD_HPP
