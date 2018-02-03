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
 * @file RHomogeneousTreeLikelihood_PIP.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 11 01 2018
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
#ifndef MINIJATI_RHOMOGENEOUSTREELIKELIHOOD_PIP_HPP
#define MINIJATI_RHOMOGENEOUSTREELIKELIHOOD_PIP_HPP

#include <Bpp/Phyl/Likelihood/AbstractHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRASRTreeLikelihoodData.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

using namespace bpp;

#include "Utilities.hpp"

namespace bpp {

    class RHomogeneousTreeLikelihood_PIP : public AbstractHomogeneousTreeLikelihood {
    private:

        mutable DRASRTreeLikelihoodData *likelihoodData_;
        mutable DRASRTreeLikelihoodData *likelihoodEmptyData_;

        mutable std::vector<Node *> likelihoodNodes_;

        mutable std::map<int, std::pair<std::vector<int>, bpp::Node*>> descCountData_;
        mutable std::map<int, std::pair<std::vector<bool>, bpp::Node*>> setAData_;
        mutable std::map<const bpp::Node*, double> iotasData_;
        mutable std::map<const bpp::Node*, double> betasData_;
        mutable std::map<const bpp::Node*, std::vector<std::vector<double>>> indicatorFun_;

        mutable double nu_;
    public:

        double getNu() const { return nu_; }


    private:
        mutable double tau_;

        mutable UtreeBppUtils::treemap treemap_;
    public:
        const UtreeBppUtils::treemap &getTreemap() const { return treemap_; }

    protected:
        double minusLogLik_;

    public:
        /**
         * @brief Build a new RHomogeneousTreeLikelihood object without data.
         *
         * This constructor only initialize the parameters.
         * To compute a likelihood, you will need to call the setData() and the computeTreeLikelihood() methods.
         *
         * @param tree The tree to use.
         * @param model The substitution model to use.
         * @param rDist The rate across sites distribution to use.
         * @param checkRooted Tell if we have to check for the tree to be unrooted.
         * If true, any rooted tree will be unrooted before likelihood computation.
         * @param verbose Should I display some info?
         * @param usePatterns Tell if recursive site compression should be performed.
         * @throw Exception in an error occured.
         */
        RHomogeneousTreeLikelihood_PIP(
                const Tree &tree,
                TransitionModel *model,
                DiscreteDistribution *rDist,
                UtreeBppUtils::treemap *tm,
                bool checkRooted = true,
                bool verbose = true,
                bool usePatterns = true)
        throw(Exception);

        /**
         * @brief Build a new RHomogeneousTreeLikelihood object with data.
         *
         * This constructor initializes all parameters, data, and likelihood arrays.
         *
         * @param tree The tree to use.
         * @param data Sequences to use.
         * @param model The substitution model to use.
         * @param rDist The rate across sites distribution to use.
         * @param checkRooted Tell if we have to check for the tree to be unrooted.
         * If true, any rooted tree will be unrooted before likelihood computation.
         * @param verbose Should I display some info?
         * @param usePatterns Tell if recursive site compression should be performed.
         * @throw Exception in an error occured.
         */
        RHomogeneousTreeLikelihood_PIP(
                const Tree &tree,
                const SiteContainer &data,
                TransitionModel *model,
                DiscreteDistribution *rDist,
                UtreeBppUtils::treemap *tm,
                bool checkRooted = true,
                bool verbose = true,
                bool usePatterns = true)
        throw(Exception);

        RHomogeneousTreeLikelihood_PIP(const RHomogeneousTreeLikelihood_PIP &lik);

        RHomogeneousTreeLikelihood_PIP &operator=(const RHomogeneousTreeLikelihood_PIP &lik);

        virtual ~RHomogeneousTreeLikelihood_PIP();

        RHomogeneousTreeLikelihood_PIP *clone() const { return new RHomogeneousTreeLikelihood_PIP(*this); }


    private:

        /**
         * @brief Method called by constructors.
         */
        void init_(bool usePatterns) throw(Exception);


        void initializeLikelihoodMatrix_(VVVdouble *_likelihoods_node);

        void initializeLikelihoodEmptyMatrix_(VVVdouble *_likelihoods_empty_node);

        void hadamardMultFvSons_(Node *node) const;

        void hadamardMultFvEmptySons_(Node *node) const ;

        void SingleRateCategoryHadamardMultFvSons_(Node *node,unsigned long site,unsigned long rate,Vdouble *fv_out) const;

        void SingleRateCategoryHadamardMultFvEmptySons_(Node *node,unsigned long rate,Vdouble *fv_out) const;

        void computePrTimesFv_(Node *node) const;

        void computePrTimesFvEmpty_(Node *node) const;

        void computePrTimesIndicator_(Node *node) const;

        void computePrTimesIndicatorEmpty_(Node *node) const;

        void InitialiseInsertionHistories() const;

    public:

        /**
         * @name The TreeLikelihood interface.
         *
         * Other methods are implemented in the AbstractHomogeneousTreeLikelihood class.
         *
         * @{
         */
        void setData(const SiteContainer &sites) throw(Exception);

        double getLikelihood() const;

        double getLogLikelihood() const;
        double getLogLikelihood(std::vector<tshlib::VirtualNode *> &listNodes) const;


        double getLikelihoodForASite(size_t site) const;

        double getLogLikelihoodForASite(size_t site) const;



        double getLogLikelihoodSubtree(const Node *node) const;
        double getLogLikelihoodSubtreeForASite(size_t site) const;
        double getLogLikelihoodSubtreeForASiteForARateClass(size_t site, size_t rateClass) const;


        /** @} */


        /**
         * @name The DiscreteRatesAcrossSites interface implementation:
         * @deprecated at the moment these methods are not used
         * @{
         */
        double getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

        double getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

        double getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const;

        double getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const;
        /** @} */

        /**
         * @brief Implements the Function interface.
         *
         * Update the parameter list and call the applyParameters() method.
         * Then compute the likelihoods at each node (computeLikelihood() method)
         * and call the getLogLikelihood() method.
         *
         * If a subset of the whole parameter list is passed to the function,
         * only these parameters are updated and the other remain constant (i.e.
         * equal to their last value).
         *
         * @param parameters The parameter list to pass to the function.
         */
        void setParameters(const ParameterList &parameters) throw(ParameterNotFoundException, ConstraintException);

        double getValue() const throw(Exception);

        size_t getSiteIndex(size_t site) const throw(IndexOutOfBoundsException) { return likelihoodData_->getRootArrayPosition(site); }

        /**
         * @name DerivableFirstOrder interface.
         *
         * @{
         */
        double getFirstOrderDerivative(const std::string &variable) const throw(Exception) { return 0; }
        /** @} */

        /**
         * @name DerivableSecondOrder interface.
         *
         * @{
         */
        double getSecondOrderDerivative(const std::string &variable) const throw(Exception) { return 0; }

        double getSecondOrderDerivative(const std::string &variable1, const std::string &variable2) const throw(Exception) { return 0; } // Not implemented for now.
        /** @} */

    public:    // Specific methods:

        std::vector<int> getNodeDescCounts(bpp::Node *node, int siteId){ return descCountData_[node->getId()].first;}

        int getNodeDescCountForASite(bpp::Node *node, int siteId) const{ return descCountData_[node->getId()].first.at(siteId);}

        bool getSetAForANodeForASite(bpp::Node *node, int siteId){ return setAData_[node->getId()].first.at(siteId);}

        DRASRTreeLikelihoodData *getLikelihoodData() { return likelihoodData_; }

        const DRASRTreeLikelihoodData *getLikelihoodData() const { return likelihoodData_; }

        /**
         * @name Interface to compute the likelihood components
         * @{
         */

        /**
         * @brief This method computes the likelihood of the tree  generating automatically postorder-traversal node list
         */
        void computeTreeLikelihood();

        /**
         * @brief This method computes the likelihood of the tree for a list of nodes computed using a postorder-traversal
         * @param nodeList
         */
        void computeTreeLikelihood(std::vector<Node *> nodeList);
        /** @} */

        //virtual double getDLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;
        //virtual double getDLikelihoodForASite(size_t site) const;
        //virtual double getDLogLikelihoodForASite(size_t site) const;
        //virtual double getDLogLikelihood() const;
        //virtual void computeTreeDLikelihood(const std::string &variable);
        //virtual double getD2LikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;
        //virtual double getD2LikelihoodForASite(size_t site) const;
        //virtual double getD2LogLikelihoodForASite(size_t site) const;
        //virtual double getD2LogLikelihood() const;
        //virtual void computeTreeD2Likelihood(const std::string &variable);

        /**
         * @brief This method computes the likelihood after a tree rearrangment
         * @return The likelihood value using the intermediate partial values
         */
        double computeLikelihoodOnTreeRearrangment(std::vector<tshlib::VirtualNode *> &listNodes) const;

        /**
         * @brief This method computes a list of nodes traversing the tree in postorder
         *
         */
        void computePostOrderNodeList(Node *startNode) const;

    protected:

        /**
         * @brief Compute the likelihood for a subtree defined by the Tree::Node <i>node</i>.
         *        This function should be called only for filling the likelihood arrays (i.e. first traversal, parameter change -- but not topology changes).
         *
         * @param node The root of the subtree.
         */
        virtual void computeSubtreeLikelihood();

        //virtual void computeDownSubtreeDLikelihood(const Node *);
        //virtual void computeDownSubtreeD2Likelihood(const Node *);

        void fireParameterChanged(const ParameterList &params);

        /**
         * @brief This method is mainly for debugging purpose.
         *
         * @param node The node at which likelihood values must be displayed.
         */
        virtual void displayLikelihood(const Node *node);

        /**
         * @brief This method sets DescCount (number of characters different from gap per column) value for all the nodes in the tree
         */
        virtual void setInsertionHistories(const SiteContainer &sites) const;

        /**
         * @brief This method sets the setA (setA=1: possible insertion on that edge) value for all the nodes in the tree
         */
        /*virtual void setAllSetAData(const SiteContainer &sites) const;*/

        /**
         * @brief This method sets the iota value for all the nodes in the tree
         */
        virtual void setAllIotas();

        /**
         * @brief This method sets the beta value for all the nodes in the tree
         */
        virtual void setAllBetas();

        /**
         * @brief This method sets to 1 all the likelihood arrays recursively from a starting node
         * @param node The node at which the likelihood arrays must be reset
         */
        virtual void resetNodeLikelihoodArrays(const Node *node);

        /**
         * @brief This method updates the likelihood arrays recursively from a starting node for a
         *        subtree
         * @param nodelist The postorder list of nodes at which the likelihood arrays must be updated
         */
        virtual void recombineFvAfterMove() const;

        virtual void recombineFvAtNode(Node *node) const;

        void setIndicatorFunction(const SiteContainer &sites) const;

        double computePhi(double lkEmptyColumn) const;

        void computeNu();

        void printFV(Node *node, VVVdouble *likelihoodvector);

        void printPrMatrix(Node *node, VVdouble *pr);

        double computeLikelihoodWholeAlignment()const;

        std::vector<Node *> remapVirtualNodeLists(std::vector<tshlib::VirtualNode *> &inputList) const;

        void ExtendNodeListOnSetA(tshlib::VirtualNode *qnode, std::vector<Node *> &listNodes, unsigned long site) const;

        //friend class RHomogeneousMixedTreeLikelihood;
        double computeLikelihoodSite(size_t i) const;

        double computeLikelihoodWholeAlignmentEmptyColumn() const;

        double computeLikelihoodWholeSites() const;

        int countNonGapCharacterInSite(const SiteContainer &sites, int siteID) const;
    };
}
#endif //MINIJATI_RHOMOGENEOUSTREELIKELIHOOD_PIP_HPP
