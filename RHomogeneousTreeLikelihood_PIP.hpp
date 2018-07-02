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
#ifndef MINIJATI_RHOMOGENEOUSTREELIKELIHOOD_PIP_HPP
#define MINIJATI_RHOMOGENEOUSTREELIKELIHOOD_PIP_HPP

#include <Bpp/Phyl/Likelihood/AbstractHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRASRTreeLikelihoodData.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>

using namespace bpp;

#include "Utilities.hpp"

namespace bpp {

    class RHomogeneousTreeLikelihood_PIP : public AbstractHomogeneousTreeLikelihood {
    protected:

        mutable DRASRTreeLikelihoodData *likelihoodData_;
        mutable DRASRTreeLikelihoodData *likelihoodEmptyData_;

        mutable std::vector<int> likelihoodNodes_;                            //The node is represented via its <int> ID

        mutable std::map<int, std::pair<std::vector<int>, int>> descCountData_;
        mutable std::map<int, std::pair<std::vector<bool>, int>> setAData_;                   // SetA flags if a node should be included in the insertion histories

        mutable std::map<int, std::pair<std::vector<int>, bool>> descGapCountData_;           // Counts the number of gaps in the descending nodes
        mutable std::map<int, int> evolutionaryEvents_;                                       // Evolutionary events counts the number of possible evolutionary events happened on the node
        mutable std::map<int, double> evolutionaryEventsWeighted_;                            // Evolutionary events counts the number of possible evolutionary events happened on the node
                                                                                              // weighted them by the number of nodes in the clade
        mutable std::map<int, double> iotasData_;
        mutable std::map<int, double> betasData_;
        mutable std::map<int, std::vector<std::vector<double>>> indicatorFun_;
        mutable std::vector<unsigned long> rootPatternLinksInverse_;

        mutable double nu_;
        mutable double tau_;

        mutable UtreeBppUtils::treemap treemap_;


    protected:

        double minusLogLik_;

        double d1bl_;
        double d2bl_;

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


    protected:

        /**
         * @brief Method called by constructors.
         */
        void init_(bool usePatterns) throw(Exception);

        void _hadamardMultFvSons(Node *node) const;

        void _hadamardMultFvEmptySons(Node *node) const;

        void _SingleRateCategoryHadamardMultFvSons(Node *node, unsigned long site, unsigned long rate, Vdouble *fv_out) const;

        void _SingleRateCategoryHadamardMultFvEmptySons(Node *node, unsigned long rate, Vdouble *fv_out) const;

        void _computePrTimesFv(Node *node) const;

        void _computePrTimesFvEmpty(Node *node) const;

        void _computePrTimesIndicator(Node *node) const;

        void _computePrTimesIndicatorEmpty(Node *node) const;

        void initialiseInsertionHistories() const;

    public:

        /**
         * @name The TreeLikelihood interface.
         *
         * Other methods are implemented in the AbstractHomogeneousTreeLikelihood class.
         *
         * @{
         */
        void setData(const SiteContainer &sites) throw(Exception);

        double getNu() const { return nu_; }

        UtreeBppUtils::treemap &getTreemap() { return treemap_; }

        double getLikelihood() const {
            std::cerr << "getLikelihood()" << std::endl;
            return 0;
        };

        double getLogLikelihood() const;


        //double getLogLikelihood(std::vector<tshlib::VirtualNode *> &listNodes) const;
        //double getLogLikelihoodSubtree(const Node *node) const;
        //double getLogLikelihoodSubtreeForASite(size_t site) const;
        //double getLogLikelihoodSubtreeForASiteForARateClass(size_t site, size_t rateClass) const;


        /** @} */


        /**
         * @name The DiscreteRatesAcrossSites interface implementation:
         * @deprecated at the moment these methods are not used
         * @{
         */
        double getLikelihoodForASite(size_t site) const {
            std::cerr << "getLikelihoodForASite()" << std::endl;
            return 0;
        };

        double getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const {
            std::cerr << "getLikelihoodForASiteForARateClass()" << std::endl;
            return 0;
        };

        double getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const {
            std::cerr << "getLikelihoodForASiteForARateClassForAState()" << std::endl;
            return 0;
        };

        double getLogLikelihoodForASite(size_t site) const {
            std::cerr << "getLogLikelihoodForASite()" << std::endl;
            return 0;
        };

        double getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const {
            std::cerr << "getLogLikelihoodForASiteForARateClass()" << std::endl;
            return 0;
        };

        double getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const {
            std::cerr << "getLogLikelihoodForASiteForARateClassForAState()" << std::endl;
            return 0;
        };

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
        double getFirstOrderDerivative(const std::string &variable) const throw(Exception);

        double computeN1DerivativeLikelihood(const std::string &variable);

        double computeN2DerivativeLikelihood(const std::string &variable);

        double evaluateLikelihoodPointForBranchDerivative(const std::string &variable, double new_branchlength);
        /** @} */

        /**
         * @name DerivableSecondOrder interface.
         *
         * @{
         */
        double getSecondOrderDerivative(const std::string &variable) const throw(Exception);

        double getSecondOrderDerivative(const std::string &variable1, const std::string &variable2) const throw(Exception) { return 0; } // Not implemented for now.
        /** @} */

    public:    // Specific methods:


        DRASRTreeLikelihoodData *getLikelihoodData() { return likelihoodData_; }

        const DRASRTreeLikelihoodData *getLikelihoodData() const { return likelihoodData_; }

        std::vector<int> getNodeDescCounts(bpp::Node *node, int siteId) { return descCountData_[node->getId()].first; }

        int getNodeDescCountForASite(bpp::Node *node, int siteId) const { return descCountData_[node->getId()].first.at(siteId); }

        bool getSetAForANodeForASite(bpp::Node *node, int siteId) { return setAData_[node->getId()].first.at(siteId); }


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
        void computeTreeLikelihood(std::vector<int> nodeList);

        /** @} */

        virtual double getDLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

        virtual double getDLikelihoodForASite(size_t site) const;

        virtual double getDLogLikelihoodForASite(size_t site) const;

        virtual double getDLogLikelihood() const;

        virtual void computeTreeDLikelihood(const std::string &variable);

        virtual double getD2LikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

        virtual double getD2LikelihoodForASite(size_t site) const;

        virtual double getD2LogLikelihoodForASite(size_t site) const;

        virtual double getD2LogLikelihood() const;

        virtual void computeTreeD2Likelihood(const std::string &variable);

        /**
         * @brief This method computes the likelihood after a tree rearrangment
         * @return The likelihood value using the intermediate partial values
         */
        void fireTopologyChange(std::vector<int> nodeList);

        double getLogLikelihoodOnTopologyChange() const;

        /**
         * @brief This method computes a list of nodes traversing the tree in postorder
         *
         */
        std::vector<int> getNodeListPostOrder(int startNodeID) const;

        void getNodeListPostOrder_(std::vector<int> &nodeList, int startNodeID) const;

        void setLikelihoodNodes(std::vector<int> &nodeList) const;

    protected:

        /**
         * @brief Compute the likelihood for a subtree defined by the Tree::Node <i>node</i>.
         *        This function should be called only for filling the likelihood arrays (i.e. first traversal, parameter change -- but not topology changes).
         *
         * @param node The root of the subtree.
         */
        virtual void computeSubtreeLikelihood();

        virtual void computeSubtreeLikelihood() const;

        virtual void computeDownSubtreeDLikelihood(const Node *);

        virtual void computeDownSubtreeD2Likelihood(const Node *);

        void fireParameterChanged(const ParameterList &params);

        /**
         * @brief This method is mainly for debugging purpose.
         *
         * @param node The node at which likelihood values must be displayed.
         */
        virtual void displayLikelihood(const Node *node);

        /**
         * @brief This method sets DescCount (number of characters different from gap per column) value for all the nodes in the tree
         * @param sites SiteContainer of the aligned sites
         */
        virtual void setInsertionHistories(const SiteContainer &sites) const;

        /**
         * @brief This method sets the indicator for the number of evolutionary events (Insertions and Deletions) at each node of the topology
         * @param sites SiteContainer reference of the aligned sites
         */

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
        //virtual void resetNodeLikelihoodArrays(const Node *node);

        /**
         * @brief This method updates the likelihood arrays recursively from a starting node for a
         *        subtree
         * @param nodelist The postorder list of nodes at which the likelihood arrays must be updated
         */
        //virtual void recombineFvAfterMove() const;

        //virtual void recombineFvAtNode(Node *node) const;

        void setIndicatorFunction(const SiteContainer &sites) const;

        double computePhi(double lkEmptyColumn) const;

        void computeNu();

        void _printFV(Node *node, VVVdouble *likelihoodvector) const;

        void _printPrMatrix(Node *node, VVdouble *pr);

        std::vector<int> remapVirtualNodeLists(std::vector<tshlib::VirtualNode *> &inputList) const;

        void _extendNodeListOnSetA(tshlib::VirtualNode *qnode, std::vector<int> &listNodes, unsigned long site) const;

        void _extendNodeListOnSetA(int qnodeID, std::vector<int> &listNodes, unsigned long site) const;

        double computeLikelihoodForASite(std::vector<int> &likelihoodNodes, size_t i) const;

        double computeLikelihoodWholeAlignmentEmptyColumn() const;

        int countNonGapCharacterInSite(const SiteContainer &sites, int siteID) const;

        void computeInDelDispersionOnTree(const SiteContainer &sites);

        SiteContainer *getSubAlignment(const SiteContainer &sites, std::vector<std::string> subsetSequences);

        void setNhNgOnNode(SiteContainer &sites, int nodeID);

        double getNodeAge(int nodeID);


        //friend class RHomogeneousMixedTreeLikelihood;

    };
}
#endif //MINIJATI_RHOMOGENEOUSTREELIKELIHOOD_PIP_HPP
