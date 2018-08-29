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
 * @file UnifiedTSHomogeneousTreeLikelihood_PIP.cpp
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

#include <glog/logging.h>
#include "UnifiedTSHomogeneousTreeLikelihood_PIP.hpp"


using namespace bpp;

/*
 * Implementation for the interface  likelihood under tree search engines (all mixed models)
 */

UnifiedTSHomogeneousTreeLikelihood_PIP::UnifiedTSHomogeneousTreeLikelihood_PIP(const Tree &tree,
                                                                               const SiteContainer &data,
                                                                               TransitionModel *model,
                                                                               DiscreteDistribution *rDist,
                                                                               tshlib::Utree *utree,
                                                                               UtreeBppUtils::treemap *tm,
                                                                               bool optNumericalDerivatives,
                                                                               std::map<std::string, std::string> &params,
                                                                               const std::string &suffix,
                                                                               bool checkRooted,
                                                                               bool verbose,
                                                                               bool usePatterns) :
        RHomogeneousTreeLikelihood_PIP(tree, data, model, rDist, tm, checkRooted, verbose, usePatterns), utree_(utree) {


    setOptimiser(static_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(this), optNumericalDerivatives, params, suffix, true, verbose, 0);

}

UnifiedTSHomogeneousTreeLikelihood_PIP::UnifiedTSHomogeneousTreeLikelihood_PIP(const Tree &tree,
                                                                               TransitionModel *model,
                                                                               DiscreteDistribution *rDist,
                                                                               tshlib::Utree *utree,
                                                                               UtreeBppUtils::treemap *tm,
                                                                               bool optNumericalDerivatives,
                                                                               std::map<std::string, std::string> &params,
                                                                               const std::string &suffix,
                                                                               bool checkRooted,
                                                                               bool verbose,
                                                                               bool usePatterns) :
        RHomogeneousTreeLikelihood_PIP(tree, model, rDist, tm, checkRooted, verbose, usePatterns), utree_(utree) {


    setOptimiser(static_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(this), optNumericalDerivatives, params, suffix, true, verbose, 0);


}


UnifiedTSHomogeneousTreeLikelihood_PIP::~UnifiedTSHomogeneousTreeLikelihood_PIP() {

}


void UnifiedTSHomogeneousTreeLikelihood_PIP::init_(bool usePatterns) {

}

void UnifiedTSHomogeneousTreeLikelihood_PIP::addTestLikelihoodData(int idxThread) {

    // allocate thread memory
    for (auto &nodeID:tree_->getNodesId()) {

        // Instantiation empty objects
        std::vector<std::vector<std::vector<double>>> sub;
        std::vector<std::vector<std::vector<double>>> empty;
        std::vector<int> descCount;
        std::vector<bool> setA;

        // Allocate memory
        descCount.resize(nbDistinctSites_);
        setA.resize(nbDistinctSites_);

        VectorTools::resize3(sub, nbDistinctSites_, nbClasses_, nbStates_);
        VectorTools::resize3(empty, 1, nbClasses_, nbStates_);

        // Store references
        testVectorLikelihoodData_[LKDataClass::sub][idxThread][nodeID] = sub;
        testVectorLikelihoodData_[LKDataClass::empty][idxThread][nodeID] = empty;
        tsTemp_descCountData_[idxThread][nodeID] = descCount;
        tsTemp_setAData_[idxThread][nodeID] = setA;
        tsTemp_node_data_origin[idxThread][nodeID] = false;

    }

}

void UnifiedTSHomogeneousTreeLikelihood_PIP::removeTestLikelihoodData(int idxThread) {

    // deallocate thread memory
    for (auto &nodeID:tree_->getNodesId()) {

        size_t dim1 = testVectorLikelihoodData_[LKDataClass::sub][idxThread][nodeID].size();
        size_t dim2 = testVectorLikelihoodData_[LKDataClass::sub][idxThread][nodeID][0].size();

        for (size_t i = 0; i < dim1; ++i) {
            for (size_t j = 0; j < dim2; ++j)
                testVectorLikelihoodData_[LKDataClass::sub][idxThread][nodeID][i][j].resize(0);

            testVectorLikelihoodData_[LKDataClass::sub][idxThread][nodeID][i].resize(0);
        }
        testVectorLikelihoodData_[LKDataClass::empty][idxThread][nodeID][0].resize(0);

        tsTemp_descCountData_[idxThread][nodeID].resize(0);
        tsTemp_setAData_[idxThread][nodeID].resize(0);

    }

}

void UnifiedTSHomogeneousTreeLikelihood_PIP::fireTopologyChange(std::vector<int> nodeList) {

    // Recompute the value of the FV 3D arrays
    //computeSubtreeLikelihood(likelihoodDataTest_, likelihoodEmptyDataTest_);
    computeSubtreeLikelihood(nodeList);
    // Compute the insertion histories set (recompute the desc_count and set A)
    setInsertionHistories(*data_, nodeList, &descCountData_, &setAData_);
}


void UnifiedTSHomogeneousTreeLikelihood_PIP::fireTopologyChange(std::vector<int> nodeList, std::map<int, VVVdouble> *ts_lkdata,
                                                                std::map<int, VVVdouble> *ts_lkemptydata,
                                                                std::map<int, std::vector<int>> *ts_desccount,
                                                                std::map<int, std::vector<bool>> *ts_setadata,
                                                                std::map<int, bool> *ts_node__data_origin,
                                                                tshlib::Utree &_utree__topology) {

    // Recompute the value of the FV 3D arrays
    //computeSubtreeLikelihood(likelihoodDataTest_, likelihoodEmptyDataTest_);
    computeSubtreeLikelihood(ts_lkdata, ts_lkemptydata, nodeList, ts_node__data_origin, _utree__topology);
    // Compute the insertion histories set (recompute the desc_count and set A)
    setInsertionHistories(*data_, nodeList, ts_desccount, ts_setadata, _utree__topology);
}


double UnifiedTSHomogeneousTreeLikelihood_PIP::updateLikelihoodOnTreeRearrangement(std::vector<int> &nodeList,
                                                                                   tshlib::Utree &_utree__topology,
                                                                                   int idxThread) {
    //fetch temporary arrays
    std::map<int, VVVdouble> *ts_lkdata = &testVectorLikelihoodData_[LKDataClass::sub][idxThread];
    std::map<int, VVVdouble> *ts_lkemptydata = &testVectorLikelihoodData_[LKDataClass::empty][idxThread];
    std::map<int, std::vector<int>> *ts_desccount = &tsTemp_descCountData_[idxThread];
    std::map<int, std::vector<bool>> *ts_setadata = &tsTemp_setAData_[idxThread];
    std::map<int, bool> *ts_node__data_origin = &tsTemp_node_data_origin[idxThread];


    // Add root to the utree structure
    _utree__topology.addVirtualRootNode();

    // 0. convert the list of tshlib::VirtualNodes into bpp::Node
    std::vector<int> _affected__nodes = remapVirtualNodeLists(nodeList, _utree__topology);

    // 1. Flag nodes to read from reference
    for (std::vector<int>::iterator it = _affected__nodes.begin(); it != _affected__nodes.end(); ++it)
        (*ts_node__data_origin)[(*it)] = true;

    // 2. Fire topology change
    fireTopologyChange(_affected__nodes, ts_lkdata, ts_lkemptydata, ts_desccount, ts_setadata, ts_node__data_origin, _utree__topology);

    // 3. Compute loglikelihood
    double logLk = getLogLikelihoodOnTreeRearrangement(_affected__nodes,
                                                       ts_lkdata,
                                                       ts_lkemptydata,
                                                       ts_desccount,
                                                       ts_setadata,
                                                       ts_node__data_origin,
                                                       _utree__topology);

    // 4. Remove root node from the utree structure
    _utree__topology.removeVirtualRootNode();

    // 5. Remove flags
    for (std::vector<int>::iterator it = _affected__nodes.begin(); it != _affected__nodes.end(); ++it)
        (*ts_node__data_origin)[(*it)] = false;


    return logLk;

}

#ifdef INTELTBB

void UnifiedTSHomogeneousTreeLikelihood_PIP::recomputeSiteLikelihoodUsingPartitions(const tbb::blocked_range<size_t>& range, std::vector<double> *lk_sites) const{
    for (size_t i = range.begin(); i < range.end(); ++i) {

        std::vector<int> tempExtendedNodeList;
        const std::vector<unsigned int> *rootWeights = &likelihoodData_->getWeights();

        // Extend it
        _extendNodeListOnSetA(likelihoodNodes_.back(), tempExtendedNodeList, i);

        // call to function which retrieves the lk value for each site
        (*lk_sites)[i] = log(computeLikelihoodForASite(tempExtendedNodeList, i)) * rootWeights->at(i);

        // Debug
        DVLOG(2) << "site log_lk[" << i << "]=" << std::setprecision(18) << (*lk_sites)[i] << std::endl;
    }
};

#endif

double UnifiedTSHomogeneousTreeLikelihood_PIP::getLogLikelihoodOnTreeRearrangement(const std::vector<int> &nodeList,
                                                                                   std::map<int, VVVdouble> *ts_lkdata,
                                                                                   std::map<int, VVVdouble> *ts_lkemptydata,
                                                                                   std::map<int, std::vector<int>> *ts_desccount,
                                                                                   std::map<int, std::vector<bool>> *ts_setadata,
                                                                                   std::map<int, bool> *ts_node__data_origin,
                                                                                   tshlib::Utree &_utree__topology) const {

    // 1. Initialise variables and contenitors
    double logLK;
    std::vector<double> lk_sites(nbDistinctSites_);

    // 2. Compute the lk of the empty column
    double lk_site_empty = computeLikelihoodWholeAlignmentEmptyColumn(ts_lkemptydata, _utree__topology);

    // 3. Compute the likelihood of each site
    const std::vector<unsigned int> *rootWeights = &likelihoodData_->getWeights();

    for (unsigned long i = 0; i < nbDistinctSites_; i++) {

        // Extend rearranged-node-list including all the nodes in the setA for each site
        std::vector<int> tempExtendedNodeList;
        _extendNodeListOnSetA(nodeList.back(), tempExtendedNodeList, i, _utree__topology);

        // call to function which retrieves the lk value for each site
        lk_sites[i] = log(computeLikelihoodForASite(tempExtendedNodeList, i, _utree__topology)) * rootWeights->at(i);


        DVLOG(2) << "site log_lk[" << i << "]=" << std::setprecision(18) << lk_sites[i] << std::endl;
    }

    // Sum all the values stored in the lk vector
    logLK = MatrixBppUtils::sumVector(&lk_sites);
    DVLOG(2) << "LK Sites [BPP] " << std::setprecision(18) << logLK;

    // compute PHi
    double log_phi_value = computePhi(lk_site_empty);
    DVLOG(2) << "PHI [BPP] " << std::setprecision(18) << log_phi_value;

    logLK += log_phi_value;


    return logLK;
}


void UnifiedTSHomogeneousTreeLikelihood_PIP::topologyChangeSuccessful(std::vector<int> listNodes) {

    // Update BPP tree using the structure in Utree
    topologyCommitTree();

    // Add virtual root to compute the likelihood
    utree_->addVirtualRootNode();

    // Fire topology change
    std::vector<int> ponl = getNodeListPostOrder(tree_->getRootNode()->getId());

    setLikelihoodNodes(ponl);

    fireTopologyChange(ponl);

    // Optimise branches involved in the tree rearrangement
    fireBranchOptimisation(UtreeBppUtils::remapNodeLists(listNodes, tree_, treemap_));

    // Remove the virtual root to allow for further tree topology improvements
    utree_->removeVirtualRootNode();

}


void UnifiedTSHomogeneousTreeLikelihood_PIP::topologyCommitTree() {

    std::vector<tshlib::VirtualNode *> nodelist;
    nodelist = utree_->listVNodes;

    std::map<int, bpp::Node *> tempMap;
    std::map<int, double> tempDistanceToFather;
    // reset inBtree
    for (auto &bnode:tree_->getNodes()) {

        tempMap.insert(std::pair<int, bpp::Node *>(bnode->getId(), bnode));
        // Empty array of sons on the node
        bnode->removeSons();

        if (bnode->hasFather()) {

            tempDistanceToFather.insert(std::pair<int, double>(bnode->getId(), bnode->getDistanceToFather()));
            // Empty father connection
            bnode->removeFather();
        }
    }

    for (auto &vnode:nodelist) {

        if (!vnode->isTerminalNode()) {

            // get corrisponding sons in inBTree
            bpp::Node *leftBNode = tempMap[treemap_.right.at(vnode->getNodeLeft()->getVnode_id())];
            bpp::Node *rightBNode = tempMap[treemap_.right.at(vnode->getNodeRight()->getVnode_id())];

            // get corrisponding parent in inBTree
            bpp::Node *pNode = tempMap[treemap_.right.at(vnode->getVnode_id())];

            leftBNode->setFather(pNode);
            rightBNode->setFather(pNode);

            leftBNode->setDistanceToFather(tempDistanceToFather[leftBNode->getId()]);
            rightBNode->setDistanceToFather(tempDistanceToFather[rightBNode->getId()]);

            //Add new sons
            pNode->setSon(0, leftBNode);
            pNode->setSon(1, rightBNode);

            pNode->setDistanceToFather(tempDistanceToFather[pNode->getId()]);

        }

        // in case the current vnode is also the pseudo-root
        if (vnode == vnode->getNodeUp()->getNodeUp()) {

            bpp::Node *leftBNode = tempMap[treemap_.right.at(vnode->getVnode_id())];
            bpp::Node *rightBNode = tempMap[treemap_.right.at(vnode->getNodeUp()->getVnode_id())];

            tree_->getRootNode()->removeSons();

            leftBNode->setFather(tree_->getRootNode());
            rightBNode->setFather(tree_->getRootNode());

            leftBNode->setDistanceToFather(tempDistanceToFather[leftBNode->getId()]);
            rightBNode->setDistanceToFather(tempDistanceToFather[rightBNode->getId()]);

            tree_->getRootNode()->setSon(0, leftBNode);
            tree_->getRootNode()->setSon(1, rightBNode);

        }

    }

}
