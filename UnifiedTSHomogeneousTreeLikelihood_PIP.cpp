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



    setOptimiser(static_cast<UnifiedTSHomogeneousTreeLikelihood_PIP*>(this), optNumericalDerivatives, params, suffix, true, verbose, 0);



}


UnifiedTSHomogeneousTreeLikelihood_PIP::~UnifiedTSHomogeneousTreeLikelihood_PIP() {

}


void UnifiedTSHomogeneousTreeLikelihood_PIP::init_(bool usePatterns) {

    //likelihoodData_ = new DRASRTreeLikelihoodData(tree_, rateDistribution_->getNumberOfCategories(), usePatterns);        // FV
    //likelihoodEmptyData_ = new DRASRTreeLikelihoodData(tree_, rateDistribution_->getNumberOfCategories(), usePatterns);         // FV empty

    likelihoodDataTest_ = new DRASRTreeLikelihoodData(tree_, rateDistribution_->getNumberOfCategories(), usePatterns);    // FV test
    likelihoodEmptyDataTest_ = new DRASRTreeLikelihoodData(tree_, rateDistribution_->getNumberOfCategories(), usePatterns);    // FV empty test

}


void UnifiedTSHomogeneousTreeLikelihood_PIP::fireTopologyChange(std::vector<int> nodeList) {
    // Store the nodes where the likelihood should be recomputed in post-order
    setLikelihoodNodes(nodeList);
    // Recompute the value of the FV 3D arrays
    computeSubtreeLikelihood();
    // Compute the insertion histories set (recompute the desc_count and set A)
    setInsertionHistories(*data_);
}


double UnifiedTSHomogeneousTreeLikelihood_PIP::updateLikelihoodOnTreeRearrangement(std::vector<tshlib::VirtualNode *> &nodeList) {

    // Add root to the utree structure
    utree_->addVirtualRootNode();

    // 0. convert the list of tshlib::VirtualNodes into bpp::Node
    std::vector<int> rearrangedNodes = remapVirtualNodeLists(nodeList);

    // 1. Fire topology change
    fireTopologyChange(rearrangedNodes);

    // 2. Compute loglikelihood
    //double logLk = getLogLikelihoodOnTopologyChange();
    double logLk = getLogLikelihoodOnTreeRearrangement();

    // Remove root node from the utree structure
    utree_->removeVirtualRootNode();

    return logLk;

}


double UnifiedTSHomogeneousTreeLikelihood_PIP::getLogLikelihoodOnTreeRearrangement() const {
    double logLK;

    // 2. Compute the lk of the empty column
    double lk_site_empty = computeLikelihoodWholeAlignmentEmptyColumn();

    // 3. Compute the likelihood of each site
    std::vector<double> lk_sites(nbDistinctSites_);

    const std::vector<unsigned int> *rootWeights = &likelihoodData_->getWeights();

    for (unsigned long i = 0; i < nbDistinctSites_; i++) {

        std::vector<int> tempExtendedNodeList;

        // Extend it
        _extendNodeListOnSetA(likelihoodNodes_.back(), tempExtendedNodeList, i);

        // call to function which retrieves the lk value for each site
        lk_sites[i] = log(computeLikelihoodForASite(tempExtendedNodeList, i)) * rootWeights->at(i);
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


void UnifiedTSHomogeneousTreeLikelihood_PIP::topologyChangeSuccessful(std::vector<tshlib::VirtualNode *> listNodes) {

    // Update BPP tree using the structure in Utree
    topologyCommitTree();

    // Add virtual root to compute the likelihood
    utree_->addVirtualRootNode();

    // remap the virtual nodes to the bpp nodes
    std::vector<Node *> extractionNodes = UtreeBppUtils::remapNodeLists(listNodes, tree_, treemap_);

    // Optimise branches involved in the tree rearrangement
    fireBranchOptimisation(extractionNodes);

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

        if (bnode->hasFather()){

            tempDistanceToFather.insert(std::pair<int, double>(bnode->getId(), bnode->getDistanceToFather()));
            // Empty father connection
            bnode->removeFather();
        }
    }

    for (auto &vnode:nodelist) {

        //std::cerr << "vnode " << vnode->getNodeName();
        if (!vnode->isTerminalNode()) {

            // get corrisponding sons in inBTree
            bpp::Node *leftBNode = tempMap[treemap_.right.at(vnode->getNodeLeft())];
            bpp::Node *rightBNode = tempMap[treemap_.right.at(vnode->getNodeRight())];

            // get corrisponding parent in inBTree
            bpp::Node *pNode = tempMap[treemap_.right.at(vnode)];

            leftBNode->setFather(pNode);
            rightBNode->setFather(pNode);

            //leftBNode->setDistanceToFather(tree_->getDistanceToFather(leftBNode->getId()));
            //rightBNode->setDistanceToFather(tree_->getDistanceToFather(rightBNode->getId()));

            leftBNode->setDistanceToFather(tempDistanceToFather[leftBNode->getId()]);
            rightBNode->setDistanceToFather(tempDistanceToFather[rightBNode->getId()]);

            //Add new sons
            pNode->setSon(0, leftBNode);
            pNode->setSon(1, rightBNode);
            //pNode->setDistanceToFather(tree_->getDistanceToFather(pNode->getId()));
            pNode->setDistanceToFather(tempDistanceToFather[pNode->getId()]);
            //std::cerr << "\t internal";

        } else {

            //std::cerr << "\t leaf";

        }
        // in case the current vnode is also the pseudo-root
        if (vnode == vnode->getNodeUp()->getNodeUp()) {
            //std::cerr << "\tvnode pseudoroot";

            bpp::Node *leftBNode = tempMap[treemap_.right.at(vnode)];
            bpp::Node *rightBNode = tempMap[treemap_.right.at(vnode->getNodeUp())];

            tree_->getRootNode()->removeSons();

            leftBNode->setFather(tree_->getRootNode());
            rightBNode->setFather(tree_->getRootNode());

            //leftBNode->setDistanceToFather(tree_->getDistanceToFather(leftBNode->getId()));
            //rightBNode->setDistanceToFather(tree_->getDistanceToFather(rightBNode->getId()));

            leftBNode->setDistanceToFather(tempDistanceToFather[leftBNode->getId()]);
            rightBNode->setDistanceToFather(tempDistanceToFather[rightBNode->getId()]);

            tree_->getRootNode()->setSon(0, leftBNode);
            tree_->getRootNode()->setSon(1, rightBNode);

        }

        //std::cerr << "\t done\n";

    }

}
