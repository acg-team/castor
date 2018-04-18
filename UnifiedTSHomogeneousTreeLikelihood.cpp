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
 * @file UnifiedTSHomogeneousTreeLikelihood.cpp
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
#include <glog/logging.h>
#include "UnifiedTSHomogeneousTreeLikelihood.hpp"
using namespace bpp;

/*
 * Implementation for the interface  likelihood under tree search engines (all canonical models)
 */

UnifiedTSHomogeneousTreeLikelihood::UnifiedTSHomogeneousTreeLikelihood(const Tree &tree,
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
                                                                       bool usePatterns) :
        RHomogeneousTreeLikelihood(tree, data, model, rDist, checkRooted, verbose, usePatterns) {

}

UnifiedTSHomogeneousTreeLikelihood::UnifiedTSHomogeneousTreeLikelihood(const Tree &tree,
                                                                       TransitionModel *model,
                                                                       DiscreteDistribution *rDist,
                                                                       tshlib::Utree *utree_,
                                                                       UtreeBppUtils::treemap *treemap_,
                                                                       bool optNumericalDerivatives,
                                                                       std::map<std::string, std::string> &params,
                                                                       const std::string &suffix,
                                                                       bool checkRooted,
                                                                       bool verbose,
                                                                       bool usePatterns) :
        RHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose, usePatterns){

}


UnifiedTSHomogeneousTreeLikelihood::~UnifiedTSHomogeneousTreeLikelihood() {}


void UnifiedTSHomogeneousTreeLikelihood::init_(bool usePatterns) {

    likelihoodDataTest_ = new DRASRTreeLikelihoodData(tree_, rateDistribution_->getNumberOfCategories(), usePatterns);    // FV

}


void UnifiedTSHomogeneousTreeLikelihood::fireTopologyChange(std::vector<int> nodeList) {

}


double UnifiedTSHomogeneousTreeLikelihood::updateLikelihoodOnTreeRearrangement(std::vector<tshlib::VirtualNode *> &nodeList) {
    return 0;
}


double UnifiedTSHomogeneousTreeLikelihood::getLogLikelihoodOnTreeRearrangement() const {
    return 0;
}


void UnifiedTSHomogeneousTreeLikelihood::topologyCommitTree() {

    std::vector<tshlib::VirtualNode *> nodelist;
    nodelist = utree_->listVNodes;
    //UtreeBppUtils::treemap tm = treemap_;

    std::map<int, bpp::Node *> tempMap;
    // reset inBtree
    for (auto &bnode:tree_->getNodes()) {

        tempMap.insert(std::pair<int, bpp::Node *>(bnode->getId(), bnode));
        // Empty array of sons on the node
        bnode->removeSons();
        // Empty father connection
        bnode->removeFather();

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

            leftBNode->setDistanceToFather(tree_->getDistanceToFather(leftBNode->getId()));
            rightBNode->setDistanceToFather(tree_->getDistanceToFather(rightBNode->getId()));
            //Add new sons
            pNode->setSon(0, leftBNode);
            pNode->setSon(1, rightBNode);
            pNode->setDistanceToFather(tree_->getDistanceToFather(pNode->getId()));
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
            leftBNode->setDistanceToFather(tree_->getDistanceToFather(leftBNode->getId()));
            rightBNode->setDistanceToFather(tree_->getDistanceToFather(rightBNode->getId()));

            tree_->getRootNode()->setSon(0, leftBNode);
            tree_->getRootNode()->setSon(1, rightBNode);

        }

        //std::cerr << "\t done\n";

    }

}


void UnifiedTSHomogeneousTreeLikelihood::topologyChangeSuccessful(std::vector<tshlib::VirtualNode *> listNodes) {

    // Update BPP tree using the structure in Utree
    topologyCommitTree();

    // Add virtual root to compute the likelihood
    utree_->addVirtualRootNode();

    // remap the virtual nodes to the bpp nodes
    std::vector<Node *> extractionNodes = UtreeBppUtils::remapNodeLists(listNodes, tree_, treemap_);

    // Optimise branches involved in the tree rearrangement
    fireBranchOptimisation(this, extractionNodes);

    // Remove the virtual root to allow for further tree topology improvements
    utree_->removeVirtualRootNode();

}


