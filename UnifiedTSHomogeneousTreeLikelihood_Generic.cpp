/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2018 by Lorenzo Gatti & Massimo Maiolo
 *******************************************************************************
 *
 * MiniJATI is a computer program whose purpose is to infer phylogentic trees
 * and multi-sequence alignments under the Poisson Indel Process.
 *
 * miniJATI is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.
 *
 * This software is based and extends the Bio++ libraries which are
 * developed by the Bio++ Development Team <http://biopp.univ-montp2.fr>
 *
 * This software is governed by the CeCILL  license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 *******************************************************************************/

/**
 * @file UnifiedTSHomogeneousTreeLikelihood.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 09 04 2018
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
#include "UnifiedTSHomogeneousTreeLikelihood_Generic.hpp"

using namespace bpp;

/*
 * Implementation for the interface  likelihood under tree search engines (all canonical models)
 */

UnifiedTSHomogeneousTreeLikelihood::UnifiedTSHomogeneousTreeLikelihood(const Tree &tree,
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
                                                                       bool usePatterns) :
        RHomogeneousTreeLikelihood_Generic(tree, data, model, rDist, checkRooted, verbose, usePatterns), data_(&data), utree_(utree_), treemap_(*tm), usePatterns_
        (usePatterns) {

    setOptimiser(static_cast<UnifiedTSHomogeneousTreeLikelihood *>(this), optNumericalDerivatives, params, suffix, true, verbose, 0);


}

UnifiedTSHomogeneousTreeLikelihood::UnifiedTSHomogeneousTreeLikelihood(const Tree &tree,
                                                                       TransitionModel *model,
                                                                       DiscreteDistribution *rDist,
                                                                       tshlib::Utree *utree_,
                                                                       UtreeBppUtils::treemap *tm,
                                                                       bool optNumericalDerivatives,
                                                                       std::map<std::string, std::string> &params,
                                                                       const std::string &suffix,
                                                                       bool checkRooted,
                                                                       bool verbose,
                                                                       bool usePatterns) :
        RHomogeneousTreeLikelihood_Generic(tree, model, rDist, checkRooted, verbose, usePatterns), utree_(utree_), treemap_(*tm), usePatterns_(usePatterns) {

    setOptimiser(static_cast<UnifiedTSHomogeneousTreeLikelihood *>(this), optNumericalDerivatives, params, suffix, true, verbose, 0);

}


UnifiedTSHomogeneousTreeLikelihood::~UnifiedTSHomogeneousTreeLikelihood() = default;


void UnifiedTSHomogeneousTreeLikelihood::init_(bool usePatterns) {

    likelihoodDataTest_ = new DRASRTreeLikelihoodData(tree_, rateDistribution_->getNumberOfCategories(), usePatterns);    // FV

}


void UnifiedTSHomogeneousTreeLikelihood::fireTopologyChange(std::vector<int> nodeList,  tshlib::Utree &_utree__topology) {

    // Store the nodes where the likelihood should be recomputed in post-order
    for (auto &nodeID:nodeList) {

        Node *node = tree_->getNode(nodeID);

        updateLikelihoodArrays(node, _utree__topology);
    }

}


double UnifiedTSHomogeneousTreeLikelihood::updateLikelihoodOnTreeRearrangement(std::vector<int> &nodeList, tshlib::Utree &_utree__topology) {
    // Add root to the utree structure
    _utree__topology.addVirtualRootNode();

    // 0. convert the list of tshlib::VirtualNodes into bpp::Node
    //std::vector<int> rearrangedNodes = remapVirtualNodeLists(nodeList);

    // 1. Fire topology change
    //fireTopologyChange(rearrangedNodes);
    fireTopologyChange(remapVirtualNodeLists(nodeList),_utree__topology);

    // 2. Compute loglikelihood
    double logLk = getLogLikelihood();

    // Remove root node from the utree structure
    _utree__topology.removeVirtualRootNode();

    return logLk;
}


double UnifiedTSHomogeneousTreeLikelihood::getLogLikelihoodOnTreeRearrangement() const {
    return 0;
}


void UnifiedTSHomogeneousTreeLikelihood::topologyCommitTree() {

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

        //std::cerr << "vnode " << vnode->getNodeName();
        if (!vnode->isTerminalNode()) {

            // get corrisponding sons in inBTree
            bpp::Node *leftBNode = tempMap[treemap_.right.at(vnode->getNodeLeft()->getVnode_id())];
            bpp::Node *rightBNode = tempMap[treemap_.right.at(vnode->getNodeRight()->getVnode_id())];

            // get corrisponding parent in inBTree
            bpp::Node *pNode = tempMap[treemap_.right.at(vnode->getVnode_id())];

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

            bpp::Node *leftBNode = tempMap[treemap_.right.at(vnode->getVnode_id())];
            bpp::Node *rightBNode = tempMap[treemap_.right.at(vnode->getNodeUp()->getVnode_id())];

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



void UnifiedTSHomogeneousTreeLikelihood::topologyChangeSuccessful(std::vector<int> listNodes) {

    // Update BPP tree using the structure in Utree
    topologyCommitTree();

    // The likelihood function components are reset to use the new topology (lk, dlk, d2lk).
    resetLikelihoodsOnTopologyChangeSuccessful();

    // Add virtual root to compute the likelihood
    utree_->addVirtualRootNode();

    // remap the virtual nodes to the bpp nodes
    //std::vector<Node *> extractionNodes = UtreeBppUtils::remapNodeLists(listNodes, tree_, treemap_);

    // Optimise branches involved in the tree rearrangement
    //fireBranchOptimisation(extractionNodes);
    fireBranchOptimisation(UtreeBppUtils::remapNodeLists(listNodes, tree_, treemap_));

    // Remove the virtual root to allow for further tree topology improvements
    utree_->removeVirtualRootNode();

}


std::vector<int> UnifiedTSHomogeneousTreeLikelihood::remapVirtualNodeLists(std::vector<int> &inputList) const {

    std::vector<int> newList;

    for (auto &vnode:inputList) {

        newList.push_back(tree_->getNode(treemap_.right.at(vnode))->getId());
    }

    return newList;
}


void UnifiedTSHomogeneousTreeLikelihood::updateLikelihoodArrays(Node *node, tshlib::Utree &_utree__topology) {

    if (node->isLeaf()) return;

    //size_t nbSites = likelihoodData_->getLikelihoodArray(node->getId()).size();
    size_t nbClasses = getLikelihoodData()->getNumberOfClasses();
    size_t nbStates = getLikelihoodData()->getNumberOfStates();
    size_t nbSites = getLikelihoodData()->getLikelihoodArray(node->getId()).size();

    // Must reset the likelihood array first (i.e. set all of them to 1):
    VVVdouble *_likelihoods_node = &getLikelihoodData()->getLikelihoodArray(node->getId());
    for (size_t i = 0; i < nbSites; i++) {
        //For each site in the sequence,
        VVdouble *_likelihoods_node_i = &(*_likelihoods_node)[i];
        for (size_t c = 0; c < nbClasses; c++) {
            //For each rate classe,
            Vdouble *_likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
            for (size_t x = 0; x < nbStates; x++) {
                //For each initial state,
                (*_likelihoods_node_i_c)[x] = 1.;
            }
        }
    }


    // Get mapped node on Utree
    std::vector<int> sonsIDs(2);

    //tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
    //tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();
    //sonsIDs.push_back(treemap_.right.at(vnode_left));
    //sonsIDs.push_back(treemap_.right.at(vnode_right));

    sonsIDs[0] = treemap_.right.at(_utree__topology.getNodeIdsMap().at(treemap_.left.at(node->getId()))->getNodeLeft()->getVnode_id());
    sonsIDs[1] = treemap_.right.at(_utree__topology.getNodeIdsMap().at(treemap_.left.at(node->getId()))->getNodeRight()->getVnode_id());


    size_t nbNodes = sonsIDs.size();

    for (size_t l = 0; l < nbNodes; l++) {
    //For each son node,

        const Node *son = tree_->getNode(sonsIDs.at(l));
        const Node *parent = son->getFather();
        //patternLinks_
        std::vector<size_t> * _patternLinks_node_son = &getLikelihoodData()->getArrayPositions(parent->getId(), son->getId());
        std::vector<size_t>_patternLinks_node_son_out = getLikelihoodData()->getArrayPositions(parent->getId(), son->getId());
        VVVdouble *_likelihoods_son = &getLikelihoodData()->getLikelihoodArray(son->getId());

        for (size_t i = 0; i < nbSites; i++) {
            //For each site in the sequence,
            VVdouble *_likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_node_son)[i]];
            VVdouble *_likelihoods_node_i = &(*_likelihoods_node)[i];
            VVVdouble pxy__son = getTransitionProbabilitiesPerRateClass(son->getId(), i);

            for (size_t c = 0; c < nbClasses; c++) {
                //For each rate classe,
                Vdouble *_likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
                Vdouble *_likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
                VVdouble *pxy__son_c = &pxy__son[c];

                for (size_t x = 0; x < nbStates; x++) {
                    //For each initial state,
                    Vdouble *pxy__son_c_x = &(*pxy__son_c)[x];
                    double likelihood = 0;
                    for (size_t y = 0; y < nbStates; y++)
                        likelihood += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];

                    (*_likelihoods_node_i_c)[x] *= likelihood;
                }
            }
        }

        //&getLikelihoodData()->getNodeData(&son).reInit();

    }

}

void UnifiedTSHomogeneousTreeLikelihood::resetLikelihoodsOnTopologyChangeSuccessful() {

    initialized_ = false;
    init_(usePatterns_);
    setData(*data_);
    initialize();

};


