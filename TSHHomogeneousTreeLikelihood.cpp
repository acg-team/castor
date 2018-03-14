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
 * @file TSHLIBHomogeneousTreeLikelihood.cpp
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
#include "TSHHomogeneousTreeLikelihood.hpp"

using namespace bpp;


TSHHomogeneousTreeLikelihood::TSHHomogeneousTreeLikelihood(AbstractHomogeneousTreeLikelihood *lk,
                                                           const SiteContainer &data,
                                                           TransitionModel *model,
                                                           DiscreteDistribution *rDist,
                                                           tshlib::Utree *inUtree,
                                                           UtreeBppUtils::treemap &inTreeMap)
throw(Exception) :
        RHomogeneousTreeLikelihood(lk->getTree(), model, rDist, false, false),
        optimiser_(0),
        utree_(inUtree),
        treemap_(inTreeMap) {

    // Import the likelihood function
    likelihoodFunc_ = lk;

    // Initialise the optimiser
    optimiser_ = new BfgsMultiDimensions(lk);
    optimiser_->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimiser_->setProfiler(new StdOut);
    optimiser_->setMessageHandler(new StdOut);
    optimiser_->setVerbose(0);

}

TSHHomogeneousTreeLikelihood::TSHHomogeneousTreeLikelihood(const TSHHomogeneousTreeLikelihood &lik) :
        RHomogeneousTreeLikelihood(lik),
        optimiser_(0),
        utree_(),
        treemap_(lik.getTreeMap()) {

    optimiser_ = dynamic_cast<BfgsMultiDimensions *>(lik.optimiser_->clone());
};

TSHHomogeneousTreeLikelihood &TSHHomogeneousTreeLikelihood::operator=(const TSHHomogeneousTreeLikelihood &lik) {
    RHomogeneousTreeLikelihood::operator=(lik);
    if (optimiser_) delete optimiser_;
    optimiser_ = dynamic_cast<BfgsMultiDimensions *>(lik.optimiser_->clone());

    return *this;
}

TSHHomogeneousTreeLikelihood::~TSHHomogeneousTreeLikelihood() {

    delete optimiser_;
}

AbstractHomogeneousTreeLikelihood *TSHHomogeneousTreeLikelihood::getLikelihoodFunction() const {
    return likelihoodFunc_;
}

void TSHHomogeneousTreeLikelihood::fixTopologyChanges(tshlib::Utree *inUTree) {

    std::vector<tshlib::VirtualNode *> nodelist;
    nodelist = inUTree->listVNodes;
    UtreeBppUtils::treemap tm = getTreeMap();

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
            bpp::Node *leftBNode = tempMap[tm.right.at(vnode->getNodeLeft())];
            bpp::Node *rightBNode = tempMap[tm.right.at(vnode->getNodeRight())];

            // get corrisponding parent in inBTree
            bpp::Node *pNode = tempMap[tm.right.at(vnode)];

            leftBNode->setFather(pNode);
            rightBNode->setFather(pNode);

            leftBNode->setDistanceToFather(likelihoodFunc_->getTree().getDistanceToFather(leftBNode->getId()));
            rightBNode->setDistanceToFather(likelihoodFunc_->getTree().getDistanceToFather(rightBNode->getId()));
            //Add new sons
            pNode->setSon(0, leftBNode);
            pNode->setSon(1, rightBNode);
            pNode->setDistanceToFather(likelihoodFunc_->getTree().getDistanceToFather(pNode->getId()));
            //std::cerr << "\t internal";

        } else {

            //std::cerr << "\t leaf";

        }
        // in case the current vnode is also the pseudo-root
        if (vnode == vnode->getNodeUp()->getNodeUp()) {
            //std::cerr << "\tvnode pseudoroot";

            bpp::Node *leftBNode = tempMap[tm.right.at(vnode)];
            bpp::Node *rightBNode = tempMap[tm.right.at(vnode->getNodeUp())];

            tree_->getRootNode()->removeSons();

            leftBNode->setFather(tree_->getRootNode());
            rightBNode->setFather(tree_->getRootNode());
            leftBNode->setDistanceToFather(likelihoodFunc_->getTree().getDistanceToFather(leftBNode->getId()));
            rightBNode->setDistanceToFather(likelihoodFunc_->getTree().getDistanceToFather(rightBNode->getId()));

            tree_->getRootNode()->setSon(0, leftBNode);
            tree_->getRootNode()->setSon(1, rightBNode);

        }

        //std::cerr << "\t done\n";

    }

}

void TSHHomogeneousTreeLikelihood::optimiseBranches(std::vector<tshlib::VirtualNode *> listNodes) {

    // remap the virtual nodes to the bpp nodes
    std::vector<Node *> extractionNodes = UtreeBppUtils::remapNodeLists(listNodes, tree_, getTreeMap());

    ParameterList parameters;
    // For each node involved in the move, get the corrisponding branch parameter (no root)
    for (auto &bnode:extractionNodes) {
        if (bnode->hasFather()) {
            Parameter brLen = likelihoodFunc_->getParameter("BrLen" + TextTools::toString(bnode->getId()));
            brLen.setName("BrLen" + TextTools::toString(bnode->getId()));
            parameters.addParameter(brLen);
        }
    }

    // set parameters on the likelihood function (inherited)
    likelihoodFunc_->setParameters(parameters);

    // Re-estimate branch length:
    optimiser_->init(parameters);
    optimiser_->doInit(parameters);
    optimiser_->getStopCondition()->setTolerance(0.001);
    optimiser_->optimize();

    // set parameters on the likelihood function (inherited)
    likelihoodFunc_->setParameters(optimiser_->getParameters());


}

double TSHHomogeneousTreeLikelihood::updateLikelihood(std::vector<tshlib::VirtualNode *> &nodeList) {

    // check if the likelihood function has been instanziated on the PIP model.
    bpp::RHomogeneousTreeLikelihood_PIP *flk_pip;
    bpp::RHomogeneousTreeLikelihood *flk;
    bool indelModel = false;

    if (dynamic_cast<bpp::RHomogeneousTreeLikelihood_PIP *>(getLikelihoodFunction())) {
        flk_pip = dynamic_cast<bpp::RHomogeneousTreeLikelihood_PIP *>(getLikelihoodFunction());
        indelModel = true;
    }

    // Add root to the utree structure
    utree_->addVirtualRootNode();

    // 0. convert the list of tshlib::VirtualNodes into bpp::Node
    std::vector<Node *> rearrangedNodes = remapVirtualNodeLists(nodeList);

    // 1. Recombine FV arrays after move
    if (indelModel) {
        flk_pip->fireTopologyChange(rearrangedNodes);
    } else {
        flk = dynamic_cast<bpp::RHomogeneousTreeLikelihood *>(getLikelihoodFunction());
        likelihoodData_ = flk->getLikelihoodData();
        for (auto &node:rearrangedNodes) {
            updateLikelihoodArrays(node);
        }
    }

    double logLk;

    // 2. Compute the log likelihood with the updated components
    if (indelModel) {
        logLk = flk_pip->getLogLikelihoodOnTopologyChange();
    } else {
        logLk = flk->getLogLikelihood();
    }

    // Remove root node from the utree structure
    utree_->removeVirtualRootNode();


    return logLk;
}

void TSHHomogeneousTreeLikelihood::updateLikelihoodArrays(const Node *node) {

    if (node->isLeaf()) return;

    auto *flk = dynamic_cast<bpp::RHomogeneousTreeLikelihood *>(getLikelihoodFunction());

    //size_t nbSites = likelihoodData_->getLikelihoodArray(node->getId()).size();
    size_t nbSites = likelihoodData_->getNumberOfDistinctSites();
    size_t nbClasses = likelihoodData_->getNumberOfClasses();
    size_t nbStates = likelihoodData_->getNumberOfStates();

    // Must reset the likelihood array first (i.e. set all of them to 1):
    VVVdouble *_likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());
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
    std::vector<int> sonsIDs;
    tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
    tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

    sonsIDs.push_back(treemap_.right.at(vnode_left));
    sonsIDs.push_back(treemap_.right.at(vnode_right));
    size_t nbNodes = sonsIDs.size();

    for (size_t l = 0; l < nbNodes; l++) {
        //For each son node,

        const Node *son = tree_->getNode(sonsIDs.at(l));

        VVVdouble *_likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

        for (size_t i = 0; i < nbSites; i++) {
            //For each site in the sequence,
            VVdouble *_likelihoods_son_i = &(*_likelihoods_son)[i];
            VVdouble *_likelihoods_node_i = &(*_likelihoods_node)[i];
            VVVdouble pxy__son = flk->getTransitionProbabilitiesPerRateClass(son->getId(), i);

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
    }

}

const UtreeBppUtils::treemap &TSHHomogeneousTreeLikelihood::getTreeMap() const {
    return treemap_;
}
