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
    UtreeBppUtils::treemap &tm = getTreeMap();

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
