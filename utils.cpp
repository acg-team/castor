/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti & Massimo Maiolo
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
 * @file utils.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 21 12 2017
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
#include "utils.hpp"



void ::UtreeBppUtils::_traverseTree(Utree *in_tree, VirtualNode *target, bpp::Node *source) {

    std::string name;

    for(auto &bppNode:source->getSons()){

        auto ichild = new VirtualNode();

        ichild->vnode_id = bppNode->getId();

        if(bppNode->hasName()){
            name = bppNode->getName();
        }else{
            name = "V"+std::to_string(bppNode->getId());
        }

        ichild->vnode_name = name;
        ichild->vnode_branchlength = bppNode->getDistanceToFather();

        if(!bppNode->isLeaf()){

            target->connectNode(ichild);
            in_tree->addMember(ichild);
            _traverseTree(in_tree, ichild, bppNode);

        }else{

            // Set all the other directions to null
            ichild->_setNodeLeft(nullptr);
            ichild->_setNodeRight(nullptr);

            // Set the LEAF flag to true
            ichild->vnode_leaf = true;
            in_tree->addMember(ichild);
            target->connectNode(ichild);

        }

    }


}


void ::UtreeBppUtils::convertUtree(bpp::TreeTemplate<bpp::Node> *in_tree, Utree *out_tree) {

    bpp::Node * RootNode = in_tree->getRootNode();
    std::string name;
    // For each node descending the root, create either a new VirtualNode
    for (auto &bppNode:RootNode->getSons()) {

        auto ichild = new VirtualNode;

        ichild->vnode_id = bppNode->getId();
        if(bppNode->hasName()){
            name = bppNode->getName();
        }else{
            name = "V"+std::to_string(bppNode->getId());
        }
        ichild->vnode_name = name;
        ichild->vnode_branchlength = bppNode->getDistanceToFather();

        _traverseTree(out_tree, ichild, bppNode);

        // Add this node as starting point of the tree
        out_tree->addMember(ichild, true);


    }

    // Collapse multiforcating trees to star tree pointing to the same pseudoroot
    // Pick root node at random within the node-vector

    out_tree->startVNodes.at(0)->_setNodeUp(out_tree->startVNodes.at(1));
    out_tree->startVNodes.at(1)->_setNodeUp(out_tree->startVNodes.at(0));

}

Eigen::MatrixXd MatrixBppUtils::Matrix2Eigen(const bpp::Matrix<double> &inMatrix) {

    size_t rows, cols;

    rows = inMatrix.getNumberOfRows();
    cols = inMatrix.getNumberOfColumns();

    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(rows, cols);

    for(int r=0; r<rows; r++){
        for (int c=0; c<cols; c++){

            m(r,c) = inMatrix(r,c);

        }

    }

    return m;
}

Eigen::VectorXd MatrixBppUtils::Vector2Eigen(std::vector<double> &inVector) {

    Eigen::VectorXd vector = Eigen::VectorXd::Zero(inVector.size());

    for(int r=0; r<inVector.size(); r++){

        vector(r) = inVector.at(r);

    }

    return vector;
}
