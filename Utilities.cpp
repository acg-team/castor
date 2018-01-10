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
#include "Utilities.hpp"



void ::UtreeBppUtils::_traverseTree_b2u(Utree *in_tree, VirtualNode *target, bpp::Node *source) {

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
            _traverseTree_b2u(in_tree, ichild, bppNode);

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


void ::UtreeBppUtils::convertTree_b2u(bpp::TreeTemplate<bpp::Node> *in_tree, Utree *out_tree) {

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

        _traverseTree_b2u(out_tree, ichild, bppNode);

        // Add this node as starting point of the tree
        out_tree->addMember(ichild, true);


    }

    // Collapse multiforcating trees to star tree pointing to the same pseudoroot
    // Pick root node at random within the node-vector

    out_tree->startVNodes.at(0)->_setNodeUp(out_tree->startVNodes.at(1));
    out_tree->startVNodes.at(1)->_setNodeUp(out_tree->startVNodes.at(0));

}


bpp::TreeTemplate<bpp::Node>* UtreeBppUtils::convertTree_u2b(tshlib::Utree *in_tree) {

    auto *RootNode = new bpp::Node;

    RootNode->setName(in_tree->rootnode->getNodeName());
    RootNode->setDistanceToFather(in_tree->rootnode->vnode_branchlength);

    _traverseTree_u2b(RootNode, in_tree->rootnode->getNodeLeft());
    _traverseTree_u2b(RootNode, in_tree->rootnode->getNodeRight());


    RootNode->setId((int) in_tree->listVNodes.size());

    // Set root node on new bpp tree
    auto *tree = new bpp::TreeTemplate<bpp::Node>();

    tree->setRootNode(RootNode);

    return tree;
}

void UtreeBppUtils::_traverseTree_u2b(bpp::Node *target, tshlib::VirtualNode *source) {

    auto *child = new bpp::Node;

    child->setName(source->getNodeName());
    child->setId(source->vnode_id);
    child->setDistanceToFather(source->vnode_branchlength);

    if(!source->isTerminalNode()) {

        _traverseTree_u2b(child, source->getNodeLeft());

        _traverseTree_u2b(child, source->getNodeRight());

    }

    target->addSon(child);


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

Eigen::VectorXd MatrixBppUtils::Vector2Eigen(const std::vector<double> &inVector) {

    Eigen::VectorXd vector = Eigen::VectorXd::Zero(inVector.size());

    for(int r=0; r<inVector.size(); r++){

        vector(r) = inVector.at(r);

    }

    return vector;
}

/*
bpp::Matrix<double> MatrixBppUtils::Eigen2Matrix(Eigen::MatrixXd &inMatrix) {

    bpp::Matrix<double> outMatrix;
    outMatrix.resize((size_t) inMatrix.cols(), (size_t) inMatrix.cols());

    for(int r=0; r<inMatrix.rows(); r++){
        for (int c=0; c<inMatrix.cols(); c++){

            outMatrix(r,c) = inMatrix(r,c);

        }

    }


    return outMatrix;
}
*/
