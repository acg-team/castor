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
 * @file pPIP.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 19 02 2018
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
 * @see For more information visit:
 */

#include <chrono>
#include <random>

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <glog/logging.h>

#include "progressivePIP.hpp"
#include "FactoryPIPnode.hpp"
#include "CompositePIPnode.hpp"

using namespace bpp;

progressivePIP::progressivePIP(tshlib::Utree *utree,
                               bpp::Tree *tree,
                               bpp::SubstitutionModel *smodel,
                               UtreeBppUtils::treemap &inTreeMap,
                               bpp::SequenceContainer *sequences,
                               bpp::DiscreteDistribution *rDist,
                               long seed) {

    utree_ = utree; // tshlib tree
    tree_ = new TreeTemplate<Node>(*tree);  // bpp tree
    substModel_ = smodel;   // substitution model
    treemap_ = inTreeMap;   // tree map
    sequences_ = sequences; // sequences
    rDist_ = rDist; // rate-variation among site distribution
    alphabet_ = substModel_->getAlphabet(); // alphabet
    alphabetSize_ = alphabet_->getSize() - 1;   // original alphabet size
    extendedAlphabetSize_ = alphabet_->getSize();   // extended alphabet size
    seed_ = seed;   // seed for random number generation

};



void progressivePIP::_setLambda(double lambda) {

    // original lambda w/o rate variation
    lambda0_ = lambda;

    // insertion rate with rate variation among site r
    lambda_.resize(numCatg_);
    for (int i = 0; i < rDist_->getNumberOfCategories(); i++) {
        // lambda(r) = lambda * r
        lambda_.at(i) = lambda * rDist_->getCategories().at(i);
    }

}

void progressivePIP::_setMu(double mu) {

    // checks division by 0 or very small value
    if (fabs(mu) < SMALL_DOUBLE) {
        PLOG(FATAL) << "ERROR: mu is too small";
    }

    // original mu w/o rate variation
    mu0_ = mu;

    // deletion rate with rate variation among site r
    mu_.resize(numCatg_);
    for (int i = 0; i < rDist_->getCategories().size(); i++) {
        // mu(r) = mu *r
        mu_.at(i) = mu * rDist_->getCategories().at(i);
    }

}

void progressivePIP::_setPi(const Vdouble &pi) {

    // copy pi (steady state frequency distribution)
    // pi is a colMatrix (column array) to simplify the matrix multiplication
    pi_.resize(pi.size(), 1);
    for (int i = 0; i < pi.size(); i++) {
        pi_(i, 0) = pi.at(i);
    }

}

void progressivePIP::_setAllIotas() {

    // compute all the insertion probabilities (iota function)
    // N.B. the insertion probability at the root is different
    // from any other internal node

    double T;

    for(auto &node:compositePIPaligner_->pip_nodes_){

        node->iotasNode_.resize(numCatg_);

        if(node->_isRootNode()){

            // formula for the root PIPnode

            for (int catg = 0; catg < numCatg_; catg++) {

                // T(r) = tau + 1/ (mu * r)
                T = tau_ + 1 / mu_.at(catg);

                // checks division by 0 or too small number
                if (fabs(T) < SMALL_DOUBLE) {
                    PLOG(WARNING) << "ERROR in _setAllIotas: T too small";
                }

                // iota(root,r) = (lambda * r)/ (mu * r) / (lambda * r * (tau + 1/ (mu * r) ) )
                //              = 1 / (mu * r) / (tau + 1/ (mu *r) )
                node->iotasNode_.at(catg) =  (1 / mu_.at(catg)) / T;

            }

        }else{

            // formula for an internal PIPnode

            for (int catg = 0; catg < numCatg_; catg++) {

                // T(r) = tau + 1/(mu * r)
                T = tau_ + 1 / mu_.at(catg);

                // checks division by 0 or too small number
                if (fabs(T) < SMALL_DOUBLE) {
                    PLOG(WARNING) << "ERROR in _setAllIotas: T too small";
                }

                //iotasNode_[node->getId()][catg] = (node->getDistanceToFather() * rDist_->getCategory(catg) ) / T;
                // iota(v,r) = ( lambda * r * b(v) ) / (lambda * r * (tau + 1/ (mu *r) ) )
                //           = b(v) / (tau + 1/ (mu *r) )
                node->iotasNode_.at(catg) = node->bnode_->getDistanceToFather() / T;

            }

        }
    }

}

void progressivePIP::_setAllAlphas() {

    // alpha(v) = sum_from_v_to_root ( iota * beta * zeta )
    // zeta = exp(- mu *b ) is the "pure" survival probability

    for(auto &node:compositePIPaligner_->pip_nodes_){

        // resize for each gamma category
        node->alphaNode_.resize(numCatg_);

        // zeta = exp(- mu * b ) is 1 at the starting node
        double zeta = 1.0;
        for(int catg=0;catg<numCatg_;catg++){
            // initialize alpha with the probability at the starting node which is
            // alpha(v) = iota * beta
            node->alphaNode_.at(catg) = node->iotasNode_.at(catg) * node->betasNode_.at(catg) * zeta;
        }

        // climb the tree from the starting node towards the root
        bpp::PIPnode *tmpPIPnode = node;

        double T; // path length from the starting node to the node below the insertion point

        if( !tmpPIPnode->_isRootNode() ) { // root node doesn't have distanceToFather
            T = tmpPIPnode->bnode_->getDistanceToFather();
        }

        while( !tmpPIPnode->_isRootNode() ){

            // get the parent node
            tmpPIPnode = tmpPIPnode->parent;

            for(int catg=0;catg<numCatg_;catg++) {

                // zeta(v) = "pure" survival probability from the starting node
                // and the node below the insertion point
                zeta = exp(-mu_.at(catg)*T);

                // alpha(v) = sum_from_v_to_root ( iota * beta * zeta )
                node->alphaNode_.at(catg) += tmpPIPnode->iotasNode_.at(catg) * tmpPIPnode->betasNode_.at(catg) * zeta;

            }

            if( !tmpPIPnode->_isRootNode() ){
                T += tmpPIPnode->bnode_->getDistanceToFather(); // increase the path length from starting node and root
            }
        }

    }

}

void progressivePIP::_setAllBetas() {

    // compute all the survival probabilities (beta function)
    // N.B. the survival probability at the root is different
    // from any other internal node

    for(auto &node:compositePIPaligner_->pip_nodes_){

        node->betasNode_.resize(numCatg_);

        if(node->_isRootNode()){

            // formula for the root PIPnode

            for (int catg = 0; catg < rDist_->getNumberOfCategories(); catg++) {
                // by definition at the root (ev. local root) the survival probability is 1
                node->betasNode_.at(catg) = 1.0;
            }

        }else{

            // formula for an internal PIPnode

            for (int catg = 0; catg < rDist_->getCategories().size(); catg++) {

                // muT(r) = r * mu * b(v)
                double muT = rDist_->getCategory(catg) * mu_.at(catg) * node->bnode_->getDistanceToFather();

                // checks division by 0 or too small value
                if (fabs(muT) < SMALL_DOUBLE) {
                    perror("ERROR mu * T is too small");
                }
                // survival probability on node v (different from (local)-root)
                // beta(v,r) = (1 - exp( -mu * r * b(v) )) / (mu * r * b(v))
                node->betasNode_.at(catg) = (1.0 - exp(-muT)) / muT;
            }

        }
    }

}

//void progressivePIP::_computeAllLkemptyRec(bpp::PIPnode *node) {
//
//    if(node->_isTerminalNode()){
//
//        node->fv_empty_sigma__.resize(numCatg_);
//
//        for(int catg=0;catg<numCatg_;catg++){
//            node->fv_empty_sigma__.at(catg) = 0.0;
//        }
//
//    }else{
//
//        _computeAllLkemptyRec(node->childL);
//        _computeAllLkemptyRec(node->childR);
//
//        double bL = node->childL->bnode_->getDistanceToFather();
//        double bR = node->childR->bnode_->getDistanceToFather();
//        double zetaL;
//        double zetaR;
//
//        node->fv_empty_sigma__.resize(numCatg_);
//
//        for(int catg=0;catg<numCatg_;catg++){
//
//            zetaL = exp(-mu_.at(catg) * bL);
//            zetaR = exp(-mu_.at(catg) * bR);
//
//            node->fv_empty_sigma__.at(catg) = \
//                    (1 - zetaL) * (1 - zetaR) + \
//                    (1 - zetaL) * zetaR * node->childR->fv_empty_sigma__.at(catg) + \
//                    zetaL * node->childL->fv_empty_sigma__.at(catg)* (1 - zetaR) + \
//                    zetaL * node->childL->fv_empty_sigma__.at(catg) * zetaR * node->childR->fv_empty_sigma__.at(catg);
//
//        }
//
//    }
//
//}

//void progressivePIP::_computeAllFvEmptySigmaRec(bpp::PIPnode *node) {
//
//    if(node->_isTerminalNode()){ // leaf
//
//        node->_setFVemptyLeaf(); // set fv_empty
//
//        node->_setFVsigmaEmptyLeaf(); // set fv_sigma_empty = fv_empty dot pi
//
//    }else{ // internal node
//
//        // recursive call
//        _computeAllFvEmptySigmaRec(node->childL);
//        _computeAllFvEmptySigmaRec(node->childR) ;
//
//        node->_setFVemptyNode(); // set fv_empty
//
//        node->_setFVsigmaEmptyNode(); // set fv_sigma_empty = fv_empty dot pi
//    }
//
//}

void progressivePIP::_computeAllFvEmptySigma() {

    // get the root PIPnode
    bpp::PIPnode *PIPnodeRoot = getPIPnodeRootNode();

    // recursive call
    PIPnodeRoot->_computeAllFvEmptySigmaRec();

}

void progressivePIP::_initializePIP(std::vector<tshlib::VirtualNode *> &list_vnode_to_root,
                                    enumDP3Dversion DPversion,
                                    int num_sb) {

    //***************************************************************************************
    // get dimensions
    numNodes_ = list_vnode_to_root.size();
    numCatg_ = rDist_->getNumberOfCategories();
    //***************************************************************************************
    // computes lambda and mu with gamma
    // set lambdas with rate variation (gamma distribution)
    _setLambda(substModel_->getParameter("lambda").getValue());
    // set mus with rate variation (gamma distribution)
    _setMu(substModel_->getParameter("mu").getValue());
    //***************************************************************************************
    // set Pi
    // local copy of steady state frequency (Pi)
    _setPi(substModel_->getFrequencies());
    //***************************************************************************************

    nodeFactory *nodeFactory = new bpp::nodeFactory(); // Factory pattern for DP3D-CPU, DP3D-RAM,...

    compositePIPaligner_ = new CompositePIPnode(numNodes_); // Composite pattern with array of PIPnodes

    for (auto &vnode:list_vnode_to_root) {

        auto bnode = tree_->getNode(treemap_.right.at(vnode), false); // get bnode from vnode through the tree-map

        // create a PIPnode
        PIPnode * pip_node = nodeFactory->getPIPnode(DPversion, // PIPnode of type CPU, RAM, SB,... to access the correct DP version
                                                     this,      // PIPnode has access to progressivePIP fields through this pointer
                                                     vnode,     // PIPnode store the correponding vnode and
                                                     bnode);    // the bnode
        //***************************************************************************************
        // get Qs
        // set substitution/deletion probabilities with rate variation (gamma distribution)
        pip_node->_getPrFromSubstitutionModel();
        //***************************************************************************************

        compositePIPaligner_->addPIPnode(pip_node); // add PIPnode to composite array of PIPnodes

    }
    //***************************************************************************************

    _buildPIPnodeTree(); // build a tree of PIPnodes

    _computeTau_(); // compute the total tree length and the left/right subtree length
                    // (length of the left/right subtree rooted at the given node) at each PIPnode

    //_computeLengthPathToRoot(); // at each node compute the length of the path from that node to
                                // the root. The length from the root to the root is 0

    _computeNu();   // compute the Poisson normalizing intensity (corresponds to the expected MSA length)

    _setAllIotas(); // set iota (survival probability) on all nodes

    _setAllBetas(); // set beta (survival probability) on all nodes

    _setAllAlphas(); // alpha(v) = sum_from_v_to_root ( iota * beta * zeta )
                     // zeta = exp(- mu *b ) is the "pure" survival probability

    _computeAllFvEmptySigma(); // compute all fv_empty and fv_empty_sigma values

}

void progressivePIP::_buildPIPnodeTree() {

    // build a binary tree of PIPnode that dictates the alignmanet order

    // this method build a PIPnode binary tree with the same structure (same relations)
    // to the bpp tree

    for (auto &node:compositePIPaligner_->pip_nodes_) {

        bpp::Node *bnode = node->bnode_;

        int bnodeId = bnode->getId();

        if(bnode->hasFather()){
            // internal node
            int bnodeFatherId = bnode->getFatherId();
            node->parent = compositePIPaligner_->pip_nodes_.at(bnodeFatherId);
        }else{
            // root node
            node->parent = nullptr;
            PIPnodeRoot = node;
        }

        if(!bnode->isLeaf()){
            // internal node
            std::vector<int> bnodeSonsId = bnode->getSonsId();
            node->childL = compositePIPaligner_->pip_nodes_.at(bnodeSonsId.at(LEFT));
            node->childR = compositePIPaligner_->pip_nodes_.at(bnodeSonsId.at(RIGHT));
        }else{
            // leaf node
            node->childL = nullptr;
            node->childR = nullptr;
        }

    }

}

void progressivePIP::_computeTauRec_(PIPnode *pipnode) {

    // recursive computation of the total tree length and
    // of the total left/right subtree length

    if(pipnode->_isTerminalNode()){

        pipnode->subTreeLenL_ = 0.0;
        pipnode->subTreeLenR_ = 0.0;

    }else{

        _computeTauRec_(pipnode->childL);
        _computeTauRec_(pipnode->childR);

        pipnode->subTreeLenL_ = pipnode->childL->subTreeLenL_ + \
                                pipnode->childL->subTreeLenR_ + \
                                pipnode->childL->bnode_->getDistanceToFather();

        pipnode->subTreeLenR_ = pipnode->childR->subTreeLenL_ + \
                                pipnode->childR->subTreeLenR_ + \
                                pipnode->childR->bnode_->getDistanceToFather();

        tau_ += pipnode->childL->bnode_->getDistanceToFather() +\
                pipnode->childR->bnode_->getDistanceToFather();

    }

}

void progressivePIP::_computeTau_() {

    // compute the total tree length and the length of the left/right subtree
    // of the tree rooted at a given PIPnode

    // get the root Id
    int rootId = getBPProotNode()->getId();

    // get the root node
    PIPnode * rootPIPnode = compositePIPaligner_->pip_nodes_.at(rootId);

    // init the total tree length
    tau_ = 0.0;

    // compute the total tree length by means of a recursive method
    progressivePIP::_computeTauRec_(rootPIPnode);

}

//void progressivePIP::_computeLengthPathToRoot(){
//
//    // get through the list of all the PIPnodes
//    for (auto &node:compositePIPaligner_->pip_nodes_) {
//
//        // get the bnode associated to the PIPnode
//        bpp::Node *bnode = node->bnode_;
//
//        // initialize the path length
//        double T = 0.0;
//
//        // climb the PIPnode tree
//        // the root skips this loop
//        while(bnode->hasFather()){
//
//            // sum the branch length
//            T += bnode->getDistanceToFather();
//
//            // get the parent bnode
//            bnode = bnode->getFather();
//        }
//
//        // save the path length at the given node
//        // the path length from root to root is 0.0
//        node->distanceToRoot = T;
//
//    }
//
//}

void progressivePIP::_computeNu() {

    // compute the normalizing Poisson intensity (expected MSA length)

    // get the number of gamma categories
    nu_.resize(numCatg_);

    for (int catg = 0; catg < numCatg_; catg++) {
        // computes the normalizing constant with discrete rate variation (gamma distribution)
        // nu(r) = lambda * r * (tau + 1/(mu *r))
        nu_.at(catg) = lambda_.at(catg) * (tau_ + 1 / mu_.at(catg));
    }

}

//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************
double progressivePIPutils::add_lns(double a_ln, double b_ln) {
    //ln(a + b) = ln{exp[ln(a) - ln(b)] + 1} + ln(b)

    double R;

    if (std::isinf(a_ln) && std::isinf(b_ln)) {
        R = -std::numeric_limits<double>::infinity();
    } else if (std::isinf(a_ln)) {
        R = b_ln;
    } else if (std::isinf(b_ln)) {
        R = a_ln;
    } else if ((abs(a_ln - b_ln) >= 36.043653389117155)) {
        //TODO:check this
        //2^52-1 = 4503599627370495.	log of that is 36.043653389117155867651465390794
        R = max(a_ln, b_ln);
    } else {
        R = log(exp(a_ln - b_ln) + 1) + b_ln;
    }

    return R;
}
//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************