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
    rDist_ = rDist; // gamma distribution
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

void progressivePIP::initializePIP(std::vector<tshlib::VirtualNode *> &list_vnode_to_root,
                                   int num_sb,
                                   enumDP3Dversion DPversion) {

    //***************************************************************************************
    // GET DIMENSIONS
    numNodes_ = list_vnode_to_root.size();
    numCatg_ = rDist_->getNumberOfCategories();
    //***************************************************************************************
    // COMPUTES LAMBDA AND MU WITH GAMMA
    // set lambdas with rate variation (gamma distribution)
    _setLambda(substModel_->getParameter("lambda").getValue());
    // set mus with rate variation (gamma distribution)
    _setMu(substModel_->getParameter("mu").getValue());
    //***************************************************************************************
    // SET PI
    // local copy of steady state frequency (Pi)
    _setPi(substModel_->getFrequencies());
    //***************************************************************************************

    nodeFactory *nodeFactory = new bpp::nodeFactory(); // Factory pattern for DP3D-CPU, DP3D-RAM,...

    compositePIPaligner_ = new CompositePIPnode(numNodes_); // Composite pattern with array of PIPnodes

    for (auto &vnode:list_vnode_to_root) {

        auto bnode = tree_->getNode(treemap_.right.at(vnode), false);

        PIPnode * pip_node = nodeFactory->getPIPnode(DPversion,this,vnode,bnode);

        pip_node->_reserve(numCatg_); // allocate memory

        //***************************************************************************************
        // GET Qs
        // set substitution/deletion probabilities with rate variation (gamma distribution)
        pip_node->_getPrFromSubstitutionModel();
        //***************************************************************************************

        compositePIPaligner_->addPIPnode(pip_node); // add PIPnode to composite array of PIPnodes

    }
    //***************************************************************************************

    // TODO:
    // - costruire albero PIPnode
    // - calcolare tau e nu globali

    //***************************************************************************************
    // COMPUTE SUBTREE LENGHTS
    //***************************************************************************************
    //pip_node->_computeSubtreeTau();
    //pip_node->_computeSubtreeTau();
    //***************************************************************************************
    // COMPUTE ALL LOCAL POISSON NORMALIZING CONSTANTS
    //***************************************************************************************
    //pip_node->_computeNu(numCatg_);

}

//void progressivePIP::_computeLocalTau() {
//
//    if(vnode_->isTerminalNode()){
//        tau_ = 0.0;
//    }else{
//
////        b0 = tree_->getNode(treemap_.right.at(vnode->getNodeLeft()), false)->getDistanceToFather() +\
//                tree_->getNode(treemap_.right.at(vnode->getNodeRight()), false)->getDistanceToFather();
//
//        //tau = _computeTauRecursive(vnode_);
//    }
//
//
//}

void progressivePIP::_buildPIPnodeTree() {

    for (auto &node:compositePIPaligner_->pip_nodes_) {

        bpp::Node *bnode = node->bnode_;

        int bnodeId = bnode->getId();

        if(bnode->hasFather()){
            int bnodeFatherId = bnode->getFatherId();
            node->parent = compositePIPaligner_->pip_nodes_.at(bnodeFatherId);
        }else{
            node->parent = nullptr;
        }

        if(!bnode->isLeaf()){
            std::vector<int> bnodeSonsId = bnode->getSonsId();
            node->childL = compositePIPaligner_->pip_nodes_.at(bnodeSonsId.at(LEFT));
            node->childR = compositePIPaligner_->pip_nodes_.at(bnodeSonsId.at(RIGHT));
        }else{
            node->childL = nullptr;
            node->childR = nullptr;
        }

    }

}

void progressivePIP::_computeTauRec_(PIPnode *pipnode) {

    if(pipnode== nullptr){

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

   int rootId = getRootNode()->getId();

   PIPnode * rootPIPnode = compositePIPaligner_->pip_nodes_.at(rootId);

   tau_ = 0.0;
   progressivePIP::_computeTauRec_(rootPIPnode);

}

//void progressivePIP::_computeTau() {
//
//    tau_ = 0.0;
//    for (auto &node:compositePIPaligner_->pip_nodes_) {
//
//        tau_ += node->bnode_->getDistanceToFather();
//
//    }
//
//}

void progressivePIP::_computeNu() {

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