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
#include "PIPnode.hpp"
#include "CompositePIPnode.hpp"

#define ERR_STATE (-999)
#define DBL_EPSILON std::numeric_limits<double>::min()
#define MATCH_STATE 1
#define GAP_X_STATE 2
#define GAP_Y_STATE 3
#define STOP_STATE 4
#define LEFT 0
#define RIGHT 1

using namespace bpp;

PIPnode::PIPnode(const progressivePIP *pPIP,
                 tshlib::VirtualNode *vnode,
                 bpp::Node *bnode){

    progressivePIP_ = pPIP;

    vnode_ = vnode;
    bnode_ = bnode;

    nodeID_ = bnode->getId();

}

void PIPnode::_reserve(int numCatg){

    // fv values for each gamma category
    fv_empty_data_.resize(numCatg);

    // fv values dot Pi for each gamma category
    fv_empty_sigma_.resize(numCatg);

    // insertion probabilities at the given PIPnode with rate variation (gamma)
    iotasNode_.resize(numCatg);

    // survival probabilities at the given PIPnode with rate variation (gamma)
    betasNode_.resize(numCatg);

    // substitution/deletion probability matrices at the given PIPnode with rate variation (gamma)
    prNode_.resize(numCatg);

};

void PIPnode::_setFVleaf(MSA_t *MSA) {

    int idx;

    size_t num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

    int lenComprSeqs;// = rev_map_compressed_seqs_.size();

    fv_data_.resize(lenComprSeqs);

    for (int i = 0; i < lenComprSeqs; i++) {
        fv_data_[i].resize(num_gamma_categories);
    }

    for (int i = 0; i < lenComprSeqs; i++) {

        idx;// = rev_map_compressed_seqs_.at(0).at(i);
        MSAcolumn_t s = MSA->at(idx);

        bpp::ColMatrix<double> fv;
        fv.resize(progressivePIP_->extendedAlphabetSize_, 1);
        bpp::MatrixTools::fill(fv, 0.0);

        idx = progressivePIP_->alphabet_->charToInt(&s[0]);
        idx = idx < 0 ? progressivePIP_->alphabetSize_ : idx;

        fv(idx, 0) = 1.0;
        for (int catg = 0; catg < num_gamma_categories; catg++) {
            fv_data_.at(i).at(catg) = fv;
        }

    }

}

void PIPnode::_setFVemptyLeaf() {

    size_t num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

    bpp::ColMatrix<double> fv;
    fv.resize(progressivePIP_->extendedAlphabetSize_, 1);
    bpp::MatrixTools::fill(fv, 0.0);
    fv(progressivePIP_->alphabetSize_, 0) = 1.0;

    fv_empty_data_.resize(num_gamma_categories);

    for (int catg = 0; catg < num_gamma_categories; catg++) {
        fv_empty_data_.at(catg) = fv;
    }

}

void PIPnode::_setFVsigmaLeaf() {

    size_t num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

    int lenComprSeqs;// = rev_map_compressed_seqs_.size();

    fv_sigma_.resize(lenComprSeqs);

    double fv0;
    for (int site = 0; site < lenComprSeqs; site++) {

        fv_sigma_.at(site).resize(num_gamma_categories);

        for(int catg = 0; catg < num_gamma_categories; catg++) {

            fv0 = MatrixBppUtils::dotProd(fv_data_.at(site).at(catg), progressivePIP_->pi_);

            fv_sigma_.at(site).at(catg) = fv0;
        }
    }

}

void PIPnode::_setFVsigmaEmptyLeaf() {

    size_t num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

    fv_empty_sigma_.resize(num_gamma_categories);

    double fv0;

    for(int catg = 0; catg < num_gamma_categories; catg++) {

        fv0 = MatrixBppUtils::dotProd(fv_empty_data_.at(catg), progressivePIP_->pi_);

        fv_empty_sigma_.at(catg) = fv0;
    }

}

double PIPnode::_computeSubtreeTau() {

//    if (vnode_->isTerminalNode()) {
//        subTreeLenL_ = 0.0;
//        subTreeLenR_ = 0.0;
////        tau_ = 0.0;
//    }else{
//        subTreeLenL_ = childL->tau_ + childL->bnode_->getDistanceToFather();
//        subTreeLenR_ = childR->tau_ + childR->bnode_->getDistanceToFather();
//        tau_ = subTreeLenL_ + subTreeLenR_;
//    }

}

void PIPnode::_getPrFromSubstitutionModel() {

    if (!bnode_->hasFather()) {
        // root PIPnode doesn't have Pr
    } else {

        for (int i = 0; i < progressivePIP_->rDist_->getNumberOfCategories(); i++) {
            // substitution/deletion probabilities with rate variation (gamma)
            // Pr = exp( branchLength * rateVariation * Q )
            double brlen = bnode_->getDistanceToFather();

            if(brlen<SMALL_DOUBLE){
                LOG(FATAL) << "\nERROR branch length too small";
            }

            prNode_.at(i) = progressivePIP_->substModel_->getPij_t(brlen * \
                            progressivePIP_->rDist_->getCategory(i));
        }
    }

}

void PIPnode::DP3D_PIP_leaf() {

/*

    //*******************************************************************************
    // ALIGNS LEAVES
    //*******************************************************************************
    std::string seqname = sequences_->getSequencesNames().at((int) vnode->vnode_seqid);

    // associates the sequence name to the leaf node
    _setMSAsequenceNames(node, seqname);

    // creates a column containing the sequence associated to the leaf node
    _setMSAleaves(node, sequences_->getSequence(seqname).toString());

    // compresses sequence at the leaves
    map_compressed_seqs_.at(node->getId()).resize(1);
    rev_map_compressed_seqs_.at(node->getId()).resize(1);
    _compressMSA(node, 0);
*/
    // computes the indicator values (fv values) at the leaves
    //_setFVleaf();

    // computes the indicator value for an empty column at the leaf
    _setFVemptyLeaf();

    // computes dotprod(pi,fv)
    _setFVsigmaLeaf();

    // computes dotprod(pi,fv) for an empty column at the leaf
    _setFVsigmaEmptyLeaf();

    /*
    // sets th etraceback path at the leaf
    _setTracebackPathleaves(node);

    _setSubMSAindexLeaves(node);

     */

}
