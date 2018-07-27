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

    vnode_ = vnode; // tshlib node
    bnode_ = bnode; // bpp node

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

//TODO: change this method to PIPmsaSingle and re-implement for PIPmsaComp
void PIPnode::_setFVleaf(MSA_t &MSA) {

    // get the number of gamma categories
    size_t num_gamma_categories = progressivePIP_->numCatg_;

    // get the number of compressed sites
    int lenComprSeqs = static_cast<PIPmsaSingle *>(MSA_)->pipmsa->rev_map_compressed_seqs_.size();

    // resize fv data([site][catg][states])
    fv_data_.resize(lenComprSeqs);
    for (int i = 0; i < lenComprSeqs; i++) {
        fv_data_[i].resize(num_gamma_categories);
    }

    int idx;
    // go through all the sites
    for (int i = 0; i < lenComprSeqs; i++) {

        // get the index in the compressed map
        idx = static_cast<PIPmsaSingle *>(MSA_)->pipmsa->rev_map_compressed_seqs_.at(i);
        MSAcolumn_t s = MSA.at(idx);

        // allocate fv column to the ext. alphabet size
        bpp::ColMatrix<double> fv;
        fv.resize(progressivePIP_->extendedAlphabetSize_, 1); // ColMatrix as Nx1 matrix
        bpp::MatrixTools::fill(fv, 0.0); // all zeros

        // check if the sequence contains a "forbidden" char
        if(s[0]=='X' || s[0]==' ' || s[0]=='-'){
            LOG(FATAL) << "\nERROR sequence contains either 'X' or ' ' or '-'";
        }

        // get the char position in the alphabet
        idx = progressivePIP_->alphabet_->charToInt(&s[0]);

        // set to 1 the indicator array at the position of the observed char
        fv(idx, 0) = 1.0;

        // assign the indicator array to all the gamma categories
        for (int catg = 0; catg < num_gamma_categories; catg++) {
            fv_data_.at(i).at(catg) = fv;
        }

    }

}

bool PIPnode::_isRootNode(){

    // true if PIPnode is the root node
    // (the last node to be aligned) otherwise false
    if(parent == nullptr){
        return true;
    }else{
        return false;
    }

}

bool PIPnode::_isTerminalNode(){

    // true if PIPnode is a leaf node
    // (a leaf node doesn't need to be aligned) otherwise false
    if(childL == nullptr || childR == nullptr){
        return true;
    }else{
        return false;
    }

}

void PIPnode::_setFVemptyLeaf() {

    // get number of gamma categories
    size_t num_gamma_categories = progressivePIP_->numCatg_;

    // indicator array (all zeros except the observed character)
    bpp::ColMatrix<double> fv;
    fv.resize(progressivePIP_->extendedAlphabetSize_, 1);
    bpp::MatrixTools::fill(fv, 0.0); // all zeros

    // get the gap position in the alphabet
    std::string ch(1, GAP_CHAR);
    int gapIndex = progressivePIP_->alphabet_->charToInt(ch);

    fv(gapIndex, 0) = 1.0; // set gap position to 1

    // for all the gamma categories an array of fv values
    fv_empty_data_.resize(num_gamma_categories);

    // assign the indicator array to all gamma categories
    for (int catg = 0; catg < num_gamma_categories; catg++) {
        fv_empty_data_.at(catg) = fv;
    }

}

void PIPnode::_setFVsigmaLeaf() {

    // get the number of gamma categories
    size_t num_gamma_categories = progressivePIP_->numCatg_;

    // get the length of the compressed input sequences
    int lenComprSeqs = static_cast<PIPmsaSingle *>(MSA_)->getCompressedMSAlength();

    // resize the array ([site][numCatg])
    fv_sigma_.resize(lenComprSeqs);

    double fv0;
    // go through all the sites
    for (int site = 0; site < lenComprSeqs; site++) {

        fv_sigma_.at(site).resize(num_gamma_categories);

        // go through all the gamma categories
        for(int catg = 0; catg < num_gamma_categories; catg++) {

            // compute fv_sigma = fv dot pi
            fv0 = MatrixBppUtils::dotProd(fv_data_.at(site).at(catg), progressivePIP_->pi_);

            fv_sigma_.at(site).at(catg) = fv0;
        }
    }

}

void PIPnode::_setFVsigmaEmptyLeaf() {

    // get the number of gamma categories
    size_t num_gamma_categories = progressivePIP_->numCatg_;

    // allocate memory ([numCatg] x 1)
    fv_empty_sigma_.resize(num_gamma_categories);

    double fv0;

    for(int catg = 0; catg < num_gamma_categories; catg++) {

        // compute fv_empty_sigma = fv_empty dot pi
        fv0 = MatrixBppUtils::dotProd(fv_empty_data_.at(catg), progressivePIP_->pi_);

        fv_empty_sigma_.at(catg) = fv0;
    }

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
