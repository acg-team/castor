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
#include <random>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <glog/logging.h>

#include "progressivePIP.hpp"
#include "FactoryPIPnode.hpp"
#include "CompositePIPmsa.hpp"

using namespace bpp;

#define EARLY_STOP_THR 10

bool nodeRAM::_index_of_max(double m,
                            double x,
                            double y,
                            double epsilon,
                            std::default_random_engine &generator,
                            std::uniform_real_distribution<double> &distribution,
                            int &index,
                            double &val) {

    // get the index and the max value among the three input values (m,x,y)

    // m:  match value
    // x:  gapx value
    // y:   gapy value
    // epsilon: small number for the comparison between to numbers
    // generator: random number generator (when two or three numbers have the same value)
    // distribution: uniform distribution
    // index: index of max (1: MATCH, 2: GAPX, 3: GAPY)
    // val: max value between the three (m,x,y)

    double random_number;

    if (std::isinf(m) & std::isinf(x) & std::isinf(y)){
        // if the three values are -inf than this cell is marked as
        // non-valid (STOP_STATE) and the max val is -inf
        index = int(STOP_STATE);
        val = -std::numeric_limits<double>::infinity();
        return true;
    }

    if (not(std::isinf(m)) & not(std::isinf(x)) & (fabs((m - x)) < epsilon)) {
        // m and x are both not -inf
        // they are identical (their difference is smaller than epsilon)
        x = m; // x is exactly equal to m
    }

    if (not(std::isinf(m)) & not(std::isinf(y)) & (fabs((m - y)) < epsilon)) {
        // y and m are both not -inf
        // they are identical (their difference is smaller than epsilon)
        y = m; // y is exactly equal to m
    }

    if (not(std::isinf(x)) & not(std::isinf(y)) & (fabs((x - y)) < epsilon)) {
        // y and x are both not -inf
        // they are identical (their difference is smaller than epsilon)
        y = x; // y is exactly equal to x
    }

    if (m > x) {
        if (m > y) {
            index = int(MATCH_STATE);
            val = m;
            return true;
        } else if (y > m) {
            index = int(GAP_Y_STATE);
            val = y;
            return true;
        } else {
            if (abs(m - y) < epsilon) {
                //m or y
                random_number = distribution(generator);
                // m and y are equal and have the same value,
                // the state is selected with a uniform random
                // distribution with 50% probability each
                if (random_number < (1.0 / 2.0)) {
                    index = int(MATCH_STATE);
                    val = m;
                    return true;
                } else {
                    index = int(GAP_Y_STATE);
                    val = y;
                    return true;
                }
            } else {
                LOG(FATAL) << "\nSomething went wrong during the comparison in function "
                              "pPIP::_index_of_max. Check call stack below.";
                return false;
            }
        }
    } else if (x > m) {
        if (x > y) {
            index = int(GAP_X_STATE);
            val = x;
            return true;
        } else if (y > x) {
            index = int(GAP_Y_STATE);
            val = y;
            return true;
        } else {
            if (abs(x - y) < epsilon) {
                //x or y
                random_number = distribution(generator);
                // x and y are equal and have the same value,
                // the state is selected with a uniform random
                // distribution with 50% probability each
                if (random_number < (1.0 / 2.0)) {
                    index = int(GAP_X_STATE);
                    val = x;
                    return true;
                } else {
                    index = int(GAP_Y_STATE);
                    val = y;
                    return true;
                }
            } else {
                LOG(FATAL) << "\nSomething went wrong during the comparison in function "
                              "pPIP::_index_of_max. Check call stack below.";
                return false;
            }
        }
    } else {

        double mx = x;
        if (mx > y) {
            //m or x
            random_number = distribution(generator);
            // m and x are equal and have the same value,
            // the state is selected with a uniform random
            // distribution with 50% probability each
            if (random_number < (1.0 / 2.0)) {
                index = int(MATCH_STATE);
                val = m;
                return true;
            } else {
                index = int(GAP_X_STATE);
                val = x;
                return true;
            }
        } else if (y > mx) {
            index = int(GAP_Y_STATE);
            val = y;
            return true;
        } else {
            if (abs(mx - y) < epsilon) {
                //m or x or y
                // m,x and y are equal and have the same value,
                // the state is selected with a uniform random
                // distribution with 1/3 probability each
                random_number = distribution(generator);
                if (random_number < (1.0 / 3.0)) {
                    index = int(MATCH_STATE);
                    val = m;
                    return true;
                } else if (random_number < (2.0 / 3.0)) {
                    index = int(GAP_X_STATE);
                    val = x;
                    return true;
                } else {
                    index = int(GAP_Y_STATE);
                    val = y;
                    return true;
                }
            } else {
                LOG(FATAL) << "\nSomething went wrong during the comparison in function "
                              "pPIP::_index_of_max. Check call stack below.";
                return false;
            }
        }
    }

}

double nodeRAM::max_of_three(double m,
                             double x,
                             double y,
                             double epsilon) {

    // get max value among the three input values (m,x,y)

    if (fabs(m) < epsilon) { // if a cell has been initialized to 0 then its lk is -inf
        m = -std::numeric_limits<double>::infinity();
    }
    if (fabs(x) < epsilon) { // if a cell has been initialized to 0 then its lk is -inf
        x = -std::numeric_limits<double>::infinity();
    }
    if (fabs(y) < epsilon) { // if a cell has been initialized to 0 then its lk is -inf
        y = -std::numeric_limits<double>::infinity();
    }

    if (std::isinf(m) && std::isinf(x) && std::isinf(y)) {
        // if all the three values are -inf then the max value is also -inf
        return -std::numeric_limits<double>::infinity();
    }

    if (m > x) {
        if (m > y) {
            // m greater than x and y
            return m;
        }
        // m greater than x but y greater than m
        return y;
    } else {
        if (x > y) {
            // x greater than m and y
            return x;
        }
        // x greater than m but y greater than x
        return y;
    }

}

//void nodeRAM::_compress_Fv(std::vector<std::vector<double>> &fv_sigma_not_compressed,
//                           std::vector<vector<bpp::ColMatrix<double> > > &fv_data_not_compressed){
//
//    // compress an array of fv values and fv_sigma values
//
//    int comprMSAlen = static_cast<PIPmsaSingle *>(MSA_)->pipmsa->rev_map_compressed_seqs_.size();
//
//    int id_map;
//
//    fv_data_.resize(comprMSAlen);
//
//    fv_sigma_.resize(comprMSAlen);
//
//    for(int i=0;i<comprMSAlen;i++){
//        id_map = static_cast<PIPmsaSingle *>(MSA_)->pipmsa->rev_map_compressed_seqs_.at(i);
//
//        fv_data_.at(i)=fv_data_not_compressed.at(id_map);
//
//        fv_sigma_.at(i)=fv_sigma_not_compressed.at(id_map);
//    }
//
//}

//void nodeRAM::_computeLK_M(std::vector<bpp::ColMatrix<double> > &fvL,
//                             std::vector<bpp::ColMatrix<double> > &fvR,
//                             std::vector<bpp::ColMatrix<double> > &Fv_M_ij,
//                             std::vector<double> &Fv_sigma_M_ij,
//                             double &pr_match,
//                             double &pr_match_full_path) {
//
//
//    // COMPUTE THE MATCH LK
//
//    // fvL: fv array of the left child
//    // fvR: fv array of the right child
//    // Fv_M_ij: result of Fv_M_ij = Pr_R * fv_R hadamard_prod Pr_L * fv_L
//    // Fv_sigma_M_ij: result of Fv_M_ij dot pi
//    // pr_match: match probability (stored for the next layer)
//    // pr_match_full_path: full match probability (used at this layer)
//    // it encompasses the probability of an insertion along
//    // the whole path between the root and this node
//
//    // number of discrete gamma categories
//    int num_gamma_categories = progressivePIP_->numCatg_;
//
//    pr_match = 0.0;
//    pr_match_full_path = 0.0;
//    for (int catg = 0; catg < num_gamma_categories; catg++) {
//
//        // PrfvL = Pr_L * fv_L
//        bpp::ColMatrix<double> PrfvL;
//        bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL.at(catg), PrfvL);
//
//        // PrfvR = Pr_R * fv_R
//        bpp::ColMatrix<double> PrfvR;
//        bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR.at(catg), PrfvR);
//
//        // fv = PrfvL * PrfvR
//        bpp::ColMatrix<double> fv;
//        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);
//
//        Fv_M_ij[catg] = fv;
//
//        // fv0 = pi * fv
//        double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);
//
//        Fv_sigma_M_ij.at(catg) = fv0;
//
//        // match probability with gamma
//        double p = progressivePIP_->rDist_->getProbability((size_t) catg) * \
//               iotasNode_[catg] * betasNode_[catg] * fv0;
//
//        pr_match_full_path += progressivePIP_->rDist_->getProbability((size_t) catg) * \
//                              alphaNode_.at(catg) * fv0;
//
//        pr_match += p;
//    }
//
//}
//
//void nodeRAM::_computeLK_X(std::vector<bpp::ColMatrix<double> > &fvL,
//                             std::vector<bpp::ColMatrix<double> > &fvR,
//                             std::vector<bpp::ColMatrix<double> > &Fv_X_ij,
//                             std::vector<double> &Fv_sigma_X_ij,
//                             double &pr_gapx,
//                             double &pr_gapx_full_path) {
//
//    // COMPUTE THE GAPX LK
//
//    // fvL: fv array of the left child
//    // fvR: fv array of the right child
//    // Fv_X_ij: result of Fv_X_ij = Pr_R * fv_R hadamard_prod Pr_L * fv_L
//    // Fv_sigma_X_ij: result of Fv_sigma_X_ij dot pi
//    // pr_gapx: gapx probability (stored for the next layer)
//    // pr_gapx_full_path: full gapx probability (used at this layer)
//    // it encompasses the probability of an insertion along
//    // the whole path between the root and this node
//
//    // number of discrete gamma categories
//    int numCatg = progressivePIP_->numCatg_;
//
//    pr_gapx = 0.0;
//    pr_gapx_full_path = 0.0;
//    for (int catg = 0; catg < numCatg; catg++) {
//
//        // PrfvL = Pr_L * fv_L
//        bpp::ColMatrix<double> PrfvL;
//        bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL.at(catg), PrfvL);
//
//        // PrfvR = Pr_R * fv_R
//        bpp::ColMatrix<double> PrfvR;
//        bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR.at(catg), PrfvR);
//
//        // fv = PrfvL * PrfvR
//        bpp::ColMatrix<double> fv;
//        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);
//
//        Fv_X_ij[catg] = fv;
//
//        // fv0 = pi * fv
//        double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);
//
//        Fv_sigma_X_ij.at(catg) = fv0;
//
//        // gapX probability with gamma
//        double p0 = progressivePIP_->rDist_->getProbability((size_t) catg) * \
//                    iotasNode_[catg] * betasNode_[catg] * fv0;
//
//        pr_gapx_full_path += progressivePIP_->rDist_->getProbability((size_t) catg) * \
//                              alphaNode_.at(catg) * fv0;
//
//        pr_gapx += p0;
//    }
//
//}
//
//void nodeRAM::_computeLK_Y(std::vector<bpp::ColMatrix<double> > &fvL,
//                             std::vector<bpp::ColMatrix<double> > &fvR,
//                             std::vector<bpp::ColMatrix<double> > &Fv_Y_ij,
//                             std::vector<double> &Fv_sigma_Y_ij,
//                             double &pr_gapy,
//                             double &pr_gapy_full_path) {
//
//    // COMPUTE THE GAPY LK
//
//    // fvL: fv array of the left child
//    // fvR: fv array of the right child
//    // Fv_Y_ij: result of Fv_Y_ij = Pr_R * fv_R hadamard_prod Pr_L * fv_L
//    // Fv_sigma_Y_ij: result of Fv_sigma_Y_ij dot pi
//    // pr_gapy: gapx probability (stored for the next layer)
//    // pr_gapy_full_path: full gapx probability (used at this layer)
//    // it encompasses the probability of an insertion along
//    // the whole path between the root and this node
//
//    // number of discrete gamma categories
//    int numCatg = progressivePIP_->numCatg_;
//
//    pr_gapy = 0.0;
//    pr_gapy_full_path = 0.0;
//    for (int catg = 0; catg < numCatg; catg++) {
//
//        // PrfvL = Pr_L * fv_L
//        bpp::ColMatrix<double> PrfvL;
//        bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL.at(catg), PrfvL);
//
//        // PrfvR = Pr_R * fv_R
//        bpp::ColMatrix<double> PrfvR;
//        bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR.at(catg), PrfvR);
//
//        // fv = PrfvL * PrfvR
//        bpp::ColMatrix<double> fv;
//        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);
//
//        Fv_Y_ij[catg] = fv;
//
//        // fv0 = pi * fv
//        double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);
//
//        Fv_sigma_Y_ij.at(catg) = fv0;
//
//        // gapY probability with gamma
//        double p0 = progressivePIP_->rDist_->getProbability((size_t) catg) * \
//                    iotasNode_[catg] * betasNode_[catg] * fv0;
//
//        pr_gapy_full_path += progressivePIP_->rDist_->getProbability((size_t) catg) * \
//                              alphaNode_.at(catg) * fv0;
//
//        pr_gapy += p0;
//
//    }
//
//}

double nodeRAM::_computeLK_MXY(double log_phi_gamma,
                               double valM,
                               double valX,
                               double valY,
                               double log_pr) {

    // compute the lk at a given matrix entry extending the previous best
    // lk (valM,valX,valY) together with the actual lk value (log_pr) and
    // the marginal lk of an empty column

    return log_phi_gamma + log_pr + max_of_three(valM, valX, valY, DBL_EPSILON);
}

//std::vector<double> nodeRAM::_computeLK_empty(std::vector<bpp::ColMatrix<double> > &fvL,
//                                              std::vector<bpp::ColMatrix<double> > &fvR,
//                                              std::vector<bpp::ColMatrix<double> > &Fv_gap,
//                                              std::vector<double> &fv_empty_sigma_,
//                                              std::vector<double> &lk_empty_down_L,
//                                              std::vector<double> &lk_empty_down_R) {
//
//    // number of discrete gamma categories
//    int num_gamma_categories = progressivePIP_->numCatg_;
//
//    double fv0;
//    double p0;
//    double pL,pR;
//
//    // array of lk (for each gamma rate) of a single column full of gaps
//    std::vector<double> pc0;
//    pc0.resize(num_gamma_categories);
//
//    for (int catg = 0; catg < num_gamma_categories; catg++) {
//
//        // PrfvL = Pr_L * fv_L
//        bpp::ColMatrix<double> PrfvL;
//        bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL.at(catg), PrfvL);
//
//        // PrfvR = Pr_R * fv_R
//        bpp::ColMatrix<double> PrfvR;
//        bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR.at(catg), PrfvR);
//
//        // fv = PrfvL * PrfvR
//        bpp::ColMatrix<double> fv;
//        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);
//
//        Fv_gap.at(catg) = fv;
//
//        // fv0 = pi * fv
//        fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);
//
//        fv_empty_sigma_.at(catg) = fv0;
//
//        if(_isRootNode()){
//            // lk at root node (beta = 1.0)
//            p0 = iotasNode_.at(catg) * fv0;
//        }else{
//            p0 = ( iotasNode_.at(catg) - \
//                   iotasNode_.at(catg) * betasNode_.at(catg) + \
//                   iotasNode_.at(catg) * betasNode_.at(catg) * fv0 );
//        }
//
//        pL = lk_empty_down_L.at(catg);
//        pR = lk_empty_down_R.at(catg);
//
//        // TODO: remove this use instead lk_empty_down_
//        pc0.at(catg) = p0 + pL + pR;
//
//        lk_empty_.at(catg) = p0 + pL + pR;
//    }
//
//    return pc0;
//}

void nodeRAM::_computeLkEmptyLeaf(){

    // compute the lk of an empty column at the leaf

    // get the number of gamma categories
    int numCatg = progressivePIP_->numCatg_;

    // allocate memory ([numCatg] x 1)
    MSA_->pipmsa->lk_empty_.resize(numCatg);

    // only 1 column
    for (int catg=0; catg<numCatg; catg++) {
        // compute the lk of an empty column at the leaf
        MSA_->pipmsa->lk_empty_.at(catg) = progressivePIP_->rDist_->getProbability((size_t) catg) * \
            iotasNode_.at(catg) * (1 - betasNode_.at(catg));
    }

}

std::vector<double> nodeRAM::_computeLkEmptyNode(){

    // number of discrete gamma categories
    int numCatg = progressivePIP_->numCatg_;

    double p0;
    double pL,pR;

    // resize array ([numCatg] x 1)
    MSA_->pipmsa->lk_empty_.resize(numCatg);

    // array of lk (for each gamma rate) of a single column full of gaps
    std::vector<double> pc0;
    pc0.resize(numCatg);

    for (int catg = 0; catg < numCatg; catg++) {

        //double fv0 = fv_empty_sigma_.at(catg);

        double pr_up = 0.0; // lk_empty UP (from the actual node to the root)

        if(_isRootNode()){ // root
            // lk at root node (beta = 1.0)
            p0 = iotasNode_.at(catg) * MSA_->pipmsa->fv_empty_sigma_.at(catg);
        }else{ // internal node
            p0 = ( iotasNode_.at(catg) - \
                   iotasNode_.at(catg) * betasNode_.at(catg) + \
                   iotasNode_.at(catg) * betasNode_.at(catg) * MSA_->pipmsa->fv_empty_sigma_.at(catg) );

            // climb the tree and compute the probability UP
            bpp::PIPnode *tmpNode = this->parent;
            while(tmpNode){
                pr_up += tmpNode->iotasNode_.at(catg) - \
                         tmpNode->iotasNode_.at(catg) * tmpNode->betasNode_.at(catg) + \
                         tmpNode->iotasNode_.at(catg) * tmpNode->betasNode_.at(catg) * MSA_->pipmsa->fv_empty_sigma_.at(catg);
                tmpNode = tmpNode->parent;
            }

        }

        pL = childL->MSA_->pipmsa->lk_empty_.at(catg); // lk_empty DOWN left
        pR = childR->MSA_->pipmsa->lk_empty_.at(catg); // lk_empty DOWN right

        pc0.at(catg) = pr_up + p0 + pL + pR; // this lk_empty is used at this layer

        MSA_->pipmsa->lk_empty_.at(catg) = p0 + pL + pR; // here store the lk for the next layer (probability UP is not added here)
    }

    return pc0;

}



//void nodeRAM::_compute_lk_empty_down_rec(std::vector<double> &lk){
//    int num_gamma_categories = rDist_->getNumberOfCategories();
//
//    for (int catg=0; catg<num_gamma_categories; catg++) {
//        lk.at(catg) = lk.at(catg) + rDist_->getProbability((size_t) catg) * \
//            (iotasNode_.at(nodeID).at(catg) - \
//            iotasNode_.at(nodeID).at(catg) * betasNode_.at(nodeID).at(catg) + \
//            iotasNode_.at(nodeID).at(catg) * betasNode_.at(nodeID).at(catg) * \
//            fv_empty_sigma_.at(nodeID).at(catg));
//    }
//
//    if (!node->isLeaf()) {
//
//        tshlib::VirtualNode *vnode_left = treemap_.left.at(nodeID)->getNodeLeft();
//        int sonLeftID = treemap_.right.at(vnode_left);
//        bpp::Node *sonLeft = tree_->getNode(sonLeftID);
//
//        tshlib::VirtualNode *vnode_right = treemap_.left.at(nodeID)->getNodeRight();
//        int sonRightID = treemap_.right.at(vnode_right);
//        bpp::Node *sonRight = tree_->getNode(sonRightID);
//
//        _compute_lk_empty_down_rec(sonLeft,lk);
//        _compute_lk_empty_down_rec(sonRight,lk);
//
//    }
//}

//void nodeRAM::_compute_lk_empty_leaf_() {
//
//    // get the number of gamma categories
//    int numCatg = progressivePIP_->numCatg_;
//
//    // allocate memory ([numCatg] x 1)
//    lk_empty_.resize(numCatg);
//
//    // only 1 column
//    for (int catg=0; catg<numCatg; catg++) {
//        // compute the lk of an empty column at the leaf
//        lk_empty_.at(catg) = progressivePIP_->rDist_->getProbability((size_t) catg) * \
//            (iotasNode_.at(catg) - \
//            iotasNode_.at(catg) * betasNode_.at(catg) + \
//            iotasNode_.at(catg) * betasNode_.at(catg) * \
//            fv_empty_sigma_.at(catg));
//    }
//
//}

//void nodeRAM::_computeLkEmptyLeaf(){
//
//    // compute the lk of an empty column at the leaf
//
//    // get the number of gamma categories
//    int numCatg = progressivePIP_->numCatg_;
//
//    // allocate memory ([numCatg] x 1)
//    static_cast<PIPmsaSingle *>(MSA_)->pipmsa->lk_empty_.resize(numCatg);
//
//    // only 1 column
//    for (int catg=0; catg<numCatg; catg++) {
//        // compute the lk of an empty column at the leaf
//        static_cast<PIPmsaSingle *>(MSA_)->pipmsa->lk_empty_.at(catg) = progressivePIP_->rDist_->getProbability((size_t) catg) * \
//            iotasNode_.at(catg) * (1 - betasNode_.at(catg));
//    }
//
//}

void nodeRAM::_computeLkLeaf(){

    // compute the lk at the leaf

    // get the number of gamma categories
    int numCatg = progressivePIP_->numCatg_;

    // get the size of the compressed sequences
    int msaLen = MSA_->getCompressedMSAlength();

    // allocate memory ([site])
    MSA_->pipmsa->log_lk_down_.resize(msaLen);

    // compute the marginal lk over all the gamma categories
    for(int site=0;site<msaLen;site++){
        // init to 0.0
        MSA_->pipmsa->log_lk_down_.at(site) = 0.0;
        for (int catg=0; catg<numCatg; catg++) {
            MSA_->pipmsa->log_lk_down_.at(site) += progressivePIP_->rDist_->getProbability((size_t) catg) * \
                                     iotasNode_.at(catg) * betasNode_.at(catg) * \
                                    MSA_->pipmsa->fv_sigma_.at(site).at(catg);
        }
        // compute the log lk
        MSA_->pipmsa->log_lk_down_.at(site) = log(MSA_->pipmsa->log_lk_down_.at(site));
    }

}

void nodeRAM::_compressLK(std::vector<double> &lk_down_not_compressed){

    // compress an array of lk values

    int comprMSAlen = MSA_->getCompressedMSAlength();

    int id_map;

    MSA_->pipmsa->log_lk_down_.resize(comprMSAlen);

    for(int i=0;i<comprMSAlen;i++){
        id_map = MSA_->pipmsa->rev_map_compressed_seqs_.at(i);
        MSA_->pipmsa->log_lk_down_.at(i)=lk_down_not_compressed.at(id_map);
    }

}

void nodeRAM::_computeAllFvEmptySigmaRec() {

    if(childL == nullptr && childR == nullptr){ // leaf

        MSA_->pipmsa->_setFVemptyLeaf(progressivePIP_->numCatg_,
                                                                   progressivePIP_->alphabet_); // set fv_empty

        MSA_->pipmsa->_setFVsigmaEmptyLeaf(progressivePIP_->numCatg_); // set fv_sigma_empty = fv_empty dot pi

    }else{ // internal node

        // recursive call
        childL->_computeAllFvEmptySigmaRec();
        childR->_computeAllFvEmptySigmaRec();

        MSA_->_setFVemptyNode(progressivePIP_->numCatg_,
                              childL->MSA_->pipmsa,
                              childR->MSA_->pipmsa,
                              childL->prNode_,
                              childR->prNode_); // set fv_empty

        MSA_->_setFVsigmaEmptyNode(progressivePIP_->numCatg_,
                                   childL->MSA_->pipmsa,
                                   childR->MSA_->pipmsa,
                                   childL->bnode_->getDistanceToFather(),
                                   childR->bnode_->getDistanceToFather(),
                                   progressivePIP_->mu_); // set fv_sigma_empty = fv_empty dot pi
    }

}

void nodeRAM::DP3D_PIP_leaf() {

    //*******************************************************************************
    // ALIGNS LEAVES
    //*******************************************************************************

    // get vnode Id
    int vnodeId = (int) vnode_->vnode_seqid;

    // get sequence name from vnodeId
    std::string seqname = progressivePIP_->sequences_->getSequencesNames().at(vnodeId);

    // create a PIPmsaSingle object
    MSA_  = new PIPmsaSingle();

    // create a new PIPmsa
    MSA_->pipmsa = new PIPmsa();

    // associates the sequence name to the leaf node
    MSA_->pipmsa->_setSeqNameLeaf(seqname);

    // get sequence from sequence name
    const bpp::Sequence *sequence = &progressivePIP_->sequences_->getSequence(seqname);

    // creates a column containing the sequence associated to the leaf node
    MSA_->pipmsa->_setMSAleaf(sequence);

    // compresses sequence at the leaves
    MSA_->_compressMSA(progressivePIP_->alphabet_);

    // computes the indicator values (fv values) at the leaves
    MSA_->pipmsa->_setFVleaf(progressivePIP_->numCatg_,progressivePIP_->alphabet_);

    // computes dotprod(pi,fv)
    MSA_->pipmsa->_setFVsigmaLeaf(MSA_->getCompressedMSAlength(),
                                  progressivePIP_->numCatg_,
                                  progressivePIP_->pi_);

    // compute the lk of an empty column
    _computeLkEmptyLeaf();

    // computes the lk for all the characters at the leaf
    _computeLkLeaf();

    // sets the traceback path at the leaf
    MSA_->pipmsa->_setTracebackPathleaves();

}

void nodeRAM::DP3D_PIP_node() {

    //*******************************************************************************
    // ALIGNS INTERNAL NODE
    //*******************************************************************************

    DVLOG(1) << "DP3D_PIP at node: "<<bnode_->getName();

    //***************************************************************************************
    // DP VARIABLES
    //***************************************************************************************
    int m_binary_this; // Level Index during computation / current
    int m_binary_prev; // Level Index during computation / old
    auto epsilon = DBL_EPSILON; // very small number
    double min_inf = -std::numeric_limits<double>::infinity(); // -inf
    //***************************************************************************************
    // TRACEBACK VARIABLES
    //***************************************************************************************
    double curr_best_score = min_inf; // best likelihood value at this node
    double prev_best_score = min_inf; // previuous best value at this node
    int level_max_lk = INT_MIN; // depth in M,X,Y with the highest lk value
    int tr_index = (int)STOP_STATE; // traceback index: 1=MATCH, 2=GAPX, 3=GAPY
    double max_lk_val = min_inf; // best lk value
    //***************************************************************************************
    // RANDOM NUMBERS GENERATOR
    //***************************************************************************************
    std::default_random_engine generator(progressivePIP_->getSeed()); // jatiapp seed
    std::uniform_real_distribution<double> distribution(0.0,1.0); // uniform distribution for the selection
                                                                  // of lks with the same value
    //***************************************************************************************
    // EARLY STOP VARIABLES
    //***************************************************************************************
    bool flag_exit = false; // early stop condition flag
    int counter_to_early_stop = 0; // current number of consecutive steps where the lk decreases
    int max_decrease_before_stop = EARLY_STOP_THR; // hardcoded to prevent early-stops
    //***************************************************************************************
    // GAMMA VARIABLES
    //***************************************************************************************
    // number of discrete gamma categories
    size_t numCatg = progressivePIP_->numCatg_;
    //***************************************************************************************
    // GET SONS
    //***************************************************************************************
    bpp::Node *sonLeft = childL->_getBnode(); // bnode of the left child
    bpp::Node *sonRight = childR->_getBnode(); // bnode of the right child

    std::vector<int> *map_compr_L = &(childL->MSA_->pipmsa->map_compressed_seqs_);
    std::vector<int> *map_compr_R = &(childR->MSA_->pipmsa->map_compressed_seqs_);

    //======= DEBUG ==================================================================================================//
//    std::vector<int> *rev_map_compr_L = &(static_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->rev_map_compressed_seqs_);
//    std::vector<int> *rev_map_compr_R = &(static_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->rev_map_compressed_seqs_);
//
//    std::cout<<"\n\n";
//    std::cout<<"map_compr_L\n";
//    for(int i=0;i<map_compr_L->size();i++){
//        std::cout<<map_compr_L->at(i)<<" , ";
//    }
//    std::cout<<"\n\n";
//    std::cout<<"map_compr_R\n";
//    for(int i=0;i<map_compr_R->size();i++){
//        std::cout<<map_compr_R->at(i)<<" , ";
//    }
//    std::cout<<"\n\n";
//    std::cout<<"rev_map_compr_L\n";
//    for(int i=0;i<rev_map_compr_L->size();i++){
//        std::cout<<rev_map_compr_L->at(i)<<" , ";
//    }
//    std::cout<<"\n\n";
//    std::cout<<"rev_map_compr_R\n";
//    for(int i=0;i<rev_map_compr_R->size();i++){
//        std::cout<<rev_map_compr_R->at(i)<<" , ";
//    }
//    std::cout<<"\n\n";
    //======= DEBUG ==================================================================================================//

    //***************************************************************************************
    MSA_ = new PIPmsaSingle(); // object that store the MSA
    //***************************************************************************************
    // DP SIZES
    //***************************************************************************************
    // Compute dimensions of the 3D block at current internal node.
    int h = childL->MSA_->getMSAlength() + 1; // dimension of the alignment on the left side
    int w = childR->MSA_->getMSAlength() + 1; // dimension of the alignment on the right side
    int d = (h - 1) + (w - 1) + 1; // third dimension of the DP matrix
    int h_compr = childL->MSA_->getCompressedMSAlength(); // dimension of the compressed alignment on the left side
    int w_compr = childR->MSA_->getCompressedMSAlength(); // dimension of the compressed alignment on the right side
    //***************************************************************************************
    // WORKING VARIABLES
    //***************************************************************************************
    int i,j,k;
    int id1m,id2m;
    int id1x,id2y;
    //***************************************************************************************
    // MEMORY ALLOCATION
    //***************************************************************************************
    // Initialisation of the data structure
    std::vector< vector< vector<double> > > Log3DM;   // DP sparse matrix for MATCH case (only 2 layer are needed)
    std::vector< vector< vector<double> > > Log3DX;   // DP sparse matrix for GAPX case (only 2 layer are needed)
    std::vector< vector< vector<double> > > Log3DY;   // DP sparse matrix for GAPY case (only 2 layer are needed)
    std::vector< vector< vector<int> > > TR;          // 3D traceback matrix
    std::vector< vector<double> > Log2DM; // 2D matrix for MATCH case
    std::vector<double> Log2DX; // 1D matrix for GAPX case
    std::vector<double> Log2DY; // 1D matrix for GAPY case
    std::vector< vector<double> > Log2DM_fp; // 2D matrix for MATCH case
    std::vector<double> Log2DX_fp; // 1D matrix for GAPX case
    std::vector<double> Log2DY_fp; // 1D matrix for GAPY case
    std::vector< vector< vector< bpp::ColMatrix<double> > > > Fv_M;
    std::vector< vector< bpp::ColMatrix<double> > > Fv_X;
    std::vector< vector< bpp::ColMatrix<double> > > Fv_Y;
    std::vector< vector< vector<double> > > Fv_sigma_M;
    std::vector< vector<double> > Fv_sigma_X;
    std::vector< vector<double> > Fv_sigma_Y;
    //*************************
    Log3DM.resize(2);
    Log3DX.resize(2);
    Log3DY.resize(2);
    // allocate memory for the 2 layers
    for (k = 0; k < 2; k++){
        Log3DM[k].resize(h);
        Log3DX[k].resize(h);
        Log3DY[k].resize(h);
        for(i = 0; i < h; i++){
            Log3DM[k][i].resize(w,min_inf);
            Log3DX[k][i].resize(w,min_inf);
            Log3DY[k][i].resize(w,min_inf);
        }
    }
    //*************************
    TR.resize(d);
    TR[0].resize(1);
    TR[0][0].resize(1,STOP_STATE);

    for (i = 1; i < d; i++) {
        TR[i].resize(h);
        for(j = 0; j < h; j++){
            TR[i][j].resize(w,0);
        }
    }

    Log2DM.resize(h_compr);
    Log2DX.resize(h_compr);
    Log2DY.resize(w_compr);

    Log2DM_fp.resize(h_compr);
    Log2DX_fp.resize(h_compr);
    Log2DY_fp.resize(w_compr);

    Fv_M.resize(h_compr);
    Fv_X.resize(h_compr);
    Fv_Y.resize(w_compr);

    for(i = 0; i < h_compr; i++){
        Log2DM[i].resize(w_compr);
        Log2DM_fp[i].resize(w_compr);
        Fv_M[i].resize(w_compr);
        for(j = 0; j < w_compr; j++){
            Fv_M[i][j].resize(numCatg);
        }
    }
    for(i = 0; i < h_compr; i++){
        Fv_X[i].resize(numCatg);
    }
    for(j = 0; j < w_compr; j++){
        Fv_Y[j].resize(numCatg);
    }

    Fv_sigma_M.resize(h_compr);
    Fv_sigma_X.resize(h_compr);
    Fv_sigma_Y.resize(w_compr);

    for(i = 0; i < h_compr; i++){
        Fv_sigma_M[i].resize(w_compr);
        for(j = 0; j < w_compr; j++){
            Fv_sigma_M[i][j].resize(numCatg);
        }
    }
    for(i = 0; i < h_compr; i++){
        Fv_sigma_X[i].resize(numCatg);
    }
    for(j = 0; j < w_compr; j++){
        Fv_sigma_Y[j].resize(numCatg);
    }
    //***************************************************************************************
    // LK COMPUTATION OF AN EMPTY COLUMNS (FULL OF GAPS)
    //***************************************************************************************
    // computes the lk of an empty column in the two subtrees
    std::vector<double> pc0 = _computeLkEmptyNode();
    //***************************************************************************************
    // COMPUTES LOG(PHI(0))
    //***************************************************************************************
    // marginal likelihood for all empty columns with rate variation (gamma distribution)
    // phi(m,pc0,r) depends on the MSA length m
    // marginal phi marginalized over gamma categories
    double nu_gamma = 0.0;
    double log_phi_gamma = 0.0;
    for (int catg = 0; catg < numCatg; catg++) {
        // log( P_gamma(r) * phi(0,pc0(r),r) ): marginal lk for all empty columns of an alignment of size 0

        nu_gamma += progressivePIP_->rDist_->getProbability((size_t) catg) * \
                    progressivePIP_->nu_.at(catg);

        log_phi_gamma += progressivePIP_->rDist_->getProbability((size_t) catg) * \
                         (progressivePIP_->nu_.at(catg) * \
                         (pc0.at(catg)-1));
    }

    double log_nu_gamma = log(nu_gamma);
    //***************************************************************************************

    //***************************************************************************************
    Log3DM[0][0][0] = log_phi_gamma;
    Log3DX[0][0][0] = log_phi_gamma;
    Log3DY[0][0][0] = log_phi_gamma;
    TR[0][0][0] = STOP_STATE;
    //***************************************************************************************
    // 2D LK COMPUTATION
    //***************************************************************************************
    // computes the lk in the two subtrees
    std::vector<double> &lk_down_L = childL->MSA_->pipmsa->log_lk_down_;
    std::vector<double> &lk_down_R = childR->MSA_->pipmsa->log_lk_down_;

    //======= DEBUG ==================================================================================================//
//    std::cout<<"\nlk_down_L\n";
//    for(int i=0;i<lk_down_L.size();i++){
//        std::cout<<lk_down_L.at(i)<<";";
//    }
//    std::cout<<"\n\n";
//    std::cout<<"\nlk_down_R\n";
//    for(int i=0;i<lk_down_R.size();i++){
//        std::cout<<lk_down_R.at(i)<<";";
//    }
//    std::cout<<"\n\n";
    //======= DEBUG ==================================================================================================//

    // MATCH2D
    double pr_m;
    double pr_m_fp;
    for (i = 0; i < h_compr; i++) {
        for (j = 0; j < w_compr; j++) {

            _computeLK_M(childL->MSA_->pipmsa->fv_data_.at(i),
                         childR->MSA_->pipmsa->fv_data_.at(j),
                         Fv_M[i][j],
                         Fv_sigma_M[i][j],
                         pr_m,
                         pr_m_fp);

            Log2DM[i][j] = log(pr_m); // stored for the next layer
            Log2DM_fp[i][j] = log(pr_m_fp); // used at this node
        }
    }
    //***************************************************************************************
    // GAPX2D
    double pr_x;
    double pr_x_fp;
    for (i = 0; i < h_compr; i++) {

        _computeLK_X(childL->MSA_->pipmsa->fv_data_.at(i),
                     childR->MSA_->pipmsa->fv_empty_data_,
                     Fv_X[i],
                     Fv_sigma_X[i],
                     pr_x,
                     pr_x_fp);

        Log2DX[i] = progressivePIPutils::add_lns(log(pr_x),lk_down_L.at(i)); // stored for the next layer
        Log2DX_fp[i] = progressivePIPutils::add_lns(log(pr_x_fp),lk_down_L.at(i)); // used at this node
    }
    //***************************************************************************************
    // GAPY2D
    double pr_y;
    double pr_y_fp;
    for (j = 0; j < w_compr; j++) {

        _computeLK_Y(childL->MSA_->pipmsa->fv_empty_data_,
                     childR->MSA_->pipmsa->fv_data_.at(j),
                     Fv_Y[j],
                     Fv_sigma_Y[j],
                     pr_y,
                     pr_y_fp);

        Log2DY[j] = progressivePIPutils::add_lns(log(pr_y),lk_down_R.at(j)); // stored for the next layer
        Log2DY_fp[j] = progressivePIPutils::add_lns(log(pr_y_fp),lk_down_R.at(j)); // used at this node
    }
    //***************************************************************************************
    // 3D DYNAMIC PROGRAMMING
    //***************************************************************************************
    // For each slice of the 3D cube, compute the values of each cell
    for (int m = 1; m < d; m++) {

        // if lk doesn't increase anymore for K steps (EARLY_STOP_THR) exit
        if (flag_exit) {
            break;
        }

        //***********************************************************************************
        // alternate the two layers
        m_binary_this = m % 2;
        m_binary_prev = (m + 1) % 2;
        //***********************************************************************************

        //***********************************************************************************
        // delta phi to add
        log_phi_gamma = - log((long double) m) + log_nu_gamma;
        //***********************************************************************************

        //***************************************************************************
        // GAPX[i][0]
        j=0;
        for (i = 1; i < h; i++) {
            id1x = map_compr_L->at(i - 1);

            Log3DM[m_binary_this][i - 1][j] = min_inf;
            Log3DY[m_binary_this][i - 1][j] = min_inf;

            double val = _computeLK_MXY(log_phi_gamma,
                                        min_inf,
                                        Log3DX[m_binary_prev][i - 1][j],
                                        min_inf,
                                        Log2DX_fp[id1x]);

            Log3DX[m_binary_this][i][j] = val;

            TR[m][i][j] = (int)GAP_X_STATE;
        }
        //***********************************************************************************
        // GAPY[0][j]
        i=0;
        for (j = 1; j < w; j++) {
            id2y = map_compr_R->at(j - 1);

            Log3DM[m_binary_this][i][j - 1] = min_inf;
            Log3DX[m_binary_this][i][j - 1] = min_inf;
            Log3DY[m_binary_this][i][j] = _computeLK_MXY(log_phi_gamma,
                                                         min_inf,
                                                         min_inf,
                                                         Log3DY[m_binary_prev][i][j - 1],
                                                         Log2DY_fp[id2y]);

            TR[m][i][j] = (int)GAP_Y_STATE;
        }
        //***********************************************************************************
        for (i = 1; i < h; i++) {
            for (j = 1; j < w; j++) {
                //***************************************************************************
                // MATCH[i][j]
                id1m = map_compr_L->at(i - 1);
                id2m = map_compr_R->at(j - 1);

                Log3DM[m_binary_this][i][j] = _computeLK_MXY(log_phi_gamma,
                                                             Log3DM[m_binary_prev][i - 1][j - 1],
                                                             Log3DX[m_binary_prev][i - 1][j - 1],
                                                             Log3DY[m_binary_prev][i - 1][j - 1],
                                                             Log2DM_fp[id1m][id2m]);
                //***************************************************************************
                // GAPX[i][j]
                id1x = map_compr_L->at(i - 1);

                Log3DX[m_binary_this][i][j] = _computeLK_MXY(log_phi_gamma,
                                                             Log3DM[m_binary_prev][i - 1][j],
                                                             Log3DX[m_binary_prev][i - 1][j],
                                                             Log3DY[m_binary_prev][i - 1][j],
                                                             Log2DX_fp[id1x]);
                //***************************************************************************
                // GAPY[i][j]
                id2y = map_compr_R->at(j - 1);

                Log3DY[m_binary_this][i][j] = _computeLK_MXY(log_phi_gamma,
                                                             Log3DM[m_binary_prev][i][j - 1],
                                                             Log3DX[m_binary_prev][i][j - 1],
                                                             Log3DY[m_binary_prev][i][j - 1],
                                                             Log2DY_fp[id2y]);
                //***************************************************************************
                // TR[i][j]
                // Find which matrix contains the best value of LK found until this point.
                _index_of_max(Log3DM[m_binary_this][i][j],
                              Log3DX[m_binary_this][i][j],
                              Log3DY[m_binary_this][i][j],
                              epsilon,
                              generator,
                              distribution,
                              tr_index,
                              max_lk_val);

                // Store the index for the traceback
                TR[m][i][j] = tr_index;

                // If we reached the corner of the 3D cube, then:
                if ( (m>=(h-1)) && (m>=(w-1)) && (i == (h - 1)) && (j == (w - 1))  ) {
                    // the algorithm is filling the last column of 3D DP matrix where
                    // all the characters are in the MSA

                    if(tr_index==(int)STOP_STATE){
                        LOG(FATAL) <<"\nSomething went wrong in reading the TR value. "
                                     "TR is neither MATCH, nor GAPX, nor GAPY. ";
                    }

                    if (max_lk_val > curr_best_score) {
                        curr_best_score = max_lk_val;
                        level_max_lk = m;
                    }

                    //***********************************************************************
                    // early stop condition
                    if (curr_best_score < prev_best_score) {
                        prev_best_score = curr_best_score;
                        counter_to_early_stop++;
                        if (counter_to_early_stop > max_decrease_before_stop) {
                            // if for max_decrease_before_stop consecutive times
                            // the lk decrease then exit, the maximum lk has been reached
                            flag_exit = true;
                        }
                    } else {
                        counter_to_early_stop = 0;
                    }
                    //***************************************************************************

                }
            }
        }
    }
    //***************************************************************************************
    // STORE THE SCORE
    //***************************************************************************************
    // level (k position) in the DP matrix that contains the highest lk value

    // create a new iPIPmsa object of type PIPmsaSingle
    MSA_  = new PIPmsaSingle();

    // create a new PIPmsa
    MSA_->pipmsa = new PIPmsa();

    MSA_->pipmsa->score_ = curr_best_score;
    //***************************************************************************************
    // TRACEBACK ALGORITHM
    //***************************************************************************************
    std::vector< vector< bpp::ColMatrix<double> > > fv_data_not_compressed;
    fv_data_not_compressed.resize(level_max_lk);

    std::vector<std::vector<double>> fv_sigma_not_compressed;
    fv_sigma_not_compressed.resize(level_max_lk);

    std::vector<double> lk_down_not_compressed;
    lk_down_not_compressed.resize(level_max_lk);

    // start backtracing the 3 matrices (MATCH, GAPX, GAPY)
    MSA_->pipmsa->traceback_path_.resize(level_max_lk);

    i = h - 1;
    j = w - 1;
    int idmL,idmR;
    int state;
    for (int lev = level_max_lk; lev > 0; lev--) {
        state = TR[lev][i][j];
        switch (state) {
            case MATCH_STATE:

                idmL = map_compr_L->at(i-1);
                idmR = map_compr_R->at(j-1);

                fv_data_not_compressed.at(lev - 1) = Fv_M[idmL][idmR];
                fv_sigma_not_compressed.at(lev -1) = Fv_sigma_M[idmL][idmR];
                lk_down_not_compressed.at(lev - 1) = Log2DM[idmL][idmR];

                i = i - 1;
                j = j - 1;

                MSA_->pipmsa->traceback_path_.at(lev - 1) = (int)MATCH_STATE;

                break;
            case GAP_X_STATE:

                idmL = map_compr_L->at(i-1);

                fv_data_not_compressed.at(lev - 1) = Fv_X[idmL];
                fv_sigma_not_compressed.at(lev -1) = Fv_sigma_X[idmL];
                lk_down_not_compressed.at(lev - 1) = Log2DX[idmL];

                i = i - 1;

                MSA_->pipmsa->traceback_path_.at(lev - 1) = (int)GAP_X_STATE;

                break;
            case GAP_Y_STATE:

                idmR = map_compr_R->at(j-1);

                fv_data_not_compressed.at(lev - 1) = Fv_Y[idmR];
                fv_sigma_not_compressed.at(lev -1) = Fv_sigma_Y[idmR];
                lk_down_not_compressed.at(lev - 1) = Log2DY[idmR];

                j = j - 1;

                MSA_->pipmsa->traceback_path_.at(lev - 1) = (int)GAP_Y_STATE;

                break;
            default:
                LOG(FATAL) << "\nSomething went wrong during the alignment reconstruction in function "
                              "pPIP::DP3D_PIP. Check call stack below.";
        }
    }

    //======== DEBUG =======================================================================//
//    std::cout<<"\n\ntraceback_path\n";
//    for(int i=0;i<static_cast<PIPmsaSingle *>(MSA_)->pipmsa->traceback_path_.size();i++){
//
//        std::cout<<static_cast<PIPmsaSingle *>(MSA_)->pipmsa->traceback_path_.at(i)<<",";
//
//    }
//    std::cout<<"\n\n";
    //======== DEBUG =======================================================================//

    //***************************************************************************************
    // BUILD NEW MSA
    //***************************************************************************************
    // converts traceback path into an MSA
    PIPmsa *msaL = childL->MSA_->_getMSA();
    PIPmsa *msaR = childR->MSA_->_getMSA();
    MSA_->_build_MSA(msaL->msa_,msaR->msa_);

    // assigns the sequence names of the new alligned sequences to the current MSA
    std::vector<string> *seqNameL = &childL->MSA_->pipmsa->seqNames_;
    std::vector<string> *seqNameR = &childR->MSA_->pipmsa->seqNames_;
    MSA_->pipmsa->_setSeqNameNode(*seqNameL,*seqNameR);

    //***************************************************************************************
    // COMPRESS INFO
    //***************************************************************************************
    // compress the MSA
    MSA_->_compressMSA(progressivePIP_->alphabet_);

    // compress fv values and lk_down
    _compressLK(lk_down_not_compressed);
    MSA_->pipmsa->_compress_Fv(fv_sigma_not_compressed, fv_data_not_compressed);

}

void nodeRAM::DP3D_PIP() {

    if (_isTerminalNode()) {
        // align leaf (prepare data)
        DP3D_PIP_leaf();
    }else{
        // align internal node
        DP3D_PIP_node();
    }

}