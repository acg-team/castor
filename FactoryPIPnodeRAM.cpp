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

using namespace bpp;

void nodeRAM::DP3D_PIP() {

    if (bnode_->isLeaf()) {
        DP3D_PIP_leaf();
    }else{
        DP3D_PIP_node();
    }

}

bool nodeRAM::_index_of_max(double m,
                            double x,
                            double y,
                            double epsilon,
                            std::default_random_engine &generator,
                            std::uniform_real_distribution<double> &distribution,
                            int &index,
                            double &val) {

    double random_number;

    if (std::isinf(m) & std::isinf(x) & std::isinf(y)){
        index = int(STOP_STATE);
        val = -std::numeric_limits<double>::infinity();
        return true;
    }

    if (not(std::isinf(m)) & not(std::isinf(x)) & (fabs((m - x)) < epsilon)) {
        x = m;
    }

    if (not(std::isinf(m)) & not(std::isinf(y)) & (fabs((m - y)) < epsilon)) {
        y = m;
    }

    if (not(std::isinf(x)) & not(std::isinf(y)) & (fabs((x - y)) < epsilon)) {
        y = x;
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
                LOG(FATAL) << "\nSomething went wrong during the comparison in function pPIP::_index_of_max. Check call stack below.";
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
                LOG(FATAL) << "\nSomething went wrong during the comparison in function pPIP::_index_of_max. Check call stack below.";
                return false;
            }
        }
    } else {

        double mx = x;
        if (mx > y) {
            //m or x
            random_number = distribution(generator);
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
                LOG(FATAL) << "\nSomething went wrong during the comparison in function pPIP::_index_of_max. Check call stack below.";
                return false;
            }
        }
    }

}

double nodeRAM::max_of_three(double a, double b, double c, double epsilon) {

    if (fabs(a) < epsilon) {
        a = -std::numeric_limits<double>::infinity();
    }
    if (fabs(b) < epsilon) {
        b = -std::numeric_limits<double>::infinity();
    }
    if (fabs(c) < epsilon) {
        c = -std::numeric_limits<double>::infinity();
    }

    if (std::isinf(a) && std::isinf(b) && std::isinf(c)) {
        return -std::numeric_limits<double>::infinity();
    }

    if (a > b) {
        if (a > c) {
            return a;
        }
        return c;
    } else {
        if (b > c) {
            return b;
        }
        return c;
    }

}

void nodeRAM::_compress_Fv(std::vector<std::vector<double>> &fv_sigma_not_compressed,
                           std::vector<vector<bpp::ColMatrix<double> > > &fv_data_not_compressed){

    int comprMSAlen = static_cast<PIPmsaSingle *>(MSA_)->pipmsa->rev_map_compressed_seqs_.size();

    int id_map;

    fv_data_.resize(comprMSAlen);

    fv_sigma_.resize(comprMSAlen);

    for(int i=0;i<comprMSAlen;i++){
        id_map = static_cast<PIPmsaSingle *>(MSA_)->pipmsa->rev_map_compressed_seqs_.at(i);

        fv_data_.at(i)=fv_data_not_compressed.at(id_map);

        fv_sigma_.at(i)=fv_sigma_not_compressed.at(id_map);
    }

}

double nodeRAM::computeLK_M_local(std::vector< bpp::ColMatrix<double> > &fvL,
                                  std::vector< bpp::ColMatrix<double> > &fvR,
                                  std::vector< bpp::ColMatrix<double> > &Fv_M_ij,
                                  std::vector<double> &Fv_sigma_M_ij) {

    // number of discrete gamma categories
    int num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

    double pr = 0.0;
    for (int catg = 0; catg < num_gamma_categories; catg++) {

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL.at(catg), PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR.at(catg), PrfvR);

        // fv = PrfvL * PrfvR
        bpp::ColMatrix<double> fv;
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

        Fv_M_ij[catg] = fv;

        // fv0 = pi * fv
        double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

        Fv_sigma_M_ij.at(catg) = fv0;

        // match probability with gamma
        double p = progressivePIP_->rDist_->getProbability((size_t) catg) * \
               iotasNode_[catg] * \
               betasNode_[catg] * \
               fv0;

        pr += p;
    }

    return pr;
}

double nodeRAM::computeLK_X_local(std::vector< bpp::ColMatrix<double> > &fvL,
                                  std::vector< bpp::ColMatrix<double> > &fvR,
                                  std::vector< bpp::ColMatrix<double> > &Fv_X_ij,
                                  std::vector<double> &Fv_sigma_X_ij ) {

    // number of discrete gamma categories
    int num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

    double pr = 0.0;
    for (int catg = 0; catg < num_gamma_categories; catg++) {

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL.at(catg), PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR.at(catg), PrfvR);

        // fv = PrfvL * PrfvR
        bpp::ColMatrix<double> fv;
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

        Fv_X_ij[catg] = fv;

        // fv0 = pi * fv
        double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

        Fv_sigma_X_ij.at(catg) = fv0;

        // gapX probability with gamma
        double p0 = progressivePIP_->rDist_->getProbability((size_t) catg) * \
                childL->iotasNode_[catg] * \
                childR->betasNode_[catg] * \
                fv0;

        pr += p0;
    }

    return pr;

}

double nodeRAM::computeLK_Y_local(std::vector< bpp::ColMatrix<double> > &fvL,
                                  std::vector< bpp::ColMatrix<double> > &fvR,
                                  std::vector< bpp::ColMatrix<double> > &Fv_Y_ij,
                                  std::vector<double> &Fv_sigma_Y_ij) {

    // number of discrete gamma categories
    int num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

    double pr = 0.0;
    for (int catg = 0; catg < num_gamma_categories; catg++) {

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL.at(catg), PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR.at(catg), PrfvR);

        // fv = PrfvL * PrfvR
        bpp::ColMatrix<double> fv;
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

        Fv_Y_ij[catg] = fv;

        // fv0 = pi * fv
        double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

        Fv_sigma_Y_ij.at(catg) = fv0;

        // gapY probability with gamma
        double p0 = progressivePIP_->rDist_->getProbability((size_t) catg) * \
                    iotasNode_[catg] * \
                    betasNode_[catg] * \
                    fv0;

        pr += p0;

    }

    return pr;

}

double nodeRAM::computeLK_MXY_local(double log_phi_gamma,
                                    double valM,
                                    double valX,
                                    double valY,
                                    double log_pr) {

    return log_phi_gamma + log_pr + max_of_three(valM, valX, valY, DBL_EPSILON);
}

std::vector<double> nodeRAM::computeLK_GapColumn_local(std::vector< bpp::ColMatrix<double> > &fvL,
                                                       std::vector< bpp::ColMatrix<double> > &fvR,
                                                       std::vector< bpp::ColMatrix<double> > &Fv_gap,
                                                       std::vector<double> &fv_empty_sigma_,
                                                       std::vector<double> &lk_empty_down_L,
                                                       std::vector<double> &lk_empty_down_R) {

    // number of discrete gamma categories
    int num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

    double fv0;
    double p0;
    double pL,pR;

    // array of lk (for each gamma rate) of a single column full of gaps
    std::vector<double> pc0;
    pc0.resize(num_gamma_categories);

    for (int catg = 0; catg < num_gamma_categories; catg++) {

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL.at(catg), PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR.at(catg), PrfvR);

        // fv = PrfvL * PrfvR
        bpp::ColMatrix<double> fv;
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

        Fv_gap.at(catg) = fv;

        // fv0 = pi * fv
        fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

        fv_empty_sigma_.at(catg) = fv0;

        // lk at the actual node (considered as root node => beta = 1.0)
        p0 = iotasNode_[catg] * fv0;

        pL = lk_empty_down_L.at(catg);
        pR = lk_empty_down_R.at(catg);

        pc0.at(catg) = p0 + pL + pR;
    }

    return pc0;
}

void nodeRAM::DP3D_PIP_node() {

    DVLOG(1) << "DP3D_PIP at node: "<<bnode_->getName();

    //***************************************************************************************
    // DP VARIABLES
    //***************************************************************************************
    int m_binary_this; // Level Index during computation / current
    int m_binary_prev; // Level Index during computation / old
    auto epsilon = DBL_EPSILON;
    double min_inf = -std::numeric_limits<double>::infinity();
    //***************************************************************************************
    // TRACEBACK VARIABLES
    //***************************************************************************************
    double curr_best_score = min_inf; // best likelihood value at this node
    double prev_best_score = min_inf; // previuous best value at this node
    int level_max_lk = INT_MIN; // Depth in M,X,Y with the highest lk value
    int tr_index = (int)STOP_STATE; // traceback index: 1=MATCH, 2=GAPX, 3=GAPY
    double max_lk_val = min_inf;
    //***************************************************************************************
    // RANDOM NUMBERS GENERATOR
    //***************************************************************************************
    std::default_random_engine generator(progressivePIP_->getSeed()); // jatiapp seed
    std::uniform_real_distribution<double> distribution(0.0,1.0); // Uniform distribution for the selection of lks with the same value
    //***************************************************************************************
    // EARLY STOP VARIABLES
    //***************************************************************************************
    bool flag_exit = false; // early stop condition flag
    int counter_to_early_stop = 0; // current number of consecutive steps where the lk decreases
    int max_decrease_before_stop = 10; // hardcoded to prevent early-stops
    //***************************************************************************************
    // GAMMA VARIABLES
    //***************************************************************************************
    // number of discrete gamma categories
    size_t num_gamma_categories = progressivePIP_->numCatg_;
    //***************************************************************************************
    // GET SONS
    //***************************************************************************************
    bpp::Node *sonLeft = childL->_getBnode();
    bpp::Node *sonRight = childR->_getBnode();
    //***************************************************************************************
    MSA_ = new PIPmsaSingle();
    //***************************************************************************************
    // DP SIZES
    //***************************************************************************************
    // Compute dimensions of the 3D block at current internal node.
    int h = static_cast<PIPmsaSingle *>(childL->MSA_)->getMSAlength() + 1; // dimension of the alignment on the left side
    int w = static_cast<PIPmsaSingle *>(childR->MSA_)->getMSAlength() + 1; // dimension of the alignment on the right side
    int d = (h - 1) + (w - 1) + 1; // third dimension of the DP matrix
    int h_compr = static_cast<PIPmsaSingle *>(childL->MSA_)->getCompressedMSAlength(); //= rev_map_compressed_seqs_.at(nodeID_L).at(0).size(); // dimension of the compressed alignment on the left side
    int w_compr = static_cast<PIPmsaSingle *>(childR->MSA_)->getCompressedMSAlength(); // dimension of the compressed alignment on the right side
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
    std::vector< vector< vector<int> > > TR;        // 3D traceback matrix
    std::vector< vector<double> > Log2DM;
    std::vector<double> Log2DX;
    std::vector<double> Log2DY;
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
        //TODO: pre-allocate only half of the enire matrix
        TR[i].resize(h);
        for(j = 0; j < h; j++){
            TR[i][j].resize(w,0);
        }
    }

    Log2DM.resize(h_compr);
    Log2DX.resize(h_compr);
    Log2DY.resize(w_compr);

    Fv_M.resize(h_compr);
    Fv_X.resize(h_compr);
    Fv_Y.resize(w_compr);

    for(i = 0; i < h_compr; i++){
        Log2DM[i].resize(w_compr);
        Fv_M[i].resize(w_compr);
        for(j = 0; j < w_compr; j++){
            Fv_M[i][j].resize(num_gamma_categories);
        }
    }
    for(i = 0; i < h_compr; i++){
        Fv_X[i].resize(num_gamma_categories);
    }
    for(j = 0; j < w_compr; j++){
        Fv_Y[j].resize(num_gamma_categories);
    }

    Fv_sigma_M.resize(h_compr);
    Fv_sigma_X.resize(h_compr);
    Fv_sigma_Y.resize(w_compr);

    for(i = 0; i < h_compr; i++){
        Fv_sigma_M[i].resize(w_compr);
        for(j = 0; j < w_compr; j++){
            Fv_sigma_M[i][j].resize(num_gamma_categories);
        }
    }
    for(i = 0; i < h_compr; i++){
        Fv_sigma_X[i].resize(num_gamma_categories);
    }
    for(j = 0; j < w_compr; j++){
        Fv_sigma_Y[j].resize(num_gamma_categories);
    }
    //***************************************************************************************
    // LK COMPUTATION OF AN EMPTY COLUMNS (FULL OF GAPS)
    //***************************************************************************************
    // computes the lk of an empty column in the two subtrees
    std::vector<double> lk_empty_down_L;// = _compute_lk_empty_down(sonLeft);
    std::vector<double> lk_empty_down_R;// = _compute_lk_empty_down(sonRight);

    // lk of a single empty column (full of gaps) with rate variation (gamma distribution)
    // compute the lk of a column full of gaps
    std::vector<double> pc0 = computeLK_GapColumn_local(childL->fv_empty_data_,
                                                        childR->fv_empty_data_,
                                                        fv_empty_data_,
                                                        fv_empty_sigma_,
                                                        lk_empty_down_L,
                                                        lk_empty_down_R);
    //***************************************************************************************
    // COMPUTES LOG(PHI(0))
    //***************************************************************************************
    // marginal likelihood for all empty columns with rate variation (gamma distribution)
    // phi(m,pc0,r) depends on the MSA length m
    // marginal phi marginalized over gamma categories
    double nu_gamma = 0.0;
    double log_phi_gamma = 0.0;
    for (int catg = 0; catg < num_gamma_categories; catg++) {
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
    //TODO: usare lk memorizzato
    std::vector<double> lk_down_L;// = _compute_lk_down(sonLeft,0);
    std::vector<double> lk_down_R;// = _compute_lk_down(sonRight,0);

    // MATCH2D
    for (i = 0; i < h_compr; i++) {
        for (j = 0; j < w_compr; j++) {
            Log2DM[i][j] = log(computeLK_M_local(childL->fv_data_.at(i),
                                                 childR->fv_data_.at(j),
                                                 Fv_M[i][j],
                                                 Fv_sigma_M[i][j]));
        }
    }
    //***************************************************************************************
    // GAPX2D
    for (i = 0; i < h_compr; i++) {
        Log2DX[i] = log(computeLK_X_local(childL->fv_data_.at(i),
                                          childR->fv_empty_data_,
                                          Fv_X[i],
                                          Fv_sigma_X[i]) + \
                                          lk_down_L.at(i));
    }
    //***************************************************************************************
    // GAPY2D
    for (j = 0; j < w_compr; j++) {
        Log2DY[j] = log(computeLK_Y_local(childL->fv_empty_data_,
                                          childR->fv_data_.at(j),
                                          Fv_Y[j],
                                          Fv_sigma_Y[j]) + \
                                          lk_down_R.at(j));
    }
    //***************************************************************************************
    // 3D DYNAMIC PROGRAMMING
    //***************************************************************************************
    // For each slice of the 3D cube, compute the values of each cell
    for (int m = 1; m < d; m++) {

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
            id1x = static_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->map_compressed_seqs_.at(i - 1);

            Log3DM[m_binary_this][i - 1][j] = min_inf;
            Log3DY[m_binary_this][i - 1][j] = min_inf;

            double val = computeLK_MXY_local(log_phi_gamma,
                                             min_inf,
                                             Log3DX[m_binary_prev][i - 1][j],
                                             min_inf,
                                             Log2DX[id1x]);

            Log3DX[m_binary_this][i][j] = val;

            _index_of_max(min_inf,
                          Log3DX[m_binary_this][i][j],
                          min_inf,
                          epsilon,
                          generator,
                          distribution,
                          tr_index,
                          max_lk_val);

            TR[m][i][j] = tr_index;
        }
        //***********************************************************************************
        // GAPY[0][j]
        i=0;
        for (j = 1; j < w; j++) {
            id2y = static_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->map_compressed_seqs_.at(j - 1);

            Log3DM[m_binary_this][i][j - 1] = min_inf;
            Log3DX[m_binary_this][i][j - 1] = min_inf;
            Log3DY[m_binary_this][i][j] = computeLK_MXY_local(log_phi_gamma,
                                                              min_inf,
                                                              min_inf,
                                                              Log3DY[m_binary_prev][i][j - 1],
                                                              Log2DY[id2y]);

            _index_of_max(min_inf,
                          min_inf,
                          Log3DY[m_binary_this][i][j],
                          epsilon,
                          generator,
                          distribution,
                          tr_index,
                          max_lk_val);

            TR[m][i][j] = tr_index;
        }
        //***********************************************************************************
        for (i = 1; i < h; i++) {
            for (j = 1; j < w; j++) {
                //***************************************************************************
                // MATCH[i][j]
                id1m = static_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->map_compressed_seqs_.at(i - 1);
                id2m = static_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->map_compressed_seqs_.at(j - 1);

                Log3DM[m_binary_this][i][j] = computeLK_MXY_local(log_phi_gamma,
                                                                  Log3DM[m_binary_prev][i-1][j-1],
                                                                  Log3DX[m_binary_prev][i-1][j-1],
                                                                  Log3DY[m_binary_prev][i-1][j-1],
                                                                  Log2DM[id1m][id2m]);
                //***************************************************************************
                // GAPX[i][j]
                id1x = static_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->map_compressed_seqs_.at(i - 1);

                Log3DX[m_binary_this][i][j] = computeLK_MXY_local(log_phi_gamma,
                                                                  Log3DM[m_binary_prev][i-1][j],
                                                                  Log3DX[m_binary_prev][i-1][j],
                                                                  Log3DY[m_binary_prev][i-1][j],
                                                                  Log2DX[id1x]);
                //***************************************************************************
                // GAPY[i][j]
                id2y = static_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->map_compressed_seqs_.at(j - 1);

                Log3DY[m_binary_this][i][j] = computeLK_MXY_local(log_phi_gamma,
                                                                  Log3DM[m_binary_prev][i][j-1],
                                                                  Log3DX[m_binary_prev][i][j-1],
                                                                  Log3DY[m_binary_prev][i][j-1],
                                                                  Log2DY[id2y]);
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
                TR[m][i][j] = tr_index; //max_val_index.index;

                // If we reached the corner of the 3D cube, then:
                if ( (m>=(h-1)) && (m>=(w-1)) && (i == (h - 1)) && (j == (w - 1))  ) {
                    // the algorithm is filling the last column of 3D DP matrix where
                    // all the characters are in the MSA

                    if(tr_index==(int)STOP_STATE){
                        LOG(FATAL) <<"\nSomething went wrong in reading the TR value. TR is neither MATCH, nor GAPX, nor GAPY. ";
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
    static_cast<PIPmsaSingle *>(MSA_)->pipmsa->score_ = curr_best_score;
    //***************************************************************************************
    // TRACEBACK ALGORITHM
    //***************************************************************************************
    std::vector< vector< bpp::ColMatrix<double> > > fv_data_not_compressed;
    fv_data_not_compressed.resize(level_max_lk);

    std::vector<std::vector<double>> fv_sigma_not_compressed;

    fv_sigma_not_compressed.resize(level_max_lk);

    // start backtracing the 3 matrices (MATCH, GAPX, GAPY)
    static_cast<PIPmsaSingle *>(MSA_)->pipmsa->traceback_path_.resize(level_max_lk);

    static_cast<PIPmsaSingle *>(MSA_)->pipmsa->traceback_mapL_.resize(level_max_lk);
    static_cast<PIPmsaSingle *>(MSA_)->pipmsa->traceback_mapR_.resize(level_max_lk);

    i = h - 1;
    j = w - 1;
    int idmL,idmR;
    int state;
    for (int lev = level_max_lk; lev > 0; lev--) {
        state = TR[lev][i][j];
        switch (state) {
            case MATCH_STATE:

                idmL = static_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->map_compressed_seqs_.at(i-1);
                idmR = static_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->map_compressed_seqs_.at(j-1);

                fv_data_not_compressed.at(lev - 1) = Fv_M[idmL][idmR];

                fv_sigma_not_compressed.at(lev -1) = Fv_sigma_M[idmL][idmR];

                i = i - 1;
                j = j - 1;

                break;
            case GAP_X_STATE:

                idmL = static_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->map_compressed_seqs_.at(i-1);

                fv_data_not_compressed.at(lev - 1) = Fv_X[idmL];

                fv_sigma_not_compressed.at(lev -1) = Fv_sigma_X[idmL];

                i = i - 1;

                break;
            case GAP_Y_STATE:

                idmR = static_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->map_compressed_seqs_.at(j-1);

                fv_data_not_compressed.at(lev - 1) = Fv_Y[idmR];

                fv_sigma_not_compressed.at(lev -1) = Fv_sigma_Y[idmR];

                j = j - 1;

                break;
            default:
                LOG(FATAL) << "\nSomething went wrong during the alignment reconstruction in function pPIP::DP3D_PIP. Check call stack below.";
        }
    }
    //***************************************************************************************
    // BUILD NEW MSA
    //***************************************************************************************
    // converts traceback path into an MSA
    PIPmsa *msaL = static_cast<PIPmsaSingle *>(childL->MSA_)->getMSA();
    PIPmsa *msaR = static_cast<PIPmsaSingle *>(childR->MSA_)->getMSA();
    static_cast<PIPmsaSingle *>(MSA_)->_build_MSA(msaL->msa_,msaR->msa_);

    // assigns the sequence names of the new alligned sequences to the current MSA
    std::vector<string> *seqNameL = &static_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->seqNames_;
    std::vector<string> *seqNameR = &static_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->seqNames_;
    static_cast<PIPmsaSingle *>(MSA_)->pipmsa->_setSeqNameNode(*seqNameL,*seqNameR);

    //***************************************************************************************
    // COMPRESS INFO
    //***************************************************************************************
    // compress the MSA
    static_cast<PIPmsaSingle *>(MSA_)->_compressMSA(progressivePIP_->alphabet_);

    // compress fv values and lk_down
    //_compress_lk_components(node, lk_down_not_compressed, fv_data_not_compressed);
    //static_cast<PIPmsaSingle *>(MSA_)->_compress_lk_components(fv_data_not_compressed)
    _compress_Fv(fv_sigma_not_compressed, fv_data_not_compressed);

}
