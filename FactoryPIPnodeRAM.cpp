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
 * @file FactoryPIPnodeRAM.cpp
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

void nodeRAM::_computeLkEmptyLeaf(){

    // compute the lk of an empty column at the leaf

    // get the number of gamma categories
    int numCatg = progressivePIP_->numCatg_;

    // allocate memory ([numCatg] x 1)
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->lk_empty_.resize(numCatg);

    // only 1 column
    for (int catg=0; catg<numCatg; catg++) {
        // compute the lk of an empty column at the leaf
        dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->lk_empty_.at(catg) = progressivePIP_->rDist_->getProbability((size_t) catg) * \
            iotasNode_.at(catg) * (1 - betasNode_.at(catg));
    }

}

void nodeRAM::_computeLkLeaf(){

    // compute the lk at the leaf

    // get the number of gamma categories
    int numCatg = progressivePIP_->numCatg_;

    // get the size of the compressed sequences
    int msaLen = dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_getCompressedMSAlength();

    // allocate memory ([site])
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->log_lk_down_.resize(msaLen);

    // compute the marginal lk over all the gamma categories
    for(int site=0;site<msaLen;site++){
        // init to 0.0
        dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->log_lk_down_.at(site) = 0.0;
        for (int catg=0; catg<numCatg; catg++) {
            dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->log_lk_down_.at(site) += progressivePIP_->rDist_->getProbability((size_t) catg) * \
                                     iotasNode_.at(catg) * betasNode_.at(catg) * \
                                    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->fv_sigma_.at(site).at(catg);
        }
        // compute the log lk
        dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->log_lk_down_.at(site) = log(dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->log_lk_down_.at(site));
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

    // associates the sequence name to the leaf node
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_setSeqNameLeaf(seqname);

    // get sequence from sequence name
    const bpp::Sequence *sequence = &progressivePIP_->sequences_->getSequence(seqname);

    // creates a vector containing the sequence associated to the leaf node
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_setMSAleaf(sequence);

    // compresses sequence at the leaves
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_compressMSA(progressivePIP_->alphabet_);

    // set fv_empty
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_setFVemptyLeaf(progressivePIP_->numCatg_,
                                                                progressivePIP_->alphabet_);

    // set fv_sigma_empty = fv_empty dot pi
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_setFVsigmaEmptyLeaf(progressivePIP_->numCatg_);

    // computes the indicator values (fv values) at the leaves
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_setFVleaf(progressivePIP_->numCatg_,
                             progressivePIP_->alphabet_);

    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_setFVsigmaLeaf(progressivePIP_->numCatg_,
                                                                progressivePIP_->pi_);

    // compute the lk of an empty column
    _computeLkEmptyLeaf();

    // computes the lk for all the characters at the leaf
    _computeLkLeaf();

    // sets the traceback path at the leaf
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_setTracebackPathLeaf();

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

    std::vector<int> *map_compr_L = &(dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->map_compressed_seqs_);
    std::vector<int> *map_compr_R = &(dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->map_compressed_seqs_);

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
    //MSA_ = new PIPmsaSingle(); // object that store the MSA
    //***************************************************************************************
    // DP SIZES
    //***************************************************************************************
    // Compute dimensions of the 3D block at current internal node.
    int h = dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->_getMSAlength() + 1; // dimension of the alignment on the left side
    int w = dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->_getMSAlength() + 1; // dimension of the alignment on the right side
    int d = (h - 1) + (w - 1) + 1; // third dimension of the DP matrix
    int h_compr = dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->_getCompressedMSAlength(); // dimension of the compressed alignment on the left side
    int w_compr = dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->_getCompressedMSAlength(); // dimension of the compressed alignment on the right side
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
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->lk_empty_.resize(numCatg);
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->fv_empty_data_.resize(numCatg);
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->fv_empty_sigma_.resize(numCatg);

    std::vector<bpp::ColMatrix<double> > &fvL = dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->fv_empty_data_;
    std::vector<bpp::ColMatrix<double> > &fvR = dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->fv_empty_data_;

    std::vector<bpp::ColMatrix<double> > &fv_empty_data = dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->fv_empty_data_;
    std::vector<double> &fv_empty_sigma = dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->fv_empty_sigma_;

    std::vector<double> &lk_emptyL = dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->lk_empty_;
    std::vector<double> &lk_emptyR = dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->lk_empty_;

    std::vector<double> &lk_empty = dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->lk_empty_;

    std::vector<double> pc0 = _computeLkEmptyNode(fvL,
                                                  fvR,
                                                  fv_empty_data,
                                                  fv_empty_sigma,
                                                  lk_emptyL,
                                                  lk_emptyR,
                                                  lk_empty);
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
    std::vector<double> &lk_down_L = dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->log_lk_down_;
    std::vector<double> &lk_down_R = dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->log_lk_down_;

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

            _computeLK_M(dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->fv_data_.at(i),
                         dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->fv_data_.at(j),
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

        _computeLK_X(dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->fv_data_.at(i),
                     dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->fv_empty_data_,
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

        _computeLK_Y(dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->fv_empty_data_,
                     dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->fv_data_.at(j),
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
    //MSA_  = new PIPmsaSingle();

    // create a new PIPmsa
    //MSA_->pipmsa = new PIPmsa();

    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->score_ = curr_best_score;
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
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->traceback_path_.resize(level_max_lk);

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

                dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->traceback_path_.at(lev - 1) = (int)MATCH_STATE;

                break;
            case GAP_X_STATE:

                idmL = map_compr_L->at(i-1);

                fv_data_not_compressed.at(lev - 1) = Fv_X[idmL];
                fv_sigma_not_compressed.at(lev -1) = Fv_sigma_X[idmL];
                lk_down_not_compressed.at(lev - 1) = Log2DX[idmL];

                i = i - 1;

                dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->traceback_path_.at(lev - 1) = (int)GAP_X_STATE;

                break;
            case GAP_Y_STATE:

                idmR = map_compr_R->at(j-1);

                fv_data_not_compressed.at(lev - 1) = Fv_Y[idmR];
                fv_sigma_not_compressed.at(lev -1) = Fv_sigma_Y[idmR];
                lk_down_not_compressed.at(lev - 1) = Log2DY[idmR];

                j = j - 1;

                dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->traceback_path_.at(lev - 1) = (int)GAP_Y_STATE;

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
    MSA_t *msaL = dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->_getMSA();
    MSA_t *msaR = dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->_getMSA();
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_build_MSA(*msaL,*msaR);

    // assigns the sequence names of the new alligned sequences to the current MSA
    std::vector<string> *seqNameL = &dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->seqNames_;
    std::vector<string> *seqNameR = &dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->seqNames_;
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_setSeqNameNode(*seqNameL,*seqNameR);

    //***************************************************************************************
    // COMPRESS INFO
    //***************************************************************************************
    // compress the MSA
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_compressMSA(progressivePIP_->alphabet_);

    // compress fv values and lk_down
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_compressLK(lk_down_not_compressed);
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_compressFv(fv_data_not_compressed);
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_compressFvSigma(fv_sigma_not_compressed);

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