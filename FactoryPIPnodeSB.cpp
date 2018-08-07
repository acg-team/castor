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
 * @file FactoryPIPnodeSB.cpp
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
#include <algorithm>
#include <vector>
#include "progressivePIP.hpp"
#include "FactoryPIPnode.hpp"

using namespace bpp;

void nodeSB::startingLevelSB(std::vector< vector< vector<double> > > &Log3DM,
                           std::vector< vector< vector<double> > > &Log3DX,
                           std::vector< vector< vector<double> > > &Log3DY,
                           double epsilon,
                           std::default_random_engine &generator,
                           std::uniform_real_distribution<double> &distribution,
                           int d,
                           int h,
                           int w,
                           int &lev,
                           double &val,
                           int &state){

    double sumM,sumX,sumY;
    // get sum of columns
    for(int k=0; k<d;k++){
        if(! std::isinf(Log3DM[k][h - 1][w - 1])){
            sumM += Log3DM[k][h - 1][w - 1];
        }
        if(! std::isinf(Log3DX[k][h - 1][w - 1])) {
            sumX += Log3DX[k][h - 1][w - 1];
        }
        if(! std::isinf(Log3DY[k][h - 1][w - 1])) {
            sumY += Log3DY[k][h - 1][w - 1];
        }
    }

    // get state and max_score
    _index_of_max(sumM,
                  sumX,
                  sumY,
                  epsilon,
                  generator,
                  distribution,
                  state,
                  val);

    std::vector< vector< vector<double> > > *mat3D;

    // get level
    double random_number = distribution(generator);
    switch (state) {
        case MATCH_STATE:
            random_number *= (-sumM);
            mat3D = &Log3DM;
            break;
        case GAP_X_STATE:
            random_number *= (-sumX);
            mat3D = &Log3DX;
            break;
        case GAP_Y_STATE:
            random_number *= (-sumY);
            mat3D = &Log3DY;
            break;
        default:
            LOG(FATAL) << "\nSomething went wrong in reading the STATE value."
                          " STATE is neither MATCH, nor GAPX, nor GAPY. ";
    }

    double cumsum=0.0;
    for(int level=0; level < d; level++){

        std::cout<<mat3D->at(level).at(h - 1).at(w - 1)<<"\n";

        if(! std::isinf(mat3D->at(level).at(h - 1).at(w - 1))){

            cumsum += abs(mat3D->at(level).at(h - 1).at(w - 1));

            if(cumsum > random_number){
                lev = level;
                break;
            }
        }
    }

}

void nodeSB::DP3D_PIP_leaf() {

    //*******************************************************************************
    // ALIGNS LEAVES
    //*******************************************************************************

    // get vnode Id
    int vnodeId = (int) vnode_->vnode_seqid;

    // get sequence name from vnodeId
    std::string seqname = progressivePIP_->sequences_->getSequencesNames().at(vnodeId);

    // associates the sequence name to the leaf node
    // leaves nodes have only 1 possible MSA which is the input sequence
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_setSeqNameLeaf(seqname);

    // get sequence from sequence name
    const bpp::Sequence *sequence = &progressivePIP_->sequences_->getSequence(seqname);

    // creates a vector containing the sequence associated to the leaf node
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_setMSAleaf(sequence);

    // compresses sequence at the leaves
    // compress MSA at position 0, the only one present at the leaf
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_compressMSA(progressivePIP_->alphabet_);

    // set fv_empty
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_setFVemptyLeaf(progressivePIP_->numCatg_,
                                                                progressivePIP_->alphabet_);

    // set fv_sigma_empty = fv_empty dot pi
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_setFVsigmaEmptyLeaf(progressivePIP_->numCatg_);

    // computes the indicator values (fv values) at the leaves
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_setFVleaf(progressivePIP_->numCatg_,
                                   progressivePIP_->alphabet_);

    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_setFVsigmaLeaf(progressivePIP_->numCatg_,
                                                                    progressivePIP_->pi_);

    // compute the lk of an empty column
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_computeLkEmptyLeaf(progressivePIP_,
                                                                        iotasNode_,
                                                                        betasNode_);

    // computes the lk for all the characters at the leaf
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_computeLkLeaf(progressivePIP_,
                                                                   iotasNode_,
                                                                   betasNode_);

    // sets the traceback path at the leaf
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_setTracebackPathLeaf();

    subMSAidxL_.resize(1);
    subMSAidxL_.at(0) = 0;

    subMSAidxR_.resize(1);
    subMSAidxR_.at(0) = 0;
}

void nodeSB::DP3D_PIP_node(int msa_idx_L,int msa_idx_R,int &position) {

    //TODO remove
    double temperature = 1.0;

    DVLOG(1) << "DP3D_PIP at node: "<<bnode_->getName();

    //***************************************************************************************
    // DP VARIABLES
    //***************************************************************************************
    double min_inf = -std::numeric_limits<double>::infinity(); // -inf
    //***************************************************************************************
    // RANDOM NUMBERS GENERATOR
    //***************************************************************************************
    std::default_random_engine generator(progressivePIP_->getSeed()); // jatiapp seed
    std::uniform_real_distribution<double> distribution(0.0,1.0); // uniform distribution for the selection
    // of lks with the same value
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
    std::vector<int> *map_compr_L = &(dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->map_compressed_seqs_);
    std::vector<int> *map_compr_R = &(dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->map_compressed_seqs_);
    //***************************************************************************************
    //MSA_ = new PIPmsaComp(); // object that store the MSA
    //***************************************************************************************
    // DP SIZES
    //***************************************************************************************
    // Compute dimensions of the 3D block at current internal node.
    int h = dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->_getMSAlength() + 1; // dimension of the alignment on the left side
    int w = dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->_getMSAlength() + 1; // dimension of the alignment on the right side
    int d = (h - 1) + (w - 1) + 1; // third dimension of the DP matrix
    int h_compr = dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->_getCompressedMSAlength(); // dimension of the compressed alignment on the left side
    int w_compr = dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->_getCompressedMSAlength(); // dimension of the compressed alignment on the right side
    //***************************************************************************************
    subMSAidxL_.at(position) = msa_idx_L;
    subMSAidxR_.at(position) = msa_idx_R;
    //***************************************************************************************
    // WORKING VARIABLES
    //***************************************************************************************
    int i,j,k;
    int id1m,id2m;
    int id1x,id2y;
    int m;
    //***************************************************************************************
    // MEMORY ALLOCATION
    //***************************************************************************************
    // Initialisation of the data structure
    std::vector< vector< vector<double> > > Log3DM;   // DP sparse matrix for MATCH case (only 2 layer are needed)
    std::vector< vector< vector<double> > > Log3DX;   // DP sparse matrix for GAPX case (only 2 layer are needed)
    std::vector< vector< vector<double> > > Log3DY;   // DP sparse matrix for GAPY case (only 2 layer are needed)
    std::vector< vector<double> > Log2DM;
    std::vector<double> Log2DX;
    std::vector<double> Log2DY;
    std::vector< vector<double> > Log2DM_fp;
    std::vector<double> Log2DX_fp;
    std::vector<double> Log2DY_fp;
    std::vector< vector< vector< bpp::ColMatrix<double> > > > Fv_M;
    std::vector< vector< bpp::ColMatrix<double> > > Fv_X;
    std::vector< vector< bpp::ColMatrix<double> > > Fv_Y;
    std::vector< vector< vector<double> > > Fv_sigma_M;
    std::vector< vector<double> > Fv_sigma_X;
    std::vector< vector<double> > Fv_sigma_Y;
    //*************************
    Log3DM.resize(d);
    Log3DX.resize(d);
    Log3DY.resize(d);
    // allocate memory for the 2 layers
    for (k = 0; k < d; k++){
        Log3DM[k].resize(h);
        Log3DX[k].resize(h);
        Log3DY[k].resize(h);
        for(int i = 0; i < h; i++){
            Log3DM[k][i].resize(w,min_inf);
            Log3DX[k][i].resize(w,min_inf);
            Log3DY[k][i].resize(w,min_inf);
        }
    }
    //*************************
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

    double nu_gamma;
    double log_phi_gamma;
    double log_nu_gamma;
    int local_position = position;
    for(int sb=0;sb<progressivePIP_->num_sb_;sb++) {

        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->lk_empty_.resize(numCatg);
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->fv_empty_data_.resize(numCatg);
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->fv_empty_sigma_.resize(numCatg);

        std::vector<bpp::ColMatrix<double> > &fvL = dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->fv_empty_data_;
        std::vector<bpp::ColMatrix<double> > &fvR = dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->fv_empty_data_;

        std::vector<bpp::ColMatrix<double> > &fv_empty_data = dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->fv_empty_data_;
        std::vector<double> &fv_empty_sigma = dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->fv_empty_sigma_;

        std::vector<double> &lk_emptyL = dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->lk_empty_;
        std::vector<double> &lk_emptyR = dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->lk_empty_;

        std::vector<double> &lk_empty = dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->lk_empty_;

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
        nu_gamma = 0.0;
        log_phi_gamma = 0.0;
        for (int catg = 0; catg < numCatg; catg++) {
            // log( P_gamma(r) * phi(0,pc0(r),r) ): marginal lk for all empty columns of an alignment of size 0
            nu_gamma += progressivePIP_->rDist_->getProbability((size_t) catg) * progressivePIP_->nu_.at(catg);
            log_phi_gamma += progressivePIP_->rDist_->getProbability((size_t) catg) * (progressivePIP_->nu_.at(catg) *\
                             (pc0.at(catg)-1));
        }

        log_nu_gamma = log(nu_gamma);
        //***************************************************************************************


        local_position ++;
    }
    //***************************************************************************************

    //***************************************************************************************
    Log3DM[0][0][0] = 0.0;
    Log3DX[0][0][0] = min_inf;
    Log3DY[0][0][0] = min_inf;
    //***************************************************************************************
    // 2D LK COMPUTATION
    //***************************************************************************************
    // computes the lk in the two subtrees
    std::vector<double> &lk_down_L = dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->log_lk_down_;
    std::vector<double> &lk_down_R = dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->log_lk_down_;

    // MATCH2D
    double pr_m;
    double pr_m_fp;
    for (i = 0; i < h_compr; i++) {
        for (j = 0; j < w_compr; j++) {

            _computeLK_M(dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->fv_data_.at(i),
                         dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->fv_data_.at(j),
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

        _computeLK_X(dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->fv_data_.at(i),
                     dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->fv_empty_data_,
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

        _computeLK_Y(dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->fv_empty_data_,
                     dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->fv_data_.at(j),
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
    double tmp_lk1,tmp_lk2,tmp_lk3;
    id1m = map_compr_L->at(0);
    id2m = map_compr_R->at(0);
    Log3DM[1][1][1] = Log2DM[id1m][id2m];
    Log3DX[1][1][0] = Log2DX[id1m];
    Log3DY[1][0][1] = Log2DY[id2m];
    for (m = 2; m < d; m++) {

        //***************************************************************************
        // GAPX[i][0]
        j=0;
        for (i = 1; i < h; i++) {
            id1x = map_compr_L->at(i - 1);

            Log3DM[m][i - 1][j] = min_inf;
            Log3DY[m][i - 1][j] = min_inf;

            if(std::isinf(Log3DX[m-1][i-1][j])){
                Log3DX[m][i][j] = min_inf;
            } else {
                tmp_lk3 = Log2DX[id1x];
                Log3DX[m][i][j] = progressivePIPutils::add_lns(Log3DX[m-1][i-1][j], tmp_lk3);
            }

        }
        //***********************************************************************************
        // GAPY[0][j]
        i=0;
        for (j = 1; j < w; j++) {
            id2y = map_compr_R->at(j - 1);

            Log3DM[m][i][j - 1] = min_inf;
            Log3DX[m][i][j - 1] = min_inf;

            if (std::isinf(Log3DY[m - 1][i][j - 1])) {
                Log3DY[m][i][j] = min_inf;
            } else {
                tmp_lk3 = Log2DY[id2y];
                Log3DY[m][i][j] = progressivePIPutils::add_lns(Log3DY[m - 1][i][j - 1], tmp_lk3);
            }

        }
        //***********************************************************************************
        for (i = 1; i < h; i++) {
            for (j = 1; j < w; j++) {
                //***************************************************************************
                // MATCH[i][j]
                id1m = map_compr_L->at(i-1);
                id2m = map_compr_R->at(j-1);

                tmp_lk1 = progressivePIPutils::add_lns(Log3DM[m - 1][i - 1][j - 1], Log3DX[m - 1][i - 1][j - 1]);
                tmp_lk2 = progressivePIPutils::add_lns(tmp_lk1, Log3DY[m - 1][i - 1][j - 1]);

                if(std::isinf(tmp_lk2)){
                    Log3DM[m][i][j] = min_inf;
                } else {
                    tmp_lk3 = Log2DM[id1m][id2m];
                    Log3DM[m][i][j] = progressivePIPutils::add_lns(tmp_lk2, tmp_lk3);
                }

                //***************************************************************************
                // GAPX[i][j]
                id1x = map_compr_L->at(i-1);

                tmp_lk1 = progressivePIPutils::add_lns(Log3DM[m - 1][i - 1][j], Log3DX[m - 1][i - 1][j]);
                tmp_lk2 = progressivePIPutils::add_lns(tmp_lk1, Log3DY[m - 1][i - 1][j]);

                if(std::isinf(tmp_lk2)){
                    Log3DX[m][i][j] = min_inf;
                } else {
                    tmp_lk3 = Log2DX[id1x];
                    Log3DX[m][i][j] = progressivePIPutils::add_lns(tmp_lk2, tmp_lk3);
                }

                //***************************************************************************
                // GAPY[i][j]
                id2y = map_compr_R->at(j-1);

                tmp_lk1 = progressivePIPutils::add_lns(Log3DM[m - 1][i][j - 1], Log3DX[m - 1][i][j - 1]);
                tmp_lk2 = progressivePIPutils::add_lns(tmp_lk1, Log3DY[m - 1][i][j - 1]);

                if(std::isinf(tmp_lk2)){
                    Log3DY[m][i][j] = min_inf;
                } else {
                    tmp_lk3 = Log2DY[id2y];
                    Log3DY[m][i][j] = progressivePIPutils::add_lns(tmp_lk2, tmp_lk3);
                }
                //***************************************************************************
            }
        }
    }
    //***************************************************************************************
    //int idx_M,idx_X,idx_Y;
    int best_level;
    double best_score;
    double pm,pmn,log_pm,max_M;
    double px,pxn,log_px,max_X;
    double py,pyn,log_py,max_Y;
    double z;
    double lk;
    double p0;
    double random_number;
    double log_P;
    int T;
    int state;
    int idmL,idmR;

    //int local_position = position;
    for(int sb=0;sb<progressivePIP_->num_sb_;sb++) {

        std::vector<int> traceback;
        //std::vector<vector<int> > traceback_map;
        //traceback_map.resize(2);

        i = h - 1;
        j = w - 1;

        std::vector<vector<bpp::ColMatrix<double> > > fv_data_not_compressed;
        std::vector<std::vector<double>> fv_sigma_not_compressed;
        std::vector<double> lk_down_not_compressed;

        //----------------------------------------------------------------------------------
        // TODO: select starting point
        startingLevelSB(Log3DM,
                        Log3DX,
                        Log3DY,
                        SMALL_DOUBLE,
                        generator,
                        distribution,
                        d,
                        h,
                        w,
                        best_level,
                        best_score,
                        state);
        //----------------------------------------------------------------------------------

        switch (state) {
            case MATCH_STATE:
                T = (int) MATCH_STATE;

                idmL = map_compr_L->at(i - 1);
                idmR = map_compr_R->at(j - 1);

                fv_data_not_compressed.push_back(Fv_M[idmL][idmR]);
                fv_sigma_not_compressed.push_back(Fv_sigma_M[idmL][idmR]);
                lk_down_not_compressed.push_back(Log2DM[idmL][idmR]);

                i = i - 1;
                j = j - 1;
                m = best_level - 1;

                traceback.push_back(T);

                //traceback_map.at(LEFT).push_back(i);
                //traceback_map.at(RIGHT).push_back(j);

                break;
            case GAP_X_STATE:
                T = (int) GAP_X_STATE;

                idmL = map_compr_L->at(i - 1);

                fv_data_not_compressed.push_back(Fv_X[idmL]);
                fv_sigma_not_compressed.push_back(Fv_sigma_X[idmL]);
                lk_down_not_compressed.push_back(Log2DX[idmL]);

                i = i - 1;
                m = best_level - 1;

                traceback.push_back(T);

                //traceback_map.at(LEFT).push_back(i);
                //traceback_map.at(RIGHT).push_back(-1);

                break;
            case GAP_Y_STATE:
                T = (int) GAP_Y_STATE;

                idmR = map_compr_R->at(j - 1);

                fv_data_not_compressed.push_back(Fv_Y[idmR]);
                fv_sigma_not_compressed.push_back(Fv_sigma_Y[idmR]);
                lk_down_not_compressed.push_back(Log2DY[idmR]);

                j = j - 1;
                m = best_level - 1;

                traceback.push_back(T);

                //traceback_map.at(LEFT).push_back(-1);
                //traceback_map.at(RIGHT).push_back(j);

                break;
            default:
                LOG(FATAL)
                        << "\nSomething went wrong in reading the STATE value. "
                           "STATE is neither MATCH, nor GAPX, nor GAPY. ";
        }

        double log_Zm = Log3DM[m][i][j];
        double log_Zx = Log3DX[m][i][j];
        double log_Zy = Log3DY[m][i][j];

        if (std::isinf(log_Zm) && std::isinf(log_Zx) && std::isinf(log_Zy)) {
            LOG(FATAL) << "\nlog_Zm,log_Zx and log_Zy are all infinite.";
        }

        double log_Zmx = progressivePIPutils::add_lns(log_Zm, log_Zx);
        double log_Z = progressivePIPutils::add_lns(log_Zmx, log_Zy);

        if (std::isinf(log_Z)) {
            LOG(FATAL) << "\nlog_Z is infinite.";
        }

        if (std::isinf(log_Zm)) {
            pm = 0;
            pmn = 0;
        } else {
            log_pm = log_Zm - log_Z;
            pm = exp(log_pm);
            pmn = exp(-(1 - pm) / temperature);
        }

        if (std::isinf(log_Zx)) {
            px = 0;
            pxn = 0;
        } else {
            log_px = log_Zx - log_Z;
            px = exp(log_px);
            pxn = exp(-(1 - px) / temperature);
        }

        if (std::isinf(log_Zy)) {
            py = 0;
            pyn = 0;
        } else {
            log_py = log_Zy - log_Z;
            py = exp(log_py);
            pyn = exp(-(1 - py) / temperature);
        }

        z = pmn + pxn + pyn;
        pm = pmn / z;
        px = pxn / z;
        py = pyn / z;


        lk = -log(m) + log_nu_gamma + log_phi_gamma;

        //m = m-1;

        while (i > 0 || j > 0 || m > 0) {

            random_number = distribution(generator);

            if (random_number < pm) {

                idmL = map_compr_L->at(i - 1);
                idmR = map_compr_R->at(j - 1);

                fv_data_not_compressed.push_back(Fv_M[idmL][idmR]);
                fv_sigma_not_compressed.push_back(Fv_sigma_M[idmL][idmR]);
                lk_down_not_compressed.push_back(Log2DM[idmL][idmR]);

                log_P = Log2DM[idmL][idmR];

                i = i - 1;
                j = j - 1;
                m = m - 1;

                T = (int) MATCH_STATE;

                traceback.push_back(T);

                //traceback_map.at(LEFT).push_back(i);
                //traceback_map.at(RIGHT).push_back(j);

            } else if (random_number < (pm + px)) {

                idmL = map_compr_L->at(i - 1);

                fv_data_not_compressed.push_back(Fv_X[idmL]);
                fv_sigma_not_compressed.push_back(Fv_sigma_X[idmL]);
                lk_down_not_compressed.push_back(Log2DX[idmL]);

                log_P = Log2DX[idmL];

                i = i - 1;
                m = m - 1;

                T = (int) GAP_X_STATE;

                traceback.push_back(T);

                //traceback_map.at(LEFT).push_back(i);
                //traceback_map.at(RIGHT).push_back(-1);

            } else {

                idmR = map_compr_R->at(j - 1);

                fv_data_not_compressed.push_back(Fv_Y[idmR]);
                fv_sigma_not_compressed.push_back(Fv_sigma_Y[idmR]);
                lk_down_not_compressed.push_back(Log2DY[idmR]);

                log_P = Log2DY[idmR];

                j = j - 1;
                m = m - 1;

                T = (int) GAP_Y_STATE;

                traceback.push_back(T);

                //traceback_map.at(LEFT).push_back(-1);
                //traceback_map.at(RIGHT).push_back(j);

            }

            if (std::isinf(log_P)) {
                LOG(FATAL) << "\nlog_P is infinite.";
            }

            lk = lk + log_P;

            log_Zm = Log3DM[m][i][j];
            log_Zx = Log3DX[m][i][j];
            log_Zy = Log3DY[m][i][j];

            if (std::isinf(log_Zm) && std::isinf(log_Zx) && std::isinf(log_Zy)) {
                LOG(FATAL) << "\nlog_Zm,log_Zx and log_Zy are all infinite.";
            }

            log_Zmx = progressivePIPutils::add_lns(log_Zm, log_Zx);
            log_Z = progressivePIPutils::add_lns(log_Zmx, log_Zy);

            if (std::isinf(log_Z)) {
                LOG(FATAL) << "\nlog_Z is infinite.";
            }

            if (std::isinf(log_Zm)) {
                pm = 0;
                pmn = 0;
            } else {
                log_pm = log_Zm - log_Z;
                pm = exp(log_pm);
                pmn = exp(-(1 - pm) / temperature);
            }

            if (std::isinf(log_Zx)) {
                px = 0;
                pxn = 0;
            } else {
                log_px = log_Zx - log_Z;
                px = exp(log_px);
                pxn = exp(-(1 - px) / temperature);
            }

            if (std::isinf(log_Zy)) {
                py = 0;
                pyn = 0;
            } else {
                log_py = log_Zy - log_Z;
                py = exp(log_py);
                pyn = exp(-(1 - py) / temperature);
            }

            z = pmn + pxn + pyn;
            pm = pmn / z;
            px = pxn / z;
            py = pyn / z;

        }

        reverse(traceback.begin(), traceback.end());

        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(position)->traceback_path_ = traceback;

//        reverse(traceback_map.at(LEFT).begin(), traceback_map.at(LEFT).end());
//        reverse(traceback_map.at(RIGHT).begin(), traceback_map.at(RIGHT).end());

        //MSA_->pipmsa.at(position)->traceback_map_.at(position) = traceback_map;

        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(position)->score_ = best_score;
        //***************************************************************************************
        // BUILD NEW MSA
        //***************************************************************************************
        // converts traceback path into an MSA
        MSA_t *msaL = dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->_getMSA();
        MSA_t *msaR = dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->_getMSA();
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(position)->_build_MSA(*msaL,*msaR);

        //if (position == 0) {
            // assigns the sequence names of the new alligned sequences to the current MSA
        std::vector<string> *seqNameL = &dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->seqNames_;
        std::vector<string> *seqNameR = &dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->seqNames_;
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(position)->_setSeqNameNode(*seqNameL,*seqNameR);
        //}
        //***************************************************************************************
        // COMPRESS INFO
        //***************************************************************************************
        // compress the MSA
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(position)->_compressMSA(progressivePIP_->alphabet_);

        // compress fv values and lk_down
        reverse(fv_data_not_compressed.begin(),fv_data_not_compressed.end());
        reverse(fv_sigma_not_compressed.begin(),fv_sigma_not_compressed.end());

        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(position)->_compressLK(lk_down_not_compressed);

        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(position)->_compressFv(fv_data_not_compressed);

        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(position)->_compressFvSigma(fv_sigma_not_compressed);
        //***************************************************************************************

        position++;

    }

}

bool sortByScore(const bpp::PIPmsa *msa1, const bpp::PIPmsa *msa2) {
    return msa1->score_ > msa2->score_;
}

void nodeSB::DP3D_PIP() {

    if (_isTerminalNode()) {
        // align leaf (prepare data)
        DP3D_PIP_leaf();
    }else{

        // align internal node

        int num_MSA_L = static_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.size();
        int num_MSA_R = static_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.size();

        subMSAidxL_.resize(progressivePIP_->num_sb_);
        subMSAidxR_.resize(progressivePIP_->num_sb_);

        this->MSA_  = new PIPmsaComp(0);

        //========================================================================
        // temporary object
        //========================================================================
        int totalSubMSAs = num_MSA_L * num_MSA_R * progressivePIP_->num_sb_;

        nodeSB *SBnode = new nodeSB(this->progressivePIP_,
                                    this->vnode_,
                                    this->bnode_);

        PIPmsaComp *MSA = new PIPmsaComp(totalSubMSAs);

        for(int i=0;i<totalSubMSAs;i++){
            // create a new PIPmsa
            dynamic_cast<PIPmsaComp *>(MSA)->pipmsa.at(i) = new PIPmsa();
        }

        SBnode->MSA_ = MSA;
        SBnode->subMSAidxL_.resize(totalSubMSAs);
        SBnode->subMSAidxR_.resize(totalSubMSAs);
        SBnode->childL = this->childL;
        SBnode->childR = this->childR;
        SBnode->parent = this->parent;

        SBnode->iotasNode_ = this->iotasNode_;
        SBnode->betasNode_ = this->betasNode_;
        SBnode->alphaNode_ = this->alphaNode_;
        SBnode->etaNode_ = this->etaNode_;
        //========================================================================

        int position = 0;
        for (unsigned int msa_idx_L = 0; msa_idx_L < num_MSA_L; msa_idx_L++) {
            for (unsigned int msa_idx_R = 0; msa_idx_R < num_MSA_R; msa_idx_R++) {
                SBnode->DP3D_PIP_node(msa_idx_L,msa_idx_R,position);
            }
        }

        if(totalSubMSAs<=progressivePIP_->num_sb_){

            for(int i=0;i<totalSubMSAs;i++){
                this->MSA_->add(static_cast<PIPmsaComp *>(SBnode->MSA_)->pipmsa.at(i));
            }

        }else{

            sort(static_cast<PIPmsaComp *>(SBnode->MSA_)->pipmsa.begin(),
                 static_cast<PIPmsaComp *>(SBnode->MSA_)->pipmsa.end(),
                 sortByScore);


            for(int i=0;i<progressivePIP_->num_sb_;i++){
                this->MSA_->add(static_cast<PIPmsaComp *>(SBnode->MSA_)->pipmsa.at(i));
            }

        }

    }

}
