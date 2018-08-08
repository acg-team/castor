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

//==============================================================================
bool sortByScore(const bpp::PIPmsa *msa1, const bpp::PIPmsa *msa2) {
    return msa1->score_ > msa2->score_;
}
//==============================================================================
double sum_3_logs(double l1,double l2,double l3){

    double tmp;

    tmp = progressivePIPutils::add_lns(l1,l2);

    return progressivePIPutils::add_lns(tmp,l3);

}
//==============================================================================
void partitionFunction(double temperature,
                       double log_Zm,
                       double log_Zx,
                       double log_Zy,
                       double &pm,
                       double &px,
                       double &py){

    double pmn,log_pm;
    double pxn,log_px;
    double pyn,log_py;

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

    double z = pmn + pxn + pyn;
    pm = pmn / z;
    px = pxn / z;
    py = pyn / z;

}
//==============================================================================
void nodeSB::selectState(LKdata &lkdata,
                          int state,
                          std::vector<int> *map_compr_L,
                          std::vector<int> *map_compr_R,
                          int &i,int &j,int &m,
                          double &log_P,
                          std::vector<int> &traceback,
                          std::vector<vector<bpp::ColMatrix<double> > > &fv_data_not_compressed,
                          std::vector<std::vector<double>> &fv_sigma_not_compressed,
                          std::vector<double> &lk_down_not_compressed){

    int idmL,idmR;

    switch (state) {
        case MATCH_STATE:

            idmL = map_compr_L->at(i - 1);
            idmR = map_compr_R->at(j - 1);

            fv_data_not_compressed.push_back(lkdata.Fv_M[idmL][idmR]);
            fv_sigma_not_compressed.push_back(lkdata.Fv_sigma_M[idmL][idmR]);
            lk_down_not_compressed.push_back(lkdata.Log2DM[idmL][idmR]);

            log_P = lkdata.Log2DM[idmL][idmR];

            i = i - 1;
            j = j - 1;
            m = m - 1;

            traceback.push_back(state);

            break;
        case GAP_X_STATE:

            idmL = map_compr_L->at(i - 1);

            fv_data_not_compressed.push_back(lkdata.Fv_X[idmL]);
            fv_sigma_not_compressed.push_back(lkdata.Fv_sigma_X[idmL]);
            lk_down_not_compressed.push_back(lkdata.Log2DX[idmL]);

            log_P = lkdata.Log2DX[idmL];

            i = i - 1;
            m = m - 1;

            traceback.push_back(state);

            break;
        case GAP_Y_STATE:

            idmR = map_compr_R->at(j - 1);

            fv_data_not_compressed.push_back(lkdata.Fv_Y[idmR]);
            fv_sigma_not_compressed.push_back(lkdata.Fv_sigma_Y[idmR]);
            lk_down_not_compressed.push_back(lkdata.Log2DY[idmR]);

            log_P = lkdata.Log2DY[idmR];

            j = j - 1;
            m = m - 1;

            traceback.push_back(state);

            break;
        default:
            LOG(FATAL)
                    << "\nSomething went wrong in reading the STATE value. "
                       "STATE is neither MATCH, nor GAPX, nor GAPY. ";
    }

}

void nodeSB::startingLevelSB(LKdata &lkdata,
                             double epsilon,
                             std::default_random_engine &generator,
                             std::uniform_real_distribution<double> &distribution,
                             int h,
                             int w,
                             int &lev,
                             double &val,
                             int &state){

    double sumM = 0.0;
    double sumX = 0.0;
    double sumY = 0.0;
    // get sum of columns
    for(int k=0; k<lkdata.d_;k++){
        if(! std::isinf(lkdata.Log3DM[k][h - 1][w - 1])){
            sumM += lkdata.Log3DM[k][h - 1][w - 1];
        }
        if(! std::isinf(lkdata.Log3DX[k][h - 1][w - 1])) {
            sumX += lkdata.Log3DX[k][h - 1][w - 1];
        }
        if(! std::isinf(lkdata.Log3DY[k][h - 1][w - 1])) {
            sumY += lkdata.Log3DY[k][h - 1][w - 1];
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
            mat3D = &lkdata.Log3DM;
            break;
        case GAP_X_STATE:
            random_number *= (-sumX);
            mat3D = &lkdata.Log3DX;
            break;
        case GAP_Y_STATE:
            random_number *= (-sumY);
            mat3D = &lkdata.Log3DY;
            break;
        default:
            LOG(FATAL) << "\nSomething went wrong in reading the STATE value."
                          " STATE is neither MATCH, nor GAPX, nor GAPY. ";
    }

    double cumsum=0.0;
    for(int level=0; level < lkdata.d_; level++){

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

void nodeSB::_forward(LKdata &lkdata,
                      int position){

    //***************************************************************************************
    // FORWARD RECURSION
    //***************************************************************************************

    //***************************************************************************************
    // WORKING VARIABLES
    //***************************************************************************************
    int id1m,id2m;
    int id1x;
    int id2y;
    int i,j,m;
    double tmp_lk;
    double min_inf = -std::numeric_limits<double>::infinity(); // -inf
    //***************************************************************************************
    // GET SONS
    //***************************************************************************************
    int msa_idx_L = subMSAidxL_.at(position);
    int msa_idx_R = subMSAidxR_.at(position);
    std::vector<int> *map_compr_L = &(dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->map_compressed_seqs_);
    std::vector<int> *map_compr_R = &(dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->map_compressed_seqs_);
    //***************************************************************************************
    // DP SIZES
    //***************************************************************************************
    int h = lkdata.h_; // dimension of the alignment on the left side
    int w = lkdata.w_; // dimension of the alignment on the right side
    int d = lkdata.d_; // third dimension of the DP matrix
    //***************************************************************************************
    // 3D DYNAMIC PROGRAMMING
    //***************************************************************************************
    // For each slice of the 3D cube, compute the values of each cell
    //***************************************************************************************
    lkdata.Log3DM[0][0][0] = -0.0;
    lkdata.Log3DX[0][0][0] = -0.0;
    lkdata.Log3DY[0][0][0] = -0.0;
    //***************************************************************************************


    //==== DEBUG ===============
//    std::cout<<"\nmapL:\n";
//    for(int ii=0;ii<map_compr_L->size();ii++){
//        std::cout<<map_compr_L->at(ii)<<" ; ";
//    }
//    std::cout<<"\nmapR:\n";
//    for(int ii=0;ii<map_compr_R->size();ii++){
//        std::cout<<map_compr_R->at(ii)<<" ; ";
//    }
//    std::cout<<"\n";
    //==== DEBUG ===============


    id1m = map_compr_L->at(0);
    id2m = map_compr_R->at(0);
    lkdata.Log3DM[1][1][1] = lkdata.Log2DM[id1m][id2m];
    lkdata.Log3DX[1][1][0] = lkdata.Log2DX[id1m];
    lkdata.Log3DY[1][0][1] = lkdata.Log2DY[id2m];
    for (m = 2; m < d; m++) {
        //***************************************************************************
        // GAPX[i][0]
        j=0;
        for (i = 1; i < h; i++) {

            lkdata.Log3DM[m][i - 1][j] = min_inf;
            lkdata.Log3DY[m][i - 1][j] = min_inf;

            if(std::isinf(lkdata.Log3DX[m-1][i-1][j])){
                lkdata.Log3DX[m][i][j] = min_inf;
            } else {
                id1x = map_compr_L->at(i - 1);
                lkdata.Log3DX[m][i][j] = progressivePIPutils::add_lns(lkdata.Log3DX[m-1][i-1][j], lkdata.Log2DX[id1x]);
            }

        }
        //***********************************************************************************
        // GAPY[0][j]
        i=0;
        for (j = 1; j < w; j++) {

            lkdata.Log3DM[m][i][j - 1] = min_inf;
            lkdata.Log3DX[m][i][j - 1] = min_inf;

            if (std::isinf(lkdata.Log3DY[m - 1][i][j - 1])) {
                lkdata.Log3DY[m][i][j] = min_inf;
            } else {
                id2y = map_compr_R->at(j - 1);
                lkdata.Log3DY[m][i][j] = progressivePIPutils::add_lns(lkdata.Log3DY[m - 1][i][j - 1], lkdata.Log2DY[id2y]);
            }

        }
        //***********************************************************************************
        for (i = 1; i < h; i++) {
            for (j = 1; j < w; j++) {
                //***************************************************************************
                // MATCH[i][j]

                tmp_lk = sum_3_logs(lkdata.Log3DM[m - 1][i - 1][j - 1],
                                    lkdata.Log3DX[m - 1][i - 1][j - 1],
                                    lkdata.Log3DY[m - 1][i - 1][j - 1]);

                if(std::isinf(tmp_lk)){
                    lkdata.Log3DM[m][i][j] = min_inf;
                } else {
                    id1m = map_compr_L->at(i-1);
                    id2m = map_compr_R->at(j-1);
                    lkdata.Log3DM[m][i][j] = progressivePIPutils::add_lns(tmp_lk, lkdata.Log2DM[id1m][id2m]);
                }

                //***************************************************************************
                // GAPX[i][j]

                tmp_lk = sum_3_logs(lkdata.Log3DM[m - 1][i - 1][j],
                                    lkdata.Log3DX[m - 1][i - 1][j],
                                    lkdata.Log3DY[m - 1][i - 1][j]);

                if(std::isinf(tmp_lk)){
                    lkdata.Log3DX[m][i][j] = min_inf;
                } else {
                    id1x = map_compr_L->at(i-1);
                    lkdata.Log3DX[m][i][j] = progressivePIPutils::add_lns(tmp_lk, lkdata.Log2DX[id1x]);
                }

                //***************************************************************************
                // GAPY[i][j]

                tmp_lk = sum_3_logs(lkdata.Log3DM[m - 1][i][j - 1],
                                    lkdata.Log3DX[m - 1][i][j - 1],
                                    lkdata.Log3DY[m - 1][i][j - 1]);

                if(std::isinf(tmp_lk)){
                    lkdata.Log3DY[m][i][j] = min_inf;
                } else {
                    id2y = map_compr_R->at(j-1);
                    lkdata.Log3DY[m][i][j] = progressivePIPutils::add_lns(tmp_lk, lkdata.Log2DY[id2y]);
                }
                //***************************************************************************
            }
        }
    }




    //==== DEBUG ===============
    std::cout<<"\n";
    for(int kk=0;kk<d;kk++){
        std::cout<<"M["<<kk<<"]\n";
        for(int ii=0;ii<h;ii++){
            for(int jj=0;jj<w;jj++){
                double lk;
                if(std::isinf(lkdata.Log3DM[kk][ii][jj])){
                    lk=-0.0;
                }else{
                    lk=lkdata.Log3DM[kk][ii][jj];
                }
                printf("%8.6lf ",lk);
            }
            std::cout<<"\n";
        }
        std::cout<<"\n\n";
    }
    std::cout<<"\n";
    for(int kk=0;kk<d;kk++){
        std::cout<<"X["<<kk<<"]\n";
        for(int ii=0;ii<h;ii++){
            for(int jj=0;jj<w;jj++){
                double lk;
                if(std::isinf(lkdata.Log3DX[kk][ii][jj])){
                    lk=-0.0;
                }else{
                    lk=lkdata.Log3DX[kk][ii][jj];
                }
                printf("%8.6lf ",lk);
            }
            std::cout<<"\n";
        }
        std::cout<<"\n\n";
    }
    std::cout<<"\n";
    for(int kk=0;kk<d;kk++){
        std::cout<<"Y["<<kk<<"]\n";
        for(int ii=0;ii<h;ii++){
            for(int jj=0;jj<w;jj++){
                double lk;
                if(std::isinf(lkdata.Log3DY[kk][ii][jj])){
                    lk=-0.0;
                }else{
                    lk=lkdata.Log3DY[kk][ii][jj];
                }
                printf("%8.6lf ",lk);
            }
            std::cout<<"\n";
        }
        std::cout<<"\n\n";
    }
    //==== DEBUG ===============

}

void nodeSB::_backward(LKdata &lkdata,
                       int position){

    //TODO remove
    double temperature = 1.0;

    //***************************************************************************************
    // BACKWARD RECURSION
    //***************************************************************************************

    //***************************************************************************************
    // WORKING VARIABLES
    //***************************************************************************************
    int i,j,m;
    int best_level;
    double best_score;
    double pm;
    double px;
    double py;
    double lk;
    double random_number;
    double log_P;
    int state;
    //***************************************************************************************
    // GET SONS
    //***************************************************************************************
    int msa_idx_L = subMSAidxL_.at(position);
    int msa_idx_R = subMSAidxR_.at(position);
    std::vector<int> *map_compr_L = &(dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->map_compressed_seqs_);
    std::vector<int> *map_compr_R = &(dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->map_compressed_seqs_);
    //***************************************************************************************
    // DP SIZES
    //***************************************************************************************
    // Compute dimensions of the 3D block at current internal node.
    int h = dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->_getMSAlength() + 1; // dimension of the alignment on the left side
    int w = dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->_getMSAlength() + 1; // dimension of the alignment on the right side
    int d = (h - 1) + (w - 1) + 1; // third dimension of the DP matrix
    //***************************************************************************************
    // RANDOM NUMBERS GENERATOR
    //***************************************************************************************
    std::default_random_engine generator(progressivePIP_->getSeed()); // jatiapp seed
    std::uniform_real_distribution<double> distribution(0.0,1.0); // uniform distribution for the selection
    //***************************************************************************************
    // GAMMA VARIABLES
    //***************************************************************************************
    // number of discrete gamma categories
    size_t numCatg = progressivePIP_->numCatg_;
    //***************************************************************************************
    // LK COMPUTATION OF AN EMPTY COLUMNS (FULL OF GAPS)
    //***************************************************************************************
    // computes the lk of an empty column in the two subtrees

    std::vector<bpp::ColMatrix<double> > &fvL = dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->fv_empty_data_;
    std::vector<bpp::ColMatrix<double> > &fvR = dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->fv_empty_data_;

    std::vector<double> &lk_emptyL = dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->lk_empty_;
    std::vector<double> &lk_emptyR = dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->lk_empty_;

    double nu_gamma;
    double log_phi_gamma;
    double log_nu_gamma;

    int local_position = position;
    for(int sb=0;sb<progressivePIP_->num_sb_;sb++) {

        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->lk_empty_.resize(numCatg);
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->fv_empty_data_.resize(numCatg);
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->fv_empty_sigma_.resize(numCatg);

        std::vector<bpp::ColMatrix<double> > &fv_empty_data = dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->fv_empty_data_;
        std::vector<double> &fv_empty_sigma = dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->fv_empty_sigma_;

        std::vector<double> &lk_empty = dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->lk_empty_;

        std::vector<double> pc0 = _computeLkEmptyNode(fvL,
                                                      fvR,
                                                      fv_empty_data,
                                                      fv_empty_sigma,
                                                      lk_emptyL,
                                                      lk_emptyR,
                                                      lk_empty);

        //***********************************************************************************
        // COMPUTES LOG(PHI(0))
        //***********************************************************************************
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
        //***********************************************************************************


        local_position ++;
    }
    //***************************************************************************************

    //***************************************************************************************
    local_position = position;
    for(int sb=0;sb<progressivePIP_->num_sb_;sb++) {

        std::vector<int> traceback;
        std::vector<vector<bpp::ColMatrix<double> > > fv_data_not_compressed;
        std::vector<std::vector<double>> fv_sigma_not_compressed;
        std::vector<double> lk_down_not_compressed;

        //***********************************************************************************
        startingLevelSB(lkdata,
                        SMALL_DOUBLE,
                        generator,
                        distribution,
                        h,
                        w,
                        best_level,
                        best_score,
                        state);

        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->score_ = best_score;

        i = h - 1;
        j = w - 1;
        m = best_level;
        //***********************************************************************************
        selectState(lkdata,
                    state,
                    map_compr_L,
                    map_compr_R,
                    i,j,m,
                    log_P,
                    traceback,
                    fv_data_not_compressed,
                    fv_sigma_not_compressed,
                    lk_down_not_compressed);
        //***********************************************************************************
        partitionFunction(temperature,
                          lkdata.Log3DM[m][i][j],
                          lkdata.Log3DX[m][i][j],
                          lkdata.Log3DY[m][i][j],
                          pm,
                          px,
                          py);

        lk = -log(m) + log_nu_gamma + log_phi_gamma;
        //***********************************************************************************
        while (i > 0 || j > 0 || m > 0) {

            random_number = distribution(generator);

            if (random_number < pm) {
                state = (int) MATCH_STATE;
            } else if (random_number < (pm + px)) {
                state = (int) GAP_X_STATE;
            } else {
                state = (int) GAP_Y_STATE;
            }

            selectState(lkdata,
                        state,
                        map_compr_L,
                        map_compr_R,
                        i,j,m,
                        log_P,
                        traceback,
                        fv_data_not_compressed,
                        fv_sigma_not_compressed,
                        lk_down_not_compressed);

            if (std::isinf(log_P)) {
                LOG(FATAL) << "\nlog_P is infinite.";
            }

            lk = lk + log_P;

            partitionFunction(temperature,
                              lkdata.Log3DM[m][i][j],
                              lkdata.Log3DX[m][i][j],
                              lkdata.Log3DY[m][i][j],
                              pm,
                              px,
                              py);


        }
        //***********************************************************************************
        reverse(traceback.begin(), traceback.end());
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->traceback_path_ = traceback;
        //***********************************************************************************
        // BUILD NEW MSA
        //***********************************************************************************
        // converts traceback path into an MSA
        MSA_t *msaL = dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->_getMSA();
        MSA_t *msaR = dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->_getMSA();
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->_build_MSA(*msaL,*msaR);

        std::vector<string> *seqNameL = &dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L)->seqNames_;
        std::vector<string> *seqNameR = &dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R)->seqNames_;
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->_setSeqNameNode(*seqNameL,*seqNameR);
        //***********************************************************************************
        // COMPRESS INFO
        //***********************************************************************************
        // compress the MSA
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->_compressMSA(progressivePIP_->alphabet_);

        // compress fv values and lk
        reverse(fv_data_not_compressed.begin(),fv_data_not_compressed.end());
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->_compressFv(fv_data_not_compressed);

        reverse(fv_sigma_not_compressed.begin(),fv_sigma_not_compressed.end());
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->_compressFvSigma(fv_sigma_not_compressed);

        reverse(lk_down_not_compressed.begin(),lk_down_not_compressed.end());
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(local_position)->_compressLK(lk_down_not_compressed);
        //***********************************************************************************

        local_position++;

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

void nodeSB::DP3D_PIP_node(int position) {

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
    int msa_idx_L = subMSAidxL_.at(position);
    int msa_idx_R = subMSAidxR_.at(position);
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
    // WORKING VARIABLES
    //***************************************************************************************
    int i,j,k;
    //***************************************************************************************
    // MEMORY ALLOCATION
    //***************************************************************************************
    // Initialisation of the data structure
    LKdata lkdata(d,h,h_compr,w,w_compr,numCatg,false);
    //***************************************************************************************
    // 2D LK COMPUTATION
    //***************************************************************************************
    PIPmsa *pipmsaL = dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_L);
    PIPmsa *pipmsaR = dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_R);

//    _DP2D(pipmsaL,
////          pipmsaR,
////          h_compr,
////          w_compr,
////          Log2DM,
////          Log2DX,
////          Log2DY,
////          Log2DM_fp,
////          Log2DX_fp,
////          Log2DY_fp,
////          Fv_M,
////          Fv_X,
////          Fv_Y,
////          Fv_sigma_M,
////          Fv_sigma_X,
////          Fv_sigma_Y);

    _DP2D(pipmsaL,
           pipmsaR,
           lkdata);
    //***************************************************************************************
    // FORWARD RECURSION
    //***************************************************************************************
    _forward(lkdata,position);
    //***************************************************************************************
    // BACKWARD RECURSION
    //***************************************************************************************
    _backward(lkdata,position);
    //***************************************************************************************

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

                SBnode->subMSAidxL_.at(position) = msa_idx_L;
                SBnode->subMSAidxR_.at(position) = msa_idx_R;

                SBnode->DP3D_PIP_node(position);

                position += progressivePIP_->num_sb_;
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
