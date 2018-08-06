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

//void progressivePIPutils::max_val_in_column(double ***M, int depth, int height, int width, double &val, int &level) {
//
//    val = -std::numeric_limits<double>::infinity();
//    level = 0;
//
//    for (int k = 0; k < depth; k++) {
//        if (M[k][height - 1][width - 1] > val) {
//            val = M[k][height - 1][width - 1];
//            level = k;
//        }
//    }
//
//
//}

void nodeSB::_computeLkEmptyLeaf(){

    // compute the lk of an empty column at the leaf

    // get the number of gamma categories
    int numCatg = progressivePIP_->numCatg_;

    // allocate memory ([numCatg] x 1)
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->lk_empty_.resize(numCatg);

    // only 1 column
    for (int catg=0; catg<numCatg; catg++) {
        // compute the lk of an empty column at the leaf
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->lk_empty_.at(catg) = progressivePIP_->rDist_->getProbability((size_t) catg) * \
            iotasNode_.at(catg) * (1 - betasNode_.at(catg));
    }

}

void nodeSB::_computeLkLeaf(){

    // compute the lk at the leaf

    // get the number of gamma categories
    int numCatg = progressivePIP_->numCatg_;

    // get the size of the compressed sequences
    int msaLen = dynamic_cast<PIPmsaComp *>(MSA_)->getCompressedMSAlength(0);

    // allocate memory ([site])
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->log_lk_down_.resize(msaLen);

    // compute the marginal lk over all the gamma categories
    for(int site=0;site<msaLen;site++){
        // init to 0.0
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->log_lk_down_.at(site) = 0.0;
        for (int catg=0; catg<numCatg; catg++) {
            dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->log_lk_down_.at(site) += progressivePIP_->rDist_->getProbability((size_t) catg) * \
                                     iotasNode_.at(catg) * betasNode_.at(catg) * \
                                     dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->fv_sigma_.at(site).at(catg);
        }
        // compute the log lk
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->log_lk_down_.at(site) = log(dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->log_lk_down_.at(site));
    }

}

/*
void nodeSB::_computeAllFvEmptySigmaRec(){

    if(childL == nullptr && childR == nullptr){ // leaf

        // set fv_empty
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_setFVemptyLeaf(progressivePIP_->numCatg_,
                                            progressivePIP_->alphabet_);

        // set fv_sigma_empty = fv_empty dot pi
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_setFVsigmaEmptyLeaf(progressivePIP_->numCatg_);

    }else{ // internal node

        // recursive call
        childL->_computeAllFvEmptySigmaRec();
        childR->_computeAllFvEmptySigmaRec();

        for(int i=0;i<dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.size();i++){
            // set fv_empty
            dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(i)->_setFVemptyNode(progressivePIP_->numCatg_,
                                                dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(0), // 0 or any of them, they are all the same
                                                dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(0), // 0 or any of them, they are all the same
                                          childL->prNode_,
                                          childR->prNode_);

            // set fv_sigma_empty = fv_empty dot pi
            dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(i)->_setFVsigmaEmptyNode(progressivePIP_->numCatg_,
                                                     dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(0), // 0 or any of them, they are all the same
                                                     dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(0), // 0 or any of them, they are all the same
                                               childL->bnode_->getDistanceToFather(),
                                               childR->bnode_->getDistanceToFather(),
                                               progressivePIP_->mu_);
        }

    }

}
*/

std::vector<double> nodeSB::_computeLkEmptyNode(){

    // number of discrete gamma categories
    int numCatg = progressivePIP_->numCatg_;

    double p0;
    double pL,pR;

    // resize array ([numCatg] x 1)
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->lk_empty_.resize(numCatg);

    // array of lk (for each gamma rate) of a single column full of gaps
    std::vector<double> pc0;
    pc0.resize(numCatg);

    for (int catg = 0; catg < numCatg; catg++) {

        //double fv0 = fv_empty_sigma_.at(catg);

        double pr_up = 0.0; // lk_empty UP (from the actual node to the root)

        if(_isRootNode()){ // root
            // lk at root node (beta = 1.0)
            p0 = iotasNode_.at(catg) * dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->fv_empty_sigma_.at(catg);
        }else{ // internal node
            p0 = ( iotasNode_.at(catg) - \
                   iotasNode_.at(catg) * betasNode_.at(catg) + \
                   iotasNode_.at(catg) * betasNode_.at(catg) * dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->fv_empty_sigma_.at(catg) );

            // climb the tree and compute the probability UP
            bpp::PIPnode *tmpNode = this->parent;
            while(tmpNode){
                pr_up += tmpNode->iotasNode_.at(catg) - \
                         tmpNode->iotasNode_.at(catg) * tmpNode->betasNode_.at(catg) + \
                         tmpNode->iotasNode_.at(catg) * tmpNode->betasNode_.at(catg) * dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->fv_empty_sigma_.at(catg);
                tmpNode = tmpNode->parent;
            }

        }

        pL = dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(0)->lk_empty_.at(catg); // lk_empty DOWN left
        pR = dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(0)->lk_empty_.at(catg); // lk_empty DOWN right

        pc0.at(catg) = pr_up + p0 + pL + pR; // this lk_empty is used at this layer

        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->lk_empty_.at(catg) = p0 + pL + pR; // here store the lk for the next layer (probability UP is not added here)
    }

    return pc0;

}

void nodeSB::DP3D_PIP_leaf() {

    //*******************************************************************************
    // ALIGNS LEAVES
    //*******************************************************************************

    // get vnode Id
    int vnodeId = (int) vnode_->vnode_seqid;

    // get sequence name from vnodeId
    std::string seqname = progressivePIP_->sequences_->getSequencesNames().at(vnodeId);

    // create a new PIPmsaComp object
    //MSA_  = new PIPmsaComp();

    // create a new PIPmsa
    //MSA_->pipmsa.resize(1);
    //MSA_->pipmsa.at(0) = new PIPmsa();

    // associates the sequence name to the leaf node
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_setSeqNameLeaf(seqname);

    // get sequence from sequence name
    const bpp::Sequence *sequence = &progressivePIP_->sequences_->getSequence(seqname);

    // creates a column containing the sequence associated to the leaf node
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_setMSAleaf(sequence);

    // compresses sequence at the leaves
    dynamic_cast<PIPmsaComp *>(MSA_)->_compressMSA(progressivePIP_->alphabet_,0);


    // computes the indicator values (fv values) at the leaves
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_setFVleaf(progressivePIP_->numCatg_,
                                   progressivePIP_->alphabet_);

    // computes dotprod(pi,fv)
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_setFVsigmaLeaf(
            dynamic_cast<PIPmsaComp *>(MSA_)->getCompressedMSAlength(0),
                                        progressivePIP_->numCatg_,
                                        progressivePIP_->pi_);

    // compute the lk of an empty column
    _computeLkEmptyLeaf();

    // computes the lk for all the characters at the leaf
    _computeLkLeaf();

    // sets the traceback path at the leaf
    dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_setTracebackPathleaves();

}

void nodeSB::DP3D_PIP_node() {


    //TODO remove
    int msa_idx_L;
    int msa_idx_R;
    int position;
    int num_sb_;
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
    int h = dynamic_cast<PIPmsaComp *>(childL->MSA_)->getMSAlength(msa_idx_L) + 1; // dimension of the alignment on the left side
    int w = dynamic_cast<PIPmsaComp *>(childR->MSA_)->getMSAlength(msa_idx_R) + 1; // dimension of the alignment on the right side
    int d = (h - 1) + (w - 1) + 1; // third dimension of the DP matrix
    int h_compr = dynamic_cast<PIPmsaComp *>(childL->MSA_)->getCompressedMSAlength(msa_idx_L); // dimension of the compressed alignment on the left side
    int w_compr = dynamic_cast<PIPmsaComp *>(childR->MSA_)->getCompressedMSAlength(msa_idx_R); // dimension of the compressed alignment on the right side
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
        nu_gamma += progressivePIP_->rDist_->getProbability((size_t) catg) * progressivePIP_->nu_.at(catg);
        log_phi_gamma += progressivePIP_->rDist_->getProbability((size_t) catg) * (progressivePIP_->nu_.at(catg) * (pc0.at(catg)-1));
    }

    double log_nu_gamma = log(nu_gamma);
    //***************************************************************************************

    //***************************************************************************************
    Log3DM[0][0][0] = 0.0;
    Log3DX[0][0][0] = min_inf;
    Log3DY[0][0][0] = min_inf;
    //***************************************************************************************
    // 2D LK COMPUTATION
    //***************************************************************************************
    // computes the lk in the two subtrees
    std::vector<double> lk_down_L;// = _compute_lk_down(sonLeft,msa_idx_L);
    std::vector<double> lk_down_R;// = _compute_lk_down(sonRight,msa_idx_R);

    // MATCH2D
    double pr_m;
    double pr_m_fp;
    for (i = 0; i < h_compr; i++) {
        for (j = 0; j < w_compr; j++) {

            _computeLK_M(dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(msa_idx_R)->fv_data_.at(i),
                         dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(msa_idx_L)->fv_data_.at(j),
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
    double l1,l2,l3;
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
                l3 = Log2DX[id1x];
                Log3DX[m][i][j] = progressivePIPutils::add_lns(Log3DX[m-1][i-1][j], l3);
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
                l3 = Log2DY[id2y];
                Log3DY[m][i][j] = progressivePIPutils::add_lns(Log3DY[m - 1][i][j - 1], l3);
            }

        }
        //***********************************************************************************
        for (i = 1; i < h; i++) {
            for (j = 1; j < w; j++) {
                //***************************************************************************
                // MATCH[i][j]
                id1m = map_compr_L->at(i-1);
                id2m = map_compr_R->at(j-1);

                l1 = progressivePIPutils::add_lns(Log3DM[m - 1][i - 1][j - 1], Log3DX[m - 1][i - 1][j - 1]);
                l2 = progressivePIPutils::add_lns(l1, Log3DY[m - 1][i - 1][j - 1]);

                if(std::isinf(l2)){
                    Log3DM[m][i][j] = min_inf;
                } else {
                    l3 = Log2DM[id1m][id2m];
                    Log3DM[m][i][j] = progressivePIPutils::add_lns(l2, l3);
                }

                //***************************************************************************
                // GAPX[i][j]
                id1x = map_compr_L->at(i-1);

                l1 = progressivePIPutils::add_lns(Log3DM[m - 1][i - 1][j], Log3DX[m - 1][i - 1][j]);
                l2 = progressivePIPutils::add_lns(l1, Log3DY[m - 1][i - 1][j]);

                if(std::isinf(l2)){
                    Log3DX[m][i][j] = min_inf;
                } else {
                    l3 = Log2DX[id1x];
                    Log3DX[m][i][j] = progressivePIPutils::add_lns(l2, l3);
                }

                //***************************************************************************
                // GAPY[i][j]
                id2y = map_compr_R->at(j-1);

                l1 = progressivePIPutils::add_lns(Log3DM[m - 1][i][j - 1], Log3DX[m - 1][i][j - 1]);
                l2 = progressivePIPutils::add_lns(l1, Log3DY[m - 1][i][j - 1]);

                if(std::isinf(l2)){
                    Log3DY[m][i][j] = min_inf;
                } else {
                    l3 = Log2DY[id2y];
                    Log3DY[m][i][j] = progressivePIPutils::add_lns(l2, l3);
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

    for(int sb=0;sb<num_sb_;sb++) {

        std::vector<int> traceback;
        std::vector<vector<int> > traceback_map;
        traceback_map.resize(2);

        i = h - 1;
        j = w - 1;

        std::vector<vector<bpp::ColMatrix<double> > > fv_data_not_compressed;
        std::vector<std::vector<double>> fv_sigma_not_compressed;

        //----------------------------------------------------------------------------------
        // TODO: select starting point
//        startingLevelSB(Log3DM,
//                        Log3DX,
//                        Log3DY,
//                        epsilon,
//                        generator,
//                        distribution,
//                        d,
//                        h,
//                        w,
//                        best_level,
//                        best_score,
//                        state);
        //----------------------------------------------------------------------------------

        switch (state) {
            case MATCH_STATE:
                T = (int) MATCH_STATE;

                idmL = map_compr_L->at(i - 1);
                idmR = map_compr_R->at(j - 1);

                fv_data_not_compressed.push_back(Fv_M[idmL][idmR]);
                fv_sigma_not_compressed.push_back(Fv_sigma_M[idmL][idmR]);

                i = i - 1;
                j = j - 1;
                m = best_level - 1;

                traceback.push_back(T);

                traceback_map.at(LEFT).push_back(i);
                traceback_map.at(RIGHT).push_back(j);

                break;
            case GAP_X_STATE:
                T = (int) GAP_X_STATE;

                idmL = map_compr_L->at(i - 1);

                fv_data_not_compressed.push_back(Fv_X[idmL]);
                fv_sigma_not_compressed.push_back(Fv_sigma_X[idmL]);

                i = i - 1;
                m = best_level - 1;

                traceback.push_back(T);

                traceback_map.at(LEFT).push_back(i);
                traceback_map.at(RIGHT).push_back(-1);

                break;
            case GAP_Y_STATE:
                T = (int) GAP_Y_STATE;

                idmR = map_compr_R->at(j - 1);

                fv_data_not_compressed.push_back(Fv_Y[idmR]);
                fv_sigma_not_compressed.push_back(Fv_sigma_Y[idmR]);

                j = j - 1;
                m = best_level - 1;

                traceback.push_back(T);

                traceback_map.at(LEFT).push_back(-1);
                traceback_map.at(RIGHT).push_back(j);

                break;
            default:
                LOG(FATAL)
                        << "\nSomething went wrong in reading the STATE value. STATE is neither MATCH, nor GAPX, nor GAPY. ";
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

                log_P = Log2DM[idmL][idmR];

                i = i - 1;
                j = j - 1;
                m = m - 1;

                T = (int) MATCH_STATE;

                traceback.push_back(T);

                traceback_map.at(LEFT).push_back(i);
                traceback_map.at(RIGHT).push_back(j);

            } else if (random_number < (pm + px)) {

                idmL = map_compr_L->at(i - 1);

                fv_data_not_compressed.push_back(Fv_X[idmL]);
                fv_sigma_not_compressed.push_back(Fv_sigma_X[idmL]);

                log_P = Log2DX[idmL];

                i = i - 1;
                m = m - 1;

                T = (int) GAP_X_STATE;

                traceback.push_back(T);

                traceback_map.at(LEFT).push_back(i);
                traceback_map.at(RIGHT).push_back(-1);

            } else {

                idmR = map_compr_R->at(j - 1);

                fv_data_not_compressed.push_back(Fv_Y[idmR]);
                fv_sigma_not_compressed.push_back(Fv_sigma_Y[idmR]);

                log_P = Log2DY[idmR];

                j = j - 1;
                m = m - 1;

                T = (int) GAP_Y_STATE;

                traceback.push_back(T);

                traceback_map.at(LEFT).push_back(-1);
                traceback_map.at(RIGHT).push_back(j);

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

        reverse(traceback_map.at(LEFT).begin(), traceback_map.at(LEFT).end());
        reverse(traceback_map.at(RIGHT).begin(), traceback_map.at(RIGHT).end());
        //MSA_->pipmsa.at(position)->traceback_map_.at(position) = traceback_map;

        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(position)->score_ = best_score;
        //***************************************************************************************
        // BUILD NEW MSA
        //***************************************************************************************
        // converts traceback path into an MSA
        //PIPmsa *msaL = childL->MSA_->_getMSA();
        //PIPmsa *msaR = childR->MSA_->_getMSA();
        //MSA_->_build_MSA(msaL->msa_,msaR->msa_);

        if (position == 0) {
            // assigns the sequence names of the new alligned sequences to the current MSA
            std::vector<string> *seqNameL = &dynamic_cast<PIPmsaComp *>(childL->MSA_)->pipmsa.at(0)->seqNames_;
            std::vector<string> *seqNameR = &dynamic_cast<PIPmsaComp *>(childR->MSA_)->pipmsa.at(0)->seqNames_;
            dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(0)->_setSeqNameNode(*seqNameL,*seqNameR);
        }
        //***************************************************************************************
        // COMPRESS INFO
        //***************************************************************************************
        // compress the MSA
        dynamic_cast<PIPmsaComp *>(MSA_)->_compressMSA(progressivePIP_->alphabet_,position);

        // compress fv values and lk_down
        reverse(fv_data_not_compressed.begin(),fv_data_not_compressed.end());
        reverse(fv_sigma_not_compressed.begin(),fv_sigma_not_compressed.end());

        //_compressLK(lk_down_not_compressed);
        dynamic_cast<PIPmsaComp *>(MSA_)->pipmsa.at(position)->_compress_Fv(fv_sigma_not_compressed,
                                                fv_data_not_compressed);
        //***************************************************************************************

        position++;

    }

}

void nodeSB::DP3D_PIP() {

    if (_isTerminalNode()) {
        // align leaf (prepare data)
        DP3D_PIP_leaf();
    }else{
        // align internal node
        DP3D_PIP_node();
    }

}
