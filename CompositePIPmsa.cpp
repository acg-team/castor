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
 * @file pPIP.hpp
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

/*
* From SeqLib:
*/
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Phyl/PatternTools.h>
#include <glog/logging.h>
#include <string>

#include "CompositePIPmsa.hpp"
#include "progressivePIP.hpp"

using namespace bpp;

void PIPmsa::_setSeqNameLeaf(std::string &seqName) {

    // leaf has only one sequence
    seqNames_.push_back(seqName);

}

void PIPmsa::_setMSAleaf(const bpp::Sequence *sequence) {

    // convert from Sequence to string
    std::string sequenceString = sequence->toString();

    /* convert a string into a vector of single char strings */
    msa_.resize(sequenceString.size());
    for (int i = 0; i < sequenceString.size(); i++) {

        const char ch = sequenceString.at(i);

        // check sequence content: if 'X', ' ' or '-' exit with error
        if(ch=='X' || ch==' ' || ch=='-'){
            LOG(FATAL) << "\nERROR sequence contains 'X' or ' ' or '-'";
        }

        // convert character to string column vector
        MSAcolumn_t msa_col(1, ch);

        // assign MSA column to MSA
        msa_.at(i) = msa_col;
    }

}

void PIPmsa::_setSeqNameNode(std::vector<std::string> &seqNamesL,
                                  std::vector<std::string> &seqNamesR) {

    for (int i = 0; i < seqNamesL.size(); i++) {
        seqNames_.push_back(seqNamesL.at(i));
    }

    for (int i = 0; i < seqNamesR.size(); i++) {
        seqNames_.push_back(seqNamesR.at(i));
    }

}

void PIPmsa::_setFVleaf(int numCatg,const bpp::Alphabet *alphabet) {

    // get the number of gamma categories
    //size_t num_gamma_categories = progressivePIP_->numCatg_;

    // get the number of compressed sites
    int lenComprSeqs = rev_map_compressed_seqs_.size();

    int lenAlphabet = alphabet->getSize();

    // resize fv data([site][catg][states])
    fv_data_.resize(lenComprSeqs);
    for (int i = 0; i < lenComprSeqs; i++) {
        fv_data_[i].resize(numCatg);
    }

    int idx;
    // go through all the sites
    for (int i = 0; i < lenComprSeqs; i++) {

        // get the index in the compressed map
        idx = rev_map_compressed_seqs_.at(i);
        MSAcolumn_t s = this->msa_.at(idx);

        // allocate fv column to the ext. alphabet size
        bpp::ColMatrix<double> fv;
        fv.resize(lenAlphabet, 1); // ColMatrix as Nx1 matrix
        bpp::MatrixTools::fill(fv, 0.0); // all zeros

        // check if the sequence contains a "forbidden" char
        if(s[0]=='X' || s[0]==' ' || s[0]=='-'){
            LOG(FATAL) << "\nERROR sequence contains either 'X' or ' ' or '-'";
        }

        // get the char position in the alphabet
        idx = alphabet->charToInt(&s[0]);

        // set to 1 the indicator array at the position of the observed char
        fv(idx, 0) = 1.0;

        // assign the indicator array to all the gamma categories
        for (int catg = 0; catg < numCatg; catg++) {
            fv_data_.at(i).at(catg) = fv;
        }

    }

}

void PIPmsa::_setFVsigmaLeaf(int lenComprSeqs,
                             int numCatg,
                             const bpp::ColMatrix<double> &pi) {

    // get the number of gamma categories
    //size_t numCatg = progressivePIP_->numCatg_;

    // get the length of the compressed input sequences
    //int lenComprSeqs = getCompressedMSAlength();

    // resize the array ([site][numCatg])
    fv_sigma_.resize(lenComprSeqs);

    double fv0;
    // go through all the sites
    for (int site = 0; site < lenComprSeqs; site++) {

        fv_sigma_.at(site).resize(numCatg);

        // go through all the gamma categories
        for(int catg = 0; catg < numCatg; catg++) {

            // compute fv_sigma = fv dot pi
            fv0 = MatrixBppUtils::dotProd(fv_data_.at(site).at(catg),pi);

            fv_sigma_.at(site).at(catg) = fv0;
        }
    }

}

void PIPmsa::_setFVemptyLeaf(int numCatg,const bpp::Alphabet *alphabet) {

    // get number of gamma categories
    //size_t numCatg = progressivePIP_->numCatg_;

    int lenAlphabet = alphabet->getSize();

    // indicator array (all zeros except the observed character)
    bpp::ColMatrix<double> fv;
    fv.resize(lenAlphabet, 1);
    bpp::MatrixTools::fill(fv, 0.0); // all zeros

    // get the gap position in the alphabet
    std::string ch(1, GAP_CHAR);
    int gapIndex = alphabet->charToInt(ch);

    fv(gapIndex, 0) = 1.0; // set gap position to 1

    // for all the gamma categories an array of fv values
    fv_empty_data_.resize(numCatg);

    // assign the indicator array to all gamma categories
    for (int catg = 0; catg < numCatg; catg++) {
        fv_empty_data_.at(catg) = fv;
    }

}

void PIPmsa::_setFVemptyNode(int numCatg,
                             PIPmsa *childL,
                             PIPmsa *childR,
                             std::vector<bpp::RowMatrix<double> > &PrL,
                             std::vector<bpp::RowMatrix<double> > &PrR){

    // number of discrete gamma categories
    //int numCatg = progressivePIP_->numCatg_;

    fv_empty_data_.resize(numCatg);

    // array of lk (for each gamma rate) of a single column full of gaps
    for (int catg = 0; catg < numCatg; catg++) {

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(PrL.at(catg), childL->fv_empty_data_.at(catg), PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(PrR.at(catg), childR->fv_empty_data_.at(catg), PrfvR);

        // fv = PrfvL * PrfvR
        bpp::ColMatrix<double> fv;
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

        fv_empty_data_.at(catg) = fv;
    }

}

void PIPmsa::_setFVsigmaEmptyLeaf(int numCatg) {

    // get the number of gamma categories
    //size_t numCatg = progressivePIP_->numCatg_;

    // allocate memory ([numCatg] x 1)
    fv_empty_sigma_.resize(numCatg);

    for(int catg=0;catg<numCatg;catg++){
        // fv_empty_sigma = fv dot pi
        // fv_empty_sigma is always 0 at the leaves
        fv_empty_sigma_.at(catg) = 0.0;
    }

}

void PIPmsaSingle::_setFVsigmaEmptyNode(int numCatg,
                                        PIPmsa *childL,
                                        PIPmsa *childR,
                                        double bL,
                                        double bR,
                                        const std::vector<double> &mu) {

    // get the number of categories
    //int numCatg = progressivePIP_->numCatg_;

    //double bL = childL->bnode_->getDistanceToFather(); // left child branch length
    //double bR = childR->bnode_->getDistanceToFather(); // right child branch length

    double zetaL;
    double zetaR;

    // resize to the number of categories
    pipmsa->fv_empty_sigma_.resize(numCatg);

    for(int catg=0;catg<numCatg;catg++){

        zetaL = exp(-mu.at(catg) * bL); // pure survival probability on the left child
        zetaR = exp(-mu.at(catg) * bR); // pure survival probability on the right child

        // fv_empty_sigma = dot(fv_empty,pi)
        // which corresponds to
        // fv_empty_sigma = not_survival_L * not_survival_R +
        //                  not_survival_L * survival_R * not_survival_below_R +
        //                  survival_L * not_survival_below_L * not_survival_R +
        //                  survival_L * not_survival_below_L * not_survival_R * not_survival_below_R
        pipmsa->fv_empty_sigma_.at(catg) = \
                            (1 - zetaL) * (1 - zetaR) + \
                            (1 - zetaL) * zetaR * childR->fv_empty_sigma_.at(catg) + \
                            zetaL * childL->fv_empty_sigma_.at(catg)* (1 - zetaR) + \
                            zetaL * childL->fv_empty_sigma_.at(catg) * zetaR * childR->fv_empty_sigma_.at(catg);

    }

}

void PIPmsa::_compress_Fv(std::vector<std::vector<double>> &fv_sigma_not_compressed,
                           std::vector<vector<bpp::ColMatrix<double> > > &fv_data_not_compressed){

    // compress an array of fv values and fv_sigma values

    int comprMSAlen = rev_map_compressed_seqs_.size();

    int id_map;

    fv_data_.resize(comprMSAlen);

    fv_sigma_.resize(comprMSAlen);

    for(int i=0;i<comprMSAlen;i++){
        id_map = rev_map_compressed_seqs_.at(i);

        fv_data_.at(i)=fv_data_not_compressed.at(id_map);

        fv_sigma_.at(i)=fv_sigma_not_compressed.at(id_map);
    }

}

void PIPmsaSingle::_compressMSA(const bpp::Alphabet *alphabet) {

    auto sequences = new bpp::VectorSequenceContainer(alphabet);

    std::vector<std::string> seqs = compositePIPmsaUtils::siteContainer2sequenceVector(pipmsa->msa_);

    for(int i = 0; i < seqs.size(); i++) {
        sequences->addSequence(*(new bpp::BasicSequence(pipmsa->seqNames_.at(i),
                                                        seqs.at(i),
                                                        alphabet)), true);


        //======= DEBUG ====================
        //std::cout<<"["<<seqs.at(i)<<"]\n";
        //======= DEBUG ====================

    }

    auto siteContainer = new bpp::VectorSiteContainer(*sequences);
    auto siteContCompr = bpp::PatternTools::shrinkSiteSet(*siteContainer);
    auto map_seqs = bpp::PatternTools::getIndexes(*siteContainer, *siteContCompr);

    pipmsa->map_compressed_seqs_ = map_seqs;

    std::vector<int> rev_map_seqs = compositePIPmsaUtils::reverse_map(map_seqs);

    pipmsa->rev_map_compressed_seqs_ = rev_map_seqs;

}

void PIPmsaComp::_compressMSA(const bpp::Alphabet *alphabet,int idx_sb) {

    auto sequences = new bpp::VectorSequenceContainer(alphabet);

    std::vector<std::string> seqs = compositePIPmsaUtils::siteContainer2sequenceVector(pipmsa.at(idx_sb)->msa_);

    for(int i = 0; i < seqs.size(); i++) {
        sequences->addSequence(*(new bpp::BasicSequence(pipmsa.at(idx_sb)->seqNames_.at(i),
                                                        seqs.at(i),
                                                        alphabet)), true);


        //======= DEBUG ====================
        //std::cout<<"["<<seqs.at(i)<<"]\n";
        //======= DEBUG ====================

    }

    auto siteContainer = new bpp::VectorSiteContainer(*sequences);
    auto siteContCompr = bpp::PatternTools::shrinkSiteSet(*siteContainer);
    auto map_seqs = bpp::PatternTools::getIndexes(*siteContainer, *siteContCompr);

    pipmsa.at(idx_sb)->map_compressed_seqs_ = map_seqs;

    std::vector<int> rev_map_seqs = compositePIPmsaUtils::reverse_map(map_seqs);

    pipmsa.at(idx_sb)->rev_map_compressed_seqs_ = rev_map_seqs;

}

void PIPmsaSingle::_build_MSA(MSA_t &msaL,MSA_t &msaR) {

    // convert traceback path into an MSA

    // get dimension of the left/right MSA column
    int lenColL = msaL.at(0).size();
    int lenColR = msaR.at(0).size();

    int idx_i = 0;
    int idx_j = 0;
    for (int j = 0; j < pipmsa->traceback_path_.size(); j++) {

        if (pipmsa->traceback_path_.at(j) == (int)MATCH_STATE) {

            // in MATCH case concatenate left_column (from seq1) with right_column (from seq2)
            pipmsa->msa_.push_back(msaL.at(idx_i) + msaR.at(idx_j));
            idx_i++;
            idx_j++;

        } else if (pipmsa->traceback_path_.at(j) == (int)GAP_X_STATE) {

            // in GAPX case concatenate left_column (from seq1) with a column full of gaps (right)
            std::string gapCol(lenColR, GAP_CHAR);
            pipmsa->msa_.push_back(msaL.at(idx_i) + gapCol);
            idx_i++;

        } else if (pipmsa->traceback_path_.at(j) == GAP_Y_STATE) {

            // in GAPY case concatenate a column (left) full of gaps with right_column (from seq2)
            std::string gapCol(lenColL, GAP_CHAR);
            pipmsa->msa_.push_back(gapCol + msaR.at(idx_j));
            idx_j++;

        } else {
            LOG(FATAL) << "\nSomething went wrong during the traceback in function "
                          "pPIP::_build_MSA. Check call stack below.";
        }
    }

    //======== DEBUG =========================================//
//    std::cout<<"\n\nMSA\n";
//    for (int j = 0; j < pipmsa->traceback_path_.size(); j++) {
//        std::cout<<pipmsa->msa_.at(j)<<"\n";
//    }
//    std::cout<<"\n\n";
    //======== DEBUG =========================================//

}

void PIPmsaComp::_build_MSA(MSA_t &msaL,MSA_t &msaR,int idx_sb) {

    // convert traceback path into an MSA

//    bpp::Node *sonLeft = childL->getBnode();
//    int sonLeftID = sonLeft->getId();
//
//    bpp::Node *sonRight = childR->getBnode();
//    int sonRightID = sonRight->getId();
//
//    MSA_t *MSA_L = &(childL->MSA_.at(0));
//    MSA_t *MSA_R = &(childR->MSA_.at(0));
//
//    int lenColL = MSA_L->at(0).size();
//    int lenColR = MSA_R->at(0).size();
//
//    MSA_t MSA;
//
//    int idx_i = 0;
//    int idx_j = 0;
//    for (int j = 0; j < traceback_path_.at(idx_sb).size(); j++) {
//
//        if (traceback_path_.at(idx_sb).at(j) == (int)MATCH_STATE) {
//            MSA.push_back(MSA_L->at(idx_i) + MSA_R->at(idx_j));
//            idx_i++;
//            idx_j++;
//
//        } else if (traceback_path_.at(idx_sb).at(j) == (int)GAP_X_STATE) {
//            std::string gapCol(lenColR, GAP_CHAR);
//            MSA.push_back(MSA_L->at(idx_i) + gapCol);
//            idx_i++;
//
//        } else if (traceback_path_.at(idx_sb).at(j) == GAP_Y_STATE) {
//
//            std::string gapCol(lenColL, GAP_CHAR);
//            MSA.push_back(gapCol + MSA_R->at(idx_j));
//            idx_j++;
//
//        } else {
//            LOG(FATAL) << "\nSomething went wrong during the traceback in function pPIP::_build_MSA. Check call stack below.";
//        }
//    }
//
//    MSA_.at(idx_sb) = MSA;
}

void PIPmsaSingle::_compress_lk_components( std::vector<double> &lk_down_not_compressed,
                                            std::vector<std::vector<bpp::ColMatrix<double> > > &fv_data_not_compressed){


    int comprMSAlen = pipmsa->rev_map_compressed_seqs_.size();

    int id_map;

    //log_lk_down_[nodeID].resize(comprMSAlen);

    //fv_data_[nodeID].resize(comprMSAlen);

    for(int i=0;i<comprMSAlen;i++){
        id_map = pipmsa->rev_map_compressed_seqs_.at(i);

        //log_lk_down_.at(i)=lk_down_not_compressed.at(id_map);

        //fv_data_at(i)=fv_data_not_compressed.at(id_map);
    }

}

void PIPmsaComp::_compress_lk_components(   std::vector<double> &lk_down_not_compressed,
                                            std::vector<std::vector<bpp::ColMatrix<double> > > &fv_data_not_compressed,
                                            int idx_sb){

//    int nodeID = node->getId();
//
//    int comprMSAlen = rev_map_compressed_seqs_.at(nodeID).at(idx_sb).size();
//
//    int id_map;
//
//    log_lk_down_[nodeID].resize(comprMSAlen);
//
//    fv_data_[nodeID].resize(comprMSAlen);
//
//    for(int i=0;i<comprMSAlen;i++){
//        id_map=rev_map_compressed_seqs_.at(nodeID).at(idx_sb).at(i);
//
//        log_lk_down_[nodeID].at(i)=lk_down_not_compressed.at(id_map);
//
//        fv_data_[nodeID].at(idx_sb).at(i)=fv_data_not_compressed.at(id_map);
//    }

}

void PIPmsa::_setTracebackPathleaves() {

    // get the MSA size
    int MSAlen = msa_.size();

    // resize traceback path
    traceback_path_.resize(MSAlen);

    // resize the traceback map
    traceback_map_.resize(MSAlen);

    for(int i = 0; i < MSAlen; i++){
        // assign MATCH STATE to all the sites
        traceback_path_.at(i) = (int)MATCH_STATE;

        // assign the corresponding position in the sequences
        traceback_map_.at(i) = i;
    }

}

int PIPmsa::getNumSites(){

    return msa_.size();

}

int PIPmsaSingle::getMSAlength(){

    return pipmsa->getNumSites();

}

int PIPmsaSingle::getCompressedMSAlength(){

    return pipmsa->rev_map_compressed_seqs_.size();

}

int PIPmsaComp::getCompressedMSAlength(int idx){

    return pipmsa.at(idx)->rev_map_compressed_seqs_.size();

}

int PIPmsaComp::getMSAlength(int idx){

    return pipmsa.at(idx)->getNumSites();

}

//int PIPmsaComp::getCompressedMSAlength(int idx){
//
//    //return pipmsa->rev_map_compressed_seqs_.size();
//
//}

void PIPmsaSingle::add(PIPmsa *x) {

    pipmsa = x;

}

void PIPmsaComp::add(PIPmsa *x) {

    pipmsa.push_back(x);

}

//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************
std::vector<std::string> compositePIPmsaUtils::siteContainer2sequenceVector(std::vector<bpp::MSAcolumn_t> &MSA){

    std::vector<std::string> seqs;

    int len = MSA.size();
    int nseq = MSA.at(0).size();

    seqs.resize(nseq);
    for(int i=0;i<nseq;i++){
        std::string s;
        s.resize(len);
        for(int j=0;j<len;j++){
            s.at(j)=MSA.at(j).at(i);
        }
        seqs[i]=s;
    }

    return seqs;

};

std::vector<int> compositePIPmsaUtils::reverse_map(std::vector<int> &m){

    std::vector<int> rev_m;

    for(int i=0;i<m.size();i++){
        if( (m.at(i)+1) > rev_m.size()){
            if(m.at(i)-rev_m.size() > 0){
                LOG(FATAL) << "\nERROR in reverse_map";
            }
            rev_m.push_back(i);
        }
    }

    return rev_m;
}

bpp::SiteContainer *compositePIPmsaUtils::pPIPmsa2Sites(const bpp::Alphabet *alphabet,
                                                        std::vector<std::string> &seqNames,
                                                        std::vector<std::string> &MSA) {

//    progressivePIP->getRootNode()
//
//    auto MSAs = pip_node->getMSA();
//
//    auto MSA = MSAs.at(idx_sb);

    auto sequences = new bpp::VectorSequenceContainer(alphabet);
//
//    auto seqNames = progressivePIP->getSeqnames(PIPnodeRoot);

    int msaLen = MSA.size();

    int numLeaves = seqNames.size();
    for (int j = 0; j < numLeaves; j++) {
        std::string seqname = seqNames.at(j);
        std::string seqdata;
        seqdata.resize(msaLen);
        for (int i = 0; i < msaLen; i++) {
            seqdata.at(i) = MSA.at(i).at(j);
        }
        sequences->addSequence(*(new bpp::BasicSequence(seqname, seqdata,alphabet)), true);
    }

    return new bpp::VectorSiteContainer(*sequences);
}
//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************