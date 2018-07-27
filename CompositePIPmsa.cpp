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

void PIPmsaSingle::_compressMSA(const bpp::Alphabet *alphabet) {

    auto sequences = new bpp::VectorSequenceContainer(alphabet);

    std::vector<std::string> seqs = compositePIPmsaUtils::siteContainer2sequenceVector(pipmsa->msa_);

    for(int i = 0; i < seqs.size(); i++) {
        sequences->addSequence(*(new bpp::BasicSequence(pipmsa->seqNames_.at(i),
                                                        seqs.at(i),
                                                        alphabet)), true);

        std::cout<<"["<<seqs.at(i)<<"]\n";

    }

    auto siteContainer = new bpp::VectorSiteContainer(*sequences);
    auto siteContCompr = bpp::PatternTools::shrinkSiteSet(*siteContainer);
    auto map_seqs = bpp::PatternTools::getIndexes(*siteContainer, *siteContCompr);

    pipmsa->map_compressed_seqs_ = map_seqs;

    std::vector<int> rev_map_seqs = compositePIPmsaUtils::reverse_map(map_seqs);

    pipmsa->rev_map_compressed_seqs_ = rev_map_seqs;

}

void PIPmsaComp::_compressMSA(const bpp::Alphabet *alphabet,int idx_sb) {

//    auto sequences = new bpp::VectorSequenceContainer(alphabet);
//
//    std::vector<std::string> seqs = compositePIPmsaUtils::siteContainer2sequenceVector(pipmsa->msa_);
//
//    for(int i = 0; i < seqs.size(); i++) {
//        sequences->addSequence(*(new bpp::BasicSequence(pipmsa->seqNames_.at(i),
//                                                        seqs.at(i),
//                                                        alphabet)), true);
//    }
//
//    auto siteContainer = new bpp::VectorSiteContainer(*sequences);
//    auto siteContCompr = bpp::PatternTools::shrinkSiteSet(*siteContainer);
//    auto map_seqs = bpp::PatternTools::getIndexes(*siteContainer, *siteContCompr);
//
//    map_compressed_seqs_.at(idx_sb) = map_seqs;
//
//    std::vector<int> rev_map_seqs = compositePIPmsaUtils::reverse_map(map_seqs);
//
//    rev_map_compressed_seqs_.at(idx_sb) = rev_map_seqs;

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
    std::cout<<"\n\nMSA\n";
    for (int j = 0; j < pipmsa->traceback_path_.size(); j++) {
        std::cout<<pipmsa->msa_.at(j)<<"\n";
    }
    std::cout<<"\n\n";
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

int PIPmsaComp::getMSAlength(int idx){

    return pipmsa.at(idx)->getNumSites();

}

int PIPmsaComp::getCompressedMSAlength(int idx){

    //return pipmsa->rev_map_compressed_seqs_.size();

}

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