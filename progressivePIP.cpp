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
 * @version 1.0
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

#include "progressivePIP.hpp"
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <glog/logging.h>

#define ERR_STATE (-999)
#define DBL_EPSILON std::numeric_limits<double>::min()
#define MATCH_STATE 1
#define GAP_X_STATE 2
#define GAP_Y_STATE 3
#define STOP_STATE 4
#define LEFT 0
#define RIGHT 1

using namespace bpp;

//***************************************************************************************
/*
void PIPnodeCPU:: DP3D_PIP() {
    std::cout<<"\n aligning with PIPnodeCPU\n";
};

void PIPnodeRAM:: DP3D_PIP() {
    std::cout<<"\n aligning with PIPnodeRAM\n";
};

void PIPnodeInterface:: DP3D_PIP() {

};

void PIPnodeInterface::PIPalignNodeNEW() {

    DP3D_PIP();

}
*/
//***************************************************************************************
void CompositePIPaligner::PIPalignNode(){

}
//***************************************************************************************
CompositePIPaligner::CompositePIPaligner(int numNodes){
    pip_nodes_.resize(numNodes);
}
//***************************************************************************************
void CompositePIPaligner::PIPalign() {

    size_t num_nodes = pip_nodes_.size();

    for (int k = 0; k < num_nodes; k++) {
        // traverses the list of nodes and aligns the MSAs on the left and right side
        // if nodeInterface is a leaf the resulting MSA is the sequence itself

        ApplicationTools::displayGauge(k, num_nodes);

        pip_nodes_[k]->PIPalignNode();

    }

}
//***************************************************************************************
PIPnode::PIPnode(const progressivePIP *pPIP,tshlib::VirtualNode *vnode,bpp::Node *bnode){

    progressivePIP_ = pPIP;

    vnode_ = vnode;
    bnode_ = bnode;

    nodeID_ = bnode->getId();

}
//***************************************************************************************
void PIPnode::_reserve(int numCatg){

    fv_empty_data_.resize(numCatg);

    fv_empty_sigma_.resize(numCatg);

    // insertion probabilities at the given nodeInterface with rate variation (gamma)
    iotasNode_.resize(numCatg);

    // survival probabilities at the given nodeInterface with rate variation (gamma)
    betasNode_.resize(numCatg);

    // substitution/deletion probability matrices at the given nodeInterface with rate variation (gamma)
    prNode_.resize(numCatg);

};
//***************************************************************************************
std::vector< std::vector<std::string> > PIPnode::getMSA() {
    return MSA_;
}
//***************************************************************************************
void PIPnode::_setMSAsequenceNames() {

    std::string seqname = progressivePIP_->sequences_->getSequencesNames().at((int) vnode_->vnode_seqid);

    std::vector<std::string> seqNames;

    seqNames.push_back(seqname);

    seqNames_ = seqNames;

}
//***************************************************************************************
void PIPnode::_setMSAleaves() {

    std::string seqname = progressivePIP_->sequences_->getSequencesNames().at((int) vnode_->vnode_seqid);
    std::string sequence = progressivePIP_->sequences_->getSequence(seqname).toString();

    /* convert a string into a vector of single char strings */
    //std::vector<std::string> msa;
    MSA_t msa;
    msa.resize(sequence.size());
    for (int i = 0; i < sequence.size(); i++) {
        //std::string s(1, MSAin.at(i));
        MSAcolumn_t msa_col(1, sequence.at(i));
        msa.at(i) = msa_col;
    }

    MSA_.at(0) = msa;

}
//***************************************************************************************
void PIPnode::_compressMSA(int idx_sb) {

    MSA_t MSA = MSA_.at(idx_sb);

    auto sequences = new bpp::VectorSequenceContainer(progressivePIP_->alphabet_);

    std::vector<std::string> seqs = progressivePIPutils::siteContainer_2_sequence_vector(MSA);

    for(int i = 0; i < seqs.size(); i++) {
        sequences->addSequence(*(new bpp::BasicSequence(seqNames_.at(i),
                                                        seqs.at(i),
                                                        progressivePIP_->alphabet_)), true);
    }

    auto siteContainer = new bpp::VectorSiteContainer(*sequences);
    auto siteContCompr = bpp::PatternTools::shrinkSiteSet(*siteContainer);
    auto map_seqs = bpp::PatternTools::getIndexes(*siteContainer, *siteContCompr);

    map_compressed_seqs_.at(idx_sb) = map_seqs;

    std::vector<int> rev_map_seqs = progressivePIPutils::reverse_map(map_seqs);

    rev_map_compressed_seqs_.at(idx_sb) = rev_map_seqs;

}
//***************************************************************************************
void PIPnode::_setFVleaf() {

    int idx;

    MSA_t MSA = MSA_.at(0);

    size_t num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

    int lenComprSeqs = rev_map_compressed_seqs_.size();

    fv_data_.resize(lenComprSeqs);

    for (int i = 0; i < lenComprSeqs; i++) {
        fv_data_[i].resize(num_gamma_categories);
    }

    for (int i = 0; i < lenComprSeqs; i++) {

        idx = rev_map_compressed_seqs_.at(0).at(i);
        MSAcolumn_t s = MSA.at(idx);

        bpp::ColMatrix<double> fv;
        fv.resize(progressivePIP_->extendedAlphabetSize_, 1);
        bpp::MatrixTools::fill(fv, 0.0);

        idx = progressivePIP_->alphabet_->charToInt(&s[0]);
        idx = idx < 0 ? progressivePIP_->alphabetSize_ : idx;

        fv(idx, 0) = 1.0;
        for (int catg = 0; catg < num_gamma_categories; catg++) {
            fv_data_.at(0).at(i).at(catg) = fv;
        }

    }

}
//***************************************************************************************
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
//***************************************************************************************
void PIPnode::_setFVsigmaLeaf() {

    size_t num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

    int lenComprSeqs = rev_map_compressed_seqs_.size();

    fv_sigma_.resize(lenComprSeqs);

    double fv0;
    for (int site = 0; site < lenComprSeqs; site++) {

        fv_sigma_.at(site).resize(num_gamma_categories);

        for(int catg = 0; catg < num_gamma_categories; catg++) {

            fv0 = MatrixBppUtils::dotProd(fv_data_.at(0).at(site).at(catg), progressivePIP_->pi_);

            fv_sigma_.at(0).at(site).at(catg) = fv0;
        }
    }

}
//***************************************************************************************
void PIPnode::_setFVsigmaEmptyLeaf() {

    size_t num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

    fv_empty_sigma_.resize(num_gamma_categories);

    double fv0;

    for(int catg = 0; catg < num_gamma_categories; catg++) {

        fv0 = MatrixBppUtils::dotProd(fv_empty_data_.at(catg), progressivePIP_->pi_);

        fv_empty_sigma_.at(catg) = fv0;
    }

}
//***************************************************************************************
void PIPnode::_setTracebackPathleaves() {

    int MSAlen = MSA_.size();

    traceback_path_.resize(1);
    traceback_path_.at(0).resize(MSAlen);

    traceback_map_.resize(1);
    traceback_map_.at(0).resize(1);
    traceback_map_.at(0).at(0).resize(MSAlen);

    for(int i = 0; i < MSAlen; i++){
        traceback_path_.at(0).at(i) = (int)MATCH_STATE;
        traceback_map_.at(0).at(0).at(i) = i;
    }

}
//***************************************************************************************
void PIPnode::PIPalignNode() {

    VLOG(1) << "[progressivePIP] Processing nodeInterface " << bnode_->getId();

    if (bnode_->isLeaf()) {
        //*******************************************************************************
        // ALIGNS LEAVES
        //*******************************************************************************
        // associates the sequence name to the leaf nodeInterface
        // TODO sequence should be on the object
        _setMSAsequenceNames();

        // creates a column containing the sequence associated to the leaf nodeInterface
        // TODO sequence should be on the object
        _setMSAleaves();

        // compresses sequence at the leaves
        _compressMSA(0);

        // computes the indicator values (fv values) at the leaves
        _setFVleaf();

        // computes the indicator value for an empty column at the leaf
        _setFVemptyLeaf();

        // computes dotprod(pi,fv)
        _setFVsigmaLeaf();

        // computes dotprod(pi,fv) for an empty column at the leaf
        _setFVsigmaEmptyLeaf();

        // sets th etraceback path at the leaf
        _setTracebackPathleaves();
        //*******************************************************************************
    } else {
/*

        //*******************************************************************************
        // ALIGNS INTERNAL NODES
        //*******************************************************************************
        // align using progressive 3D DP PIP
        if(flag_fv){
            DP3D_PIP_RAM_FAST(nodeInterface);
        }else {
            if (flag_RAM) {
                // DP3D_PIP_RAM(nodeInterface, local, flag_map, flag_pattern); // local: tree rooted at the given nodeInterface
            } else {
                DP3D_PIP(nodeInterface, local, flag_map); // local: tree rooted at the given nodeInterface
                // DP3D_PIP_no_gamma(nodeInterface, local, flag_map); // local: tree rooted at the given nodeInterface
            }
*/
        }
        //*******************************************************************************

}
//***************************************************************************************
int CompositePIPaligner::getId(){
    //TODO: understand this
}
//***************************************************************************************
void CompositePIPaligner::addPIPcomponent(PIPcomponent *pip_node){

    pip_nodes_.at(pip_node->getId()) = pip_node;

};
//***************************************************************************************
void PIPnode::_computeLocalTau() {

    if(vnode_->isTerminalNode()){
        tau_ = 0.0;
    }else{

//        b0 = tree_->getNode(treemap_.right.at(vnode->getNodeLeft()), false)->getDistanceToFather() +\
                tree_->getNode(treemap_.right.at(vnode->getNodeRight()), false)->getDistanceToFather();

        //tau = _computeTauRecursive(vnode_);
    }


}
//***************************************************************************************
void PIPnode::_computeLocalNu(int numCategories) {

    nu_.resize(numCategories);

    for (int catg = 0; catg < numCategories; catg++) {
        // computes the normalizing constant with discrete rate variation (gamma distribution)
        // nu(r) = lambda * r * (tau + 1/(mu *r))
        nu_.at(catg) = progressivePIP_->lambda_.at(catg) * (tau_ + 1 / progressivePIP_->mu_.at(catg));
    }

}
//***************************************************************************************
void PIPnode::_getPrFromSubstutionModel() {

    if (!bnode_->hasFather()) {
        // root nodeInterface doesn't have Pr
    } else {

        for (int i = 0; i < progressivePIP_->rDist_->getNumberOfCategories(); i++) {
            // substitution/deletion probabilities with rate variation (gamma)
            // Pr = exp( branchLength * rateVariation * Q )
            prNode_.at(i) = progressivePIP_->substModel_->getPij_t(bnode_->getDistanceToFather() * \
                            progressivePIP_->rDist_->getCategory(i));
        }
    }

}
//***************************************************************************************
progressivePIP::progressivePIP(tshlib::Utree *utree,
                               bpp::Tree *tree,
                               bpp::SubstitutionModel *smodel,
                               UtreeBppUtils::treemap &inTreeMap,
                               bpp::SequenceContainer *sequences,
                               bpp::DiscreteDistribution *rDist,
                               long seed) {

    utree_ = utree;
    _setTree(tree);
    substModel_ = smodel;
    treemap_ = inTreeMap;
    sequences_ = sequences;
    rDist_ = rDist;
    alphabet_ = substModel_->getAlphabet();
    alphabetSize_ = alphabet_->getSize() - 1;
    extendedAlphabetSize_ = alphabetSize_ + 1;
    seed_ = seed;

};

void progressivePIP::initializePIP(std::vector<tshlib::VirtualNode *> &list_vnode_to_root,
                                   int num_sb) {

    //***************************************************************************************
    int numNodes = list_vnode_to_root.size();
    int numCatg = rDist_->getNumberOfCategories();
    //***************************************************************************************
    progressivePIP::_reserve(numCatg);
    //***************************************************************************************
    // COMPUTES LAMBDA AND MU WITH GAMMA
    //***************************************************************************************
    // set lambdas with rate variation (gamma distribution)
    _setLambda(substModel_->getParameter("lambda").getValue());
    // set mus with rate variation (gamma distribution)
    _setMu(substModel_->getParameter("mu").getValue());
    //***************************************************************************************
    // SET PI
    //***************************************************************************************
    // copy pi
    _setPi(substModel_->getFrequencies());
    //***************************************************************************************

    compositePIPaligner_ = new CompositePIPaligner(numNodes);





    nodeFactory *nodeFactory = new bpp::nodeFactory();

    nodeInterface *node1 = nodeFactory->getNode(CPU,1);
    nodeInterface *node2 = nodeFactory->getNode(CPU,2);
    nodeInterface *node3 = nodeFactory->getNode(RAM,3);

    CompositeInterface* clientComposite = new Composite(3);

    clientComposite->Add(node1,0);
    clientComposite->Add(node2,1);
    clientComposite->Add(node3,2);

    clientComposite->Align();





    for (auto &vnode:list_vnode_to_root) {

        auto bnode = tree_->getNode(treemap_.right.at(vnode), false);

        PIPnode * pip_node = new PIPnode(this,vnode,bnode);

        pip_node->_reserve(numCatg);

        //***************************************************************************************
        // COMPUTE ALL LOCAL TREE LENGHTS
        //***************************************************************************************
        pip_node->_computeLocalTau();
        //***************************************************************************************
        // COMPUTE ALL LOCAL POISSON NORMALIZING CONSTANTS
        //***************************************************************************************
        pip_node->_computeLocalNu(numCatg);
        //***************************************************************************************
        // GET Qs
        //***************************************************************************************
        // set substitution/deletion probabilities with rate variation (gamma distribution)
        pip_node->_getPrFromSubstutionModel();
        //***************************************************************************************

        compositePIPaligner_->addPIPcomponent(pip_node);

    }
    //***************************************************************************************

    compositePIPaligner_->PIPalign();

    //compositePIPalignerNEW_->PIPalignNEW();


}

void progressivePIP::_setTree(const Tree *tree) {
    tree_ = new TreeTemplate<Node>(*tree);
}

void progressivePIP::_reserve(int numCatg) {

    // insertion rate with rate variation (gamma)
    lambda_.resize(numCatg);

    // deletion rate with rate variation (gamma)
    mu_.resize(numCatg);

}

void progressivePIP::_setLambda(double lambda) {

    // original lambda w/o rate variation
    lambda0_ = lambda;

    // insertion rate with rate variation among site r
    for (int i = 0; i < rDist_->getNumberOfCategories(); i++) {
        // lambda(r) = lambda * r
        lambda_.at(i) = lambda * rDist_->getCategories().at(i);
    }

}

void progressivePIP::_setMu(double mu) {

    // checks division by 0 or very small value
    if (fabs(mu) < SMALL_DOUBLE) {
        PLOG(FATAL) << "ERROR: mu is too small";
    }

    // original mu w/o rate variation
    mu0_ = mu;

    // deletion rate with rate variation among site r
    for (int i = 0; i < rDist_->getCategories().size(); i++) {
        // mu(r) = mu *r
        mu_.at(i) = mu * rDist_->getCategories().at(i);
    }

}

void progressivePIP::_setPi(const Vdouble &pi) {

    // copy pi (steady state frequency distribution)
    // pi is a colMatrix (column array) to simplify the matrix multiplication
    pi_.resize(pi.size(), 1);
    for (int i = 0; i < pi.size(); i++) {
        pi_(i, 0) = pi.at(i);
    }

}

const Alphabet *progressivePIP::getAlphabet() const {
    return alphabet_;
}

bpp::Node *progressivePIP::getRootNode() {
    return tree_->getRootNode();
}

bpp::SiteContainer *progressivePIPutils::pPIPmsa2Sites(bpp::progressivePIP *progressivePIP,bpp::PIPnode *pip_node,int idx_sb) {

//    progressivePIP->getRootNode()
//
//    auto MSAs = pip_node->getMSA();
//
//    auto MSA = MSAs.at(idx_sb);
//
//    auto sequences = new bpp::VectorSequenceContainer(progressivePIP->getAlphabet());
//
//    auto seqNames = progressivePIP->getSeqnames(progressivePIP->getRootNode());
//
//    int msaLen = MSA.size();
//
//    int numLeaves = seqNames.size();
//    for (int j = 0; j < numLeaves; j++) {
//        std::string seqname = seqNames.at(j);
//        std::string seqdata;
//        seqdata.resize(msaLen);
//        for (int i = 0; i < msaLen; i++) {
//            seqdata.at(i) = MSA.at(i).at(j);
//        }
//        sequences->addSequence(*(new bpp::BasicSequence(seqname, seqdata, progressivePIP->getAlphabet())), true);
//    }
//
//    return new bpp::VectorSiteContainer(*sequences);
}

double progressivePIPutils::add_lns(double a_ln, double b_ln) {
    //ln(a + b) = ln{exp[ln(a) - ln(b)] + 1} + ln(b)

    double R;

    if (std::isinf(a_ln) && std::isinf(b_ln)) {
        R = -std::numeric_limits<double>::infinity();
    } else if (std::isinf(a_ln)) {
        R = b_ln;
    } else if (std::isinf(b_ln)) {
        R = a_ln;
    } else if ((abs(a_ln - b_ln) >= 36.043653389117155)) {
        //TODO:check this
        //2^52-1 = 4503599627370495.	log of that is 36.043653389117155867651465390794
        R = max(a_ln, b_ln);
    } else {
        R = log(exp(a_ln - b_ln) + 1) + b_ln;
    }

    return R;
}

void progressivePIPutils::max_val_in_column(double ***M, int depth, int height, int width, double &val, int &level) {

    val = -std::numeric_limits<double>::infinity();
    level = 0;

    for (int k = 0; k < depth; k++) {
        if (M[k][height - 1][width - 1] > val) {
            val = M[k][height - 1][width - 1];
            level = k;
        }
    }


}

std::vector<std::string> progressivePIPutils::siteContainer_2_sequence_vector(std::vector<bpp::MSAcolumn_t> &MSA){

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

std::vector<int> progressivePIPutils::reverse_map(std::vector<int> &m){

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

