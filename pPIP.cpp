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

#include "pPIP.hpp"
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <glog/logging.h>

#define ERR_STATE (-999)
#define DBL_EPSILON std::numeric_limits<double>::min()
#define MATCH_STATE 1
#define GAP_X_STATE 2
#define GAP_Y_STATE 3
#define STOP_STATE 4

#define MATCH_CHAR '1'
#define GAP_X_CHAR '2'
#define GAP_Y_CHAR '3'

using namespace bpp;

pPIP::pPIP(tshlib::Utree *utree,
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

void pPIP::_reserve(std::vector<tshlib::VirtualNode *> &nodeList) {
    int numNodes = nodeList.size();
    int numCatg = rDist_->getNumberOfCategories();

    // lk score at each node
    score_.resize(numNodes);                                                        // best likelihood score at each node (final)
    score_.assign(numNodes, -std::numeric_limits<double>::infinity());              // initialisation of the best score

    // traceback path at each node
    traceback_path_.resize(numNodes);                                               // vector of strings with the traceback to the path that generates the best MSA

    // sequence names in the MSA at each node
    seqNames_.resize(numNodes);                                                     // it stores the order of the sequences added to the MSA at each node

    // MSA at each node
    MSA_.resize(numNodes);

    // insertion rate with rate variation (gamma)
    lambda_.resize(numCatg);

    // deletion rate with rate variation (gamma)
    mu_.resize(numCatg);

    // normalizing constant with rate variation (gamma)
    nu_.resize(numCatg);

    //
    lk_down_.resize(numNodes);
    lk_empty_down_.resize(numNodes);
    fv_data_.resize(numNodes);

    //======================================
    fv_empty_data_.resize(numNodes);
    for(int i = 0; i < numNodes; i++){
        fv_empty_data_[i].resize(numCatg);
    }
    //======================================

    map_compressed_seqs_.resize(numNodes);
    rev_map_compressed_seqs_.resize(numNodes);

    // Initialise iotas and betas maps
    for (auto &vnode:nodeList) {

        // get node ID
        int nodeID = treemap_.right.at(vnode);

        // insertion probabilities at the given node with rate variation (gamma)
        iotasNode_[nodeID].resize(numCatg);

        // survival probabilities at the given node with rate variation (gamma)
        betasNode_[nodeID].resize(numCatg);

        // substitution/deletion probability matrices at the given node with rate variation (gamma)
        prNode_[nodeID].resize(numCatg);
    }

}

void pPIP::_setTree(const Tree *tree) {
    tree_ = new TreeTemplate<Node>(*tree);
}

std::vector<std::string> pPIP::getMSA(bpp::Node *node) {
    return MSA_.at(node->getId());
}

double pPIP::getScore(bpp::Node *node) {
    return score_.at(node->getId());
}

std::vector<std::string> pPIP::getSeqnames(bpp::Node *node) {
    return seqNames_.at(node->getId());
}

bpp::Node *pPIP::getRootNode() {
    return tree_->getRootNode();
}

const Alphabet *pPIP::getAlphabet() const {
    return alphabet_;
}

bool pPIP::is_inside(int x0, int y0, int xf, int yf, int xt, int yt) {

    if ((xt < x0) || (yt > y0) || (xt > xf) || (yt < yf)) {
        return false;
    }

    if ((y0 - yt) > (xt - x0)) {
        return false;
    }

    return true;
}

void pPIP::set_indeces_M(int &up_corner_i,
                         int &up_corner_j,
                         int &bot_corner_i,
                         int &bot_corner_j,
                         int level,
                         int h,
                         int w) {

    if (level == 0) {
        up_corner_i = 0;
        up_corner_j = 0;
        bot_corner_i = 0;
        bot_corner_j = 0;
    } else {
        up_corner_i = 1 + level - std::min(w - 1, level);
        up_corner_j = std::min(w - 1, level);
        bot_corner_i = std::min(h - 1, level);
        bot_corner_j = 1 + level - std::min(h - 1, level);
    }

}

void pPIP::set_indeces_X(int &up_corner_i,
                         int &up_corner_j,
                         int &bot_corner_i,
                         int &bot_corner_j,
                         int level,
                         int h,
                         int w) {

    if (level == 0) {
        up_corner_i = 0;
        up_corner_j = 0;
        bot_corner_i = 0;
        bot_corner_j = 0;
    } else {
        up_corner_i = 1 + level - 1 - std::min(w - 1, level - 1);
        up_corner_j = std::min(w - 1, level - 1);
        bot_corner_i = std::min(h - 1, level);
        bot_corner_j = level - std::min(h - 1, level);
    }

}

void pPIP::set_indeces_Y(int &up_corner_i,
                         int &up_corner_j,
                         int &bot_corner_i,
                         int &bot_corner_j,
                         int level,
                         int h,
                         int w) {

    if (level == 0) {
        up_corner_i = 0;
        up_corner_j = 0;
        bot_corner_i = 0;
        bot_corner_j = 0;
    } else {
        up_corner_i = level - std::min(w - 1, level);
        up_corner_j = std::min(w - 1, level);
        bot_corner_i = std::min(h - 1, level - 1);
        bot_corner_j = 1 + level - 1 - std::min(h - 1, level - 1);
    }

}

signed long pPIP::get_indices_M(int nx,
                                int ny,
                                int up_corner_i,
                                int up_corner_j,
                                int bot_corner_i,
                                int bot_corner_j,
                                int m,
                                int h,
                                int w) {

    signed long idx;

    set_indeces_M(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, m, h, w);

    if (is_inside(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, nx, ny)) {

        int dx, sx;

        dx = nx - up_corner_i + 1;

        sx = ((dx + 1) * dx / 2) - 1;

        idx = sx + (ny - up_corner_j);
    } else {
        idx = ERR_STATE;
    }

    return idx;

}

signed long pPIP::get_indices_X(int nx,
                                int ny,
                                int up_corner_i,
                                int up_corner_j,
                                int bot_corner_i,
                                int bot_corner_j,
                                int m,
                                int h,
                                int w) {

    signed long idx;

    set_indeces_X(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, m, h, w);

    if (is_inside(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, nx, ny)) {

        int dx, sx;

        dx = nx - up_corner_i + 1;

        sx = ((dx + 1) * dx / 2) - 1;

        idx = sx + (ny - up_corner_j);
    } else {
        idx = ERR_STATE;
    }

    return idx;

}

signed long pPIP::get_indices_Y(int nx,
                                int ny,
                                int up_corner_i,
                                int up_corner_j,
                                int bot_corner_i,
                                int bot_corner_j,
                                int m,
                                int h,
                                int w) {

    signed long idx;

    set_indeces_Y(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, m, h, w);

    if (is_inside(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, nx, ny)) {

        int dx, sx;

        dx = nx - up_corner_i + 1;

        sx = ((dx + 1) * dx / 2) - 1;

        idx = sx + (ny - up_corner_j);
    } else {
        idx = ERR_STATE;
    }

    return idx;

}

void pPIP::set_indeces_T(int &up_corner_i,
                         int &up_corner_j,
                         int &bot_corner_i,
                         int &bot_corner_j,
                         int level,
                         int h,
                         int w) {

    int up_corner_i_x;
    int up_corner_i_y;

    int up_corner_j_x;
    int up_corner_j_y;

    int bot_corner_i_x;
    int bot_corner_i_y;

    int bot_corner_j_x;
    int bot_corner_j_y;

    set_indeces_X(up_corner_i_x, up_corner_j_x, bot_corner_i_x, bot_corner_j_x, level, h, w);

    set_indeces_Y(up_corner_i_y, up_corner_j_y, bot_corner_i_y, bot_corner_j_y, level, h, w);

    int delta_i, delta_j;

    delta_i = bot_corner_i_x - up_corner_i_y;
    delta_j = up_corner_j_y - bot_corner_j_x;

    if (delta_i > delta_j) {
        up_corner_i = up_corner_i_y;
        up_corner_j = up_corner_j_y;
        bot_corner_i = up_corner_i_y + delta_i;
        bot_corner_j = up_corner_j_y - delta_i;
    } else {
        up_corner_i = bot_corner_i_x - delta_j;
        up_corner_j = bot_corner_j_x + delta_j;
        bot_corner_i = bot_corner_i_x;
        bot_corner_j = bot_corner_j_x;
    }

}

void pPIP::reset_corner(int &up_corner_i,
                        int &up_corner_j,
                        int &bot_corner_i,
                        int &bot_corner_j,
                        int h,
                        int w) {

    int delta;

    if (up_corner_j >= w) {
        delta = up_corner_j - w + 1;
        up_corner_j -= delta;
        up_corner_i += delta;
    }
    if (bot_corner_i >= h) {
        delta = bot_corner_i - h + 1;
        bot_corner_i -= delta;
        bot_corner_j += delta;
    }

}

int pPIP::get_indices_T(int nx,
                        int ny,
                        int up_corner_i,
                        int up_corner_j,
                        int bot_corner_i,
                        int bot_corner_j,
                        int m,
                        int h,
                        int w) {

    set_indeces_T(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, m, h, w);

    reset_corner(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, h, w);

    int idx;
    int dx, sx;

    dx = nx - up_corner_i + 1;

    sx = ((dx + 1) * dx / 2) - 1;

    idx = sx + (ny - up_corner_j);

    return idx;

}

bool pPIP::index_of_max(double m,
                        double x,
                        double y,
                        double epsilon,
                        std::default_random_engine &generator,
                        std::uniform_real_distribution<double> &distribution,
                        max_val_str &max_val,
                        bool flag_RAM,
                        int &index,
                        double &val) {

    double random_number;

    if (std::isinf(m) & std::isinf(x) & std::isinf(y))
        if (flag_RAM) {
            max_val.index = int(STOP_STATE);
            max_val.val = -std::numeric_limits<double>::infinity();
            index = int(STOP_STATE);
            val = -std::numeric_limits<double>::infinity();
            return true;
        } else{
            LOG(FATAL)
                    << "\nSomething went wrong during the comparison of m,x,y variables in function pPIP::index_of_max. Check call stack below. ";
            return false;
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
            max_val.index = int(MATCH_STATE);
            max_val.val = m;
            index = int(MATCH_STATE);
            val = m;
            return true;
            //return int(MATCH_STATE);
        } else if (y > m) {
            max_val.index = int(GAP_Y_STATE);
            max_val.val = y;
            index = int(GAP_Y_STATE);
            val = y;
            return true;
            //return int(GAP_Y_STATE);
        } else {
            if (abs(m - y) < epsilon) {
                //m or y
                random_number = distribution(generator);
                if (random_number < (1.0 / 2.0)) {
                    max_val.index = int(MATCH_STATE);
                    max_val.val = m;
                    index = int(MATCH_STATE);
                    val = m;
                    return true;
                    //return int(MATCH_STATE);
                } else {
                    max_val.index = int(GAP_Y_STATE);
                    max_val.val = y;
                    index = int(GAP_Y_STATE);
                    val = y;
                    return true;
                    //return int(GAP_Y_STATE);
                }
            } else {
                LOG(FATAL) << "\nSomething went wrong during the comparison in function pPIP::index_of_max. Check call stack below.";
                return false;
            }
        }
    } else if (x > m) {
        if (x > y) {
            max_val.index = int(GAP_X_STATE);
            max_val.val = x;
            index = int(GAP_X_STATE);
            val = x;
            return true;
            //return int(GAP_X_STATE);
        } else if (y > x) {
            max_val.index = int(GAP_Y_STATE);
            max_val.val = y;
            index = int(GAP_Y_STATE);
            val = y;
            return true;
            //return int(GAP_Y_STATE);
        } else {
            if (abs(x - y) < epsilon) {
                //x or y
                random_number = distribution(generator);
                if (random_number < (1.0 / 2.0)) {
                    max_val.index = int(GAP_X_STATE);
                    max_val.val = x;
                    index = int(GAP_X_STATE);
                    val = x;
                    return true;
                    //return int(GAP_X_STATE);
                } else {
                    max_val.index = int(GAP_Y_STATE);
                    max_val.val = y;
                    index = int(GAP_Y_STATE);
                    val = y;
                    return true;
                    //return int(GAP_Y_STATE);
                }
            } else {
                LOG(FATAL) << "\nSomething went wrong during the comparison in function pPIP::index_of_max. Check call stack below.";
                return false;
            }
        }
    } else {

        double mx = x;
        if (mx > y) {
            //m or x
            random_number = distribution(generator);
            if (random_number < (1.0 / 2.0)) {
                max_val.index = int(MATCH_STATE);
                max_val.val = m;
                index = int(MATCH_STATE);
                val = m;
                return true;
                //return int(MATCH_STATE);
            } else {
                max_val.index = int(GAP_X_STATE);
                max_val.val = x;
                index = int(GAP_X_STATE);
                val = x;
                return true;
                //return int(GAP_X_STATE);
            }
        } else if (y > mx) {
            max_val.index = int(GAP_Y_STATE);
            max_val.val = y;
            index = int(GAP_Y_STATE);
            val = y;
            return true;
            //return int(GAP_Y_STATE);
        } else {
            if (abs(mx - y) < epsilon) {
                //m or x or y
                random_number = distribution(generator);
                if (random_number < (1.0 / 3.0)) {
                    max_val.index = int(MATCH_STATE);
                    max_val.val = m;
                    index = int(MATCH_STATE);
                    val = m;
                    return true;
                    //return int(MATCH_STATE);
                } else if (random_number < (2.0 / 3.0)) {
                    max_val.index = int(GAP_X_STATE);
                    max_val.val = x;
                    index = int(GAP_X_STATE);
                    val = x;
                    return true;
                    //return int(GAP_X_STATE);
                } else {
                    max_val.index = int(GAP_Y_STATE);
                    max_val.val = y;
                    index = int(GAP_Y_STATE);
                    val = y;
                    return true;
                    //return int(GAP_Y_STATE);
                }
            } else {
                LOG(FATAL) << "\nSomething went wrong during the comparison in function pPIP::index_of_max. Check call stack below.";
                return false;
            }
        }
    }

}

double pPIP::max_of_three(double a, double b, double c, double epsilon,bool flag_RAM) {

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
        if(flag_RAM){
            return -std::numeric_limits<double>::infinity();
        }else{
            LOG(FATAL) << "\nSomething went wrong during the comparison in function pPIP::max_of_three. Check call stack below.";
        }

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

bool pPIP::checkboundary(int up_corner_i,
                         int up_corner_j,
                         int bot_corner_i,
                         int bot_corner_j,
                         int h,
                         int w) {

    if ((up_corner_i >= 0) & (up_corner_i < h) & \
       (up_corner_j >= 0) & (up_corner_j < w) & \
       (bot_corner_i >= 0) & (bot_corner_i < h) & \
       (bot_corner_j >= 0) & (bot_corner_j < w)) {
        return true;
    }

    return false;
}

std::string pPIP::createGapCol(int len) {

    // create an MSA column full of gaps
    std::string colMSA(len, '-');

    return colMSA;
}

void pPIP::build_MSA(bpp::Node *node, TracebackPath_t traceback_path) {

    // convert traceback path into an MSA

    tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
    tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

    int sonLeftID = treemap_.right.at(vnode_left);
    int sonRightID = treemap_.right.at(vnode_right);

    MSA_t *MSA_L = &(MSA_.at(sonLeftID));
    MSA_t *MSA_R = &(MSA_.at(sonRightID));

    int lenColL = MSA_L->at(0).size();
    int lenColR = MSA_R->at(0).size();

    MSA_t MSA;

    int idx_i = 0;
    int idx_j = 0;
    for (int j = 0; j < traceback_path.size(); j++) {

        if (traceback_path.at(j) == MATCH_CHAR) {

            MSA.push_back(MSA_L->at(idx_i) + MSA_R->at(idx_j));
            idx_i++;
            idx_j++;

        } else if (traceback_path.at(j) == GAP_X_CHAR) {

            std::string gapCol(lenColR, GAP_CHAR);
            MSA.push_back(MSA_L->at(idx_i) + gapCol);
            idx_i++;

        } else if (traceback_path.at(j) == GAP_Y_CHAR) {

            std::string gapCol(lenColL, GAP_CHAR);
            MSA.push_back(gapCol + MSA_R->at(idx_j));
            idx_j++;

        } else {
            LOG(FATAL) << "\nSomething went wrong during the traceback in function pPIP::build_MSA. Check call stack below.";
        }
    }

    MSA_.at(node->getId()) = MSA;
}

void pPIP::setMSAsequenceNames(bpp::Node *node) {

    tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
    tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

    int sonLeftID = treemap_.right.at(vnode_left);
    int sonRightID = treemap_.right.at(vnode_right);

    std::vector<std::string> seqNames;

    for (int i = 0; i < seqNames_.at(sonLeftID).size(); i++) {
        seqNames.push_back(seqNames_.at(sonLeftID).at(i));
    }

    for (int i = 0; i < seqNames_.at(sonRightID).size(); i++) {
        seqNames.push_back(seqNames_.at(sonRightID).at(i));
    }

    seqNames_.at(node->getId()) = seqNames;

}

void pPIP::setMSAsequenceNames(bpp::Node *node, std::string seqname) {

    std::vector<std::string> seqNames;

    seqNames.push_back(seqname);

    seqNames_.at(node->getId()) = seqNames;

}

void pPIP::setMSAleaves(bpp::Node *node, const std::string &sequence) {

    /* convert a string into a vector of single char strings */
    //std::vector<std::string> msa;
    MSA_t msa;
    msa.resize(sequence.size());
    for (int i = 0; i < sequence.size(); i++) {
        //std::string s(1, MSAin.at(i));
        MSAcolumn_t msa_col(1, sequence.at(i));
        msa.at(i) = msa_col;
    }

    MSA_.at(node->getId()) = msa;

}

void pPIP::_setNu() {

    for (int i = 0; i < rDist_->getNumberOfCategories(); i++) {
        // computes the normalizing constant with discrete rate variation (gamma distribution)
        // nu(r) = lambda * r * (tau + 1/(mu *r))
        nu_.at(i) = lambda_.at(i) * (tau_ + 1 / mu_.at(i));
    }

}

void pPIP::_setLambda(double lambda) {

    // original lambda w/o rate variation
    lambda0 = lambda;

    // insertion rate with rate variation among site r
    for (int i = 0; i < rDist_->getNumberOfCategories(); i++) {
        // lambda(r) = lambda * r
        lambda_.at(i) = lambda * rDist_->getCategories().at(i);
    }

}

void pPIP::_setMu(double mu) {

    // checks division by 0 or very small value
    if (fabs(mu) < SMALL_DOUBLE) {
        PLOG(FATAL) << "ERROR: mu is too small";
    }

    // original mu w/o rate variation
    mu0 = mu;

    // deletion rate with rate variation among site r
    for (int i = 0; i < rDist_->getCategories().size(); i++) {
        // mu(r) = mu *r
        mu_.at(i) = mu * rDist_->getCategories().at(i);
    }

}

void pPIP::_setPi(const Vdouble &pi) {

    // copy pi (steady state frequency distribution)
    // pi is a colMatrix (column array) to simplify the matrix multiplication
    pi_.resize(pi.size(), 1);
    for (int i = 0; i < pi.size(); i++) {
        pi_(i, 0) = pi.at(i);
    }

}

double pPIP::_setTauRecursive(tshlib::VirtualNode *vnode) {

    if (vnode->isTerminalNode()) {
        // return the branch length
        return tree_->getNode(treemap_.right.at(vnode), false)->getDistanceToFather();
    } else {
        // recursive call
        double bl = _setTauRecursive(vnode->getNodeLeft());
        double br = _setTauRecursive(vnode->getNodeRight());

        // actual branch length (at the given node)
        double b0 = tree_->getNode(treemap_.right.at(vnode), false)->getDistanceToFather();

        // sum of actual branch length + total branch length of left subtree + total branch length right subtree
        return b0 + bl + br;
    }

}

void pPIP::_setTau(tshlib::VirtualNode *vnode) {

    // computes the total tree length of the subtree rooted at vnode
    tau_ = _setTauRecursive(vnode->getNodeLeft()) + _setTauRecursive(vnode->getNodeRight());

}

void pPIP::_setAllIotas(bpp::Node *node, bool local_root) {
    double T;

    // recursive function that computes the insertion probability on the actual branch
    //local_root: flag true only the at the first recursive call (local root) then always false

    if (local_root) {

        for (int catg = 0; catg < rDist_->getNumberOfCategories(); catg++) {

            // T(r) = tau + 1/ (mu * r)
            T = tau_ + 1 / mu_.at(catg);

            // checks division by 0 or too small number
            if (fabs(T) < SMALL_DOUBLE) {
                PLOG(WARNING) << "ERROR in set_iota: T too small";
            }

            // iota(root,r) = (lambda * r)/ (mu * r) / (lambda * r * (tau + 1/ (mu * r) ) )
            //           = 1 / (mu * r) / (tau + 1/ (mu *r) )
            iotasNode_[node->getId()][catg] = (1 / mu_.at(catg)) / T;

        }

    } else {

        for (int catg = 0; catg < rDist_->getNumberOfCategories(); catg++) {

            // T(r) = tau + 1/(mu * r)
            T = tau_ + 1 / mu_.at(catg);

            // checks division by 0 or too small number
            if (fabs(T) < SMALL_DOUBLE) {
                PLOG(WARNING) << "ERROR in set_iota: T too small";
            }

            //iotasNode_[node->getId()][catg] = (node->getDistanceToFather() * rDist_->getCategory(catg) ) / T;
            // iota(v,r) = ( lambda * r * b(v) ) / (lambda * r * (tau + 1/ (mu *r) ) )
            //           = b(v) / (tau + 1/ (mu *r) )
            iotasNode_[node->getId()][catg] = node->getDistanceToFather() / T;
        }

    }

    if (!node->isLeaf()) {
        tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
        tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

        int sonLeftID = treemap_.right.at(vnode_left);
        bpp::Node *sonLeft = tree_->getNode(sonLeftID);

        int sonRightID = treemap_.right.at(vnode_right);
        bpp::Node *sonRight = tree_->getNode(sonRightID);

        _setAllIotas(sonLeft, false);  // false: only at the first call local_root=true (first node is the actual root)
        _setAllIotas(sonRight, false); // false: only at the first call local_root=true (first node is the actual root)

    }

}

void pPIP::_setAllBetas(bpp::Node *node, bool local_root) {

    // recursive function that computes the survival probability on the actual branch
    //local_root: flag true only the at the first recursive call (local root) then always false

    if (local_root) {

        for (int catg = 0; catg < rDist_->getNumberOfCategories(); catg++) {
            // by definition at the root (ev. local root) the survival probabiloty is 1
            betasNode_[node->getId()][catg] = 1.0;
        }

    } else {

        for (int catg = 0; catg < rDist_->getCategories().size(); catg++) {

            // muT(r) = r * mu * b(v)
            double muT = rDist_->getCategory(catg) * mu_.at(catg) * node->getDistanceToFather();

            // checks division by 0 or too small value
            if (fabs(muT) < SMALL_DOUBLE) {
                perror("ERROR mu * T is too small");
            }
            // survival probability on node v (different from (local)-root)
            // beta(v,r) = (1 - exp( -mu * r * b(v) )) / (mu * r * b(v))
            betasNode_[node->getId()][catg] = (1.0 - exp(-muT)) / muT;
        }
    }


    if (!node->isLeaf()) {

        tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
        tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

        int sonLeftID = treemap_.right.at(vnode_left);
        bpp::Node *sonLeft = tree_->getNode(sonLeftID);

        int sonRightID = treemap_.right.at(vnode_right);
        bpp::Node *sonRight = tree_->getNode(sonRightID);

        _setAllBetas(sonLeft, false); // false: only at the first call local_root=true (first node is the actual root)
        _setAllBetas(sonRight, false); // false: only at the first call local_root=true (first node is the actual root)

    }

}

void pPIP::_getPrFromSubstutionModel(std::vector<tshlib::VirtualNode *> &listNodes) {

    for (auto &vnode:listNodes) {

        auto node = tree_->getNode(treemap_.right.at(vnode), false);

        if (!node->hasFather()) {
            // root node doesn't have Pr
        } else {
            for (int i = 0; i < rDist_->getNumberOfCategories(); i++) {
                // substitution/deletion probabilities with rate variation (gamma)
                // Pr = exp( branchLength * rateVariation * Q )
                prNode_[node->getId()].at(i) = substModel_->getPij_t(node->getDistanceToFather() * rDist_->getCategory(i));
            }
        }

    }

}

bpp::ColMatrix<double> pPIP::fv_observed(MSAcolumn_t &s, int &idx) {

    // TODO: indicator functions in the initialization, here only retrieval of values

    // fills the indicator function (array) I
    // I is an array of zeros and a 1 only at the position of the observed char

    bpp::ColMatrix<double> fv;
    int ii;
    char ch = s[idx];

    fv.resize(extendedAlphabetSize_, 1);
    bpp::MatrixTools::fill(fv, 0.0);

    ii = alphabet_->charToInt(&ch);
    ii = ii < 0 ? alphabetSize_ : ii;

    fv(ii, 0) = 1.0;
    idx++;

    return fv;
}

bpp::ColMatrix<double> pPIP::computeFVrec(bpp::Node *node, MSAcolumn_t &s, int &idx, int catg) {

    bpp::ColMatrix<double> fv;

    if (node->isLeaf()) {

        fv = fv_observed(s, idx);

    } else {

        tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
        tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

        int sonLeftID = treemap_.right.at(vnode_left);
        bpp::Node *sonLeft = tree_->getNode(sonLeftID);

        int sonRightID = treemap_.right.at(vnode_right);
        bpp::Node *sonRight = tree_->getNode(sonRightID);

        // computes the recursive Felsenstein's peeling weight on the left subtree
        bpp::ColMatrix<double> fvL = computeFVrec(sonLeft, s, idx, catg);

        // computes the recursive Felsenstein's peeling weight on the right subtree
        bpp::ColMatrix<double> fvR = computeFVrec(sonRight, s, idx, catg);

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(prNode_[sonLeftID].at(catg), fvL, PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(prNode_[sonRightID].at(catg), fvR, PrfvR);

        // fv = PrfvL * PrfvR
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

    }

    return fv;
}

void pPIP::allgaps(bpp::Node *node, std::string &s, int &idx, bool &flag) {

    // flag is true if all the leaves of the subtree rooted in node contain a gap

    if (node->isLeaf()) {
        char ch = s[idx];

        idx++;

        if (ch != '-') {
            flag = false;
        }

    } else {

        tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
        tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

        int sonLeftID = treemap_.right.at(vnode_left);
        bpp::Node *sonLeft = tree_->getNode(sonLeftID);

        int sonRightID = treemap_.right.at(vnode_right);
        bpp::Node *sonRight = tree_->getNode(sonRightID);

        allgaps(sonLeft, s, idx, flag);
        allgaps(sonRight, s, idx, flag);
    }

}

double pPIP::compute_lk_gap_down(bpp::Node *node, MSAcolumn_t &s, int catg) {

    int idx;

    int nodeID = node->getId();

    if (node->isLeaf()) {

        idx = 0;
        bpp::ColMatrix<double> fv = computeFVrec(node, s, idx, catg);

        // fv0 = pi * fv
        double fv0 = MatrixBppUtils::dotProd(fv, pi_);

        // lk of non survival till the leaves
        // lk = iota(v,r) - iota(v,r)*beta(v,r) + iota(v,r)*beta(v,r)*fv
        double pr = iotasNode_[nodeID][catg] - \
             iotasNode_[nodeID][catg] * betasNode_[nodeID][catg] + \
             iotasNode_[nodeID][catg] * betasNode_[nodeID][catg] * fv0;

        return pr;

    }

    tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
    tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

    int sonLeftID = treemap_.right.at(vnode_left);
    bpp::Node *sonLeft = tree_->getNode(sonLeftID);

    int sonRightID = treemap_.right.at(vnode_right);
    bpp::Node *sonRight = tree_->getNode(sonRightID);

    idx = 0;
    bpp::ColMatrix<double> fv = computeFVrec(node, s, idx, catg);

    // fv0 = pi * fv
    double fv0 = MatrixBppUtils::dotProd(fv, pi_);

    // lk of non survival till the leaves
    // lk = iota(v,r) - iota(v,r)*beta(v,r) + iota(v,r)*beta(v,r)*fv
    double pr = iotasNode_[nodeID][catg] - \
         iotasNode_[nodeID][catg] * betasNode_[nodeID][catg] + \
         iotasNode_[nodeID][catg] * betasNode_[nodeID][catg] * fv0;


    bool flagL = true;
    bool flagR = true;
    idx = 0;
    allgaps(sonLeft, s, idx, flagL);

    int ixx = idx;
    allgaps(sonRight, s, idx, flagR);

    int len;

    MSAcolumn_t sL;
    len = ixx;
    sL = s.substr(0, len);
    double pL = compute_lk_gap_down(sonLeft, sL, catg);

    MSAcolumn_t sR;
    sR = s.substr(ixx);
    double pR = compute_lk_gap_down(sonRight, sR, catg);

    return pr + pL + pR;
}

double pPIP::compute_lk_down(bpp::Node *node, MSAcolumn_t &s, int catg) {

    int idx;
    //bpp::ColMatrix<double> fvL;
    //bpp::ColMatrix<double> fvR;

    int nodeID = node->getId();

    if (node->isLeaf()) {

        idx = 0;
        bpp::ColMatrix<double> fv = computeFVrec(node, s, idx, catg);

        // fv0 = pi * fv
        double fv0 = MatrixBppUtils::dotProd(fv, pi_);

        double pr = iotasNode_[nodeID][catg] * betasNode_[nodeID][catg] * fv0;

        return pr;

    }

    tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
    tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

    int sonLeftID = treemap_.right.at(vnode_left);
    bpp::Node *sonLeft = tree_->getNode(sonLeftID);

    int sonRightID = treemap_.right.at(vnode_right);
    bpp::Node *sonRight = tree_->getNode(sonRightID);

    idx = 0;
    bpp::ColMatrix<double> fv = computeFVrec(node, s, idx, catg);

    // fv0 = pi * fv
    double fv0 = MatrixBppUtils::dotProd(fv, pi_);

    double pr = iotasNode_[nodeID][catg] * betasNode_[nodeID][catg] * fv0;

    bool flagL = true;
    bool flagR = true;
    idx = 0;
    allgaps(sonLeft, s, idx, flagL);
    int ixx = idx;
    allgaps(sonRight, s, idx, flagR);

    int len;
    if (flagR) {
        std::string sL;
        len = ixx;
        sL = s.substr(0, len);
        return pr + compute_lk_down(sonLeft, sL, catg);
    }

    if (flagL) {
        std::string sR;
        sR = s.substr(ixx);
        return pr + compute_lk_down(sonRight, sR, catg);
    }

    return pr;
}

std::vector<double> pPIP::computeLK_GapColumn_local(int nodeID,
                                                    int sonLeftID,
                                                    int sonRightID,
                                                    std::vector< bpp::ColMatrix<double> > &fvL,
                                                    std::vector< bpp::ColMatrix<double> > &fvR,
                                                    std::vector< bpp::ColMatrix<double> > &Fv_gap) {

    // number of discrete gamma categories
    int num_gamma_categories = rDist_->getNumberOfCategories();

    // array of lk (for each gamma rate) of a single column full of gaps
    std::vector<double> pc0;
    pc0.resize(num_gamma_categories);

//    tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
//    tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();
//
//    int sonLeftID = treemap_.right.at(vnode_left);
//    bpp::Node *sonLeft = tree_->getNode(sonLeftID);
//
//    int sonRightID = treemap_.right.at(vnode_right);
//    bpp::Node *sonRight = tree_->getNode(sonRightID);

    int idx;

    for (int catg = 0; catg < num_gamma_categories; catg++) {

        // computes the recursive Felsenstein's peeling weight on the left subtree
//        idx = 0;
//        bpp::ColMatrix<double> fvL = computeFVrec(sonLeft, sL, idx, catg);
//
//        // computes the recursive Felsenstein's peeling weight on the right subtree
//        idx = 0;
//        bpp::ColMatrix<double> fvR = computeFVrec(sonRight, sR, idx, catg);

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(prNode_[sonLeftID].at(catg), fvL.at(catg), PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(prNode_[sonRightID].at(catg), fvR.at(catg), PrfvR);

        // fv = PrfvL * PrfvR
        bpp::ColMatrix<double> fv;
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

        Fv_gap.at(catg) = fv;

        // fv0 = pi * fv
        double fv0 = MatrixBppUtils::dotProd(fv, pi_);

        // lk at the actual node (considered as root node => beta = 1.0)
        double p0 = iotasNode_[nodeID][catg] * fv0;

        double pL,pR;
////        if (flag_RAM) {
//            if(sonLeft->isLeaf()){
//                pL = compute_lk_gap_down(sonLeft, sL, catg);
//            }else{
        pL=lk_empty_down_[sonLeftID][catg];
        //}
//            if(sonRight->isLeaf()){
//                pR = compute_lk_gap_down(sonRight, sR, catg);
//            }else{
        pR=lk_empty_down_[sonRightID][catg];
        //}
//        }else{
//            pL = compute_lk_gap_down(sonLeft, sL, catg);
//            pR = compute_lk_gap_down(sonRight, sR, catg);
//        }

        pc0.at(catg) = p0 + pL + pR;
    }

    return pc0;
}

std::vector<double> pPIP::computeLK_GapColumn_local(bpp::Node *node,
                                                    MSAcolumn_t &sL,
                                                    MSAcolumn_t &sR,
                                                    bool flag_RAM) {

    // number of discrete gamma categories
    int num_gamma_categories = rDist_->getNumberOfCategories();

    // array of lk (for each gamma rate) of a single column full of gaps
    std::vector<double> pc0;
    pc0.resize(num_gamma_categories);

    tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
    tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

    int sonLeftID = treemap_.right.at(vnode_left);
    bpp::Node *sonLeft = tree_->getNode(sonLeftID);

    int sonRightID = treemap_.right.at(vnode_right);
    bpp::Node *sonRight = tree_->getNode(sonRightID);

    for (int catg = 0; catg < num_gamma_categories; catg++) {

        int idx;

        // computes the recursive Felsenstein's peeling weight on the left subtree
        idx = 0;
        bpp::ColMatrix<double> fvL = computeFVrec(sonLeft, sL, idx, catg);

        // computes the recursive Felsenstein's peeling weight on the right subtree
        idx = 0;
        bpp::ColMatrix<double> fvR = computeFVrec(sonRight, sR, idx, catg);

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(prNode_[sonLeftID].at(catg), fvL, PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(prNode_[sonRightID].at(catg), fvR, PrfvR);

        // fv = PrfvL * PrfvR
        bpp::ColMatrix<double> fv;
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

        // fv0 = pi * fv
        double fv0 = MatrixBppUtils::dotProd(fv, pi_);

        // lk at the actual node (considered as root node => beta = 1.0)
        double p0 = iotasNode_[node->getId()][catg] * fv0;

        double pL,pR;
        if (flag_RAM) {
            if(sonLeft->isLeaf()){
                pL = compute_lk_gap_down(sonLeft, sL, catg);
            }else{
                pL=lk_empty_down_[sonLeftID][catg];
            }
            if(sonRight->isLeaf()){
                pR = compute_lk_gap_down(sonRight, sR, catg);
            }else{
                pR=lk_empty_down_[sonRightID][catg];
            }
        }else{
            pL = compute_lk_gap_down(sonLeft, sL, catg);
            pR = compute_lk_gap_down(sonRight, sR, catg);
        }

        pc0.at(catg) = p0 + pL + pR;
    }

    return pc0;
}

double pPIP::computeLK_M_local(int nodeID,
                               int sonLeftID,
                               int sonRightID,
                               MSAcolumn_t &sL,
                               MSAcolumn_t &sR,
                               std::vector< bpp::ColMatrix<double> > &fvL,
                               std::vector< bpp::ColMatrix<double> > &fvR,
                               std::vector< bpp::ColMatrix<double> > &Fv_M_ij) {

    int idx;

    // number of discrete gamma categories
    int num_gamma_categories = rDist_->getNumberOfCategories();

    double pr = 0.0;
    for (int catg = 0; catg < num_gamma_categories; catg++) {

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(prNode_[sonLeftID].at(catg), fvL.at(catg), PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(prNode_[sonRightID].at(catg), fvR.at(catg), PrfvR);

        // fv = PrfvL * PrfvR
        bpp::ColMatrix<double> fv;
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

        Fv_M_ij[catg] = fv;

        // fv0 = pi * fv
        double fv0 = MatrixBppUtils::dotProd(fv, pi_);

        // match probability with gamma
        double p = rDist_->getProbability((size_t) catg) * \
               iotasNode_[nodeID][catg] * \
               betasNode_[nodeID][catg] * \
               fv0;

        // marginal lk, marginalized over N gamma discrete classes
        pr += p;
    }

    return log((long double) pr);
}

double pPIP::computeLK_MXY_local(double NU,
                                 double valM,
                                 double valX,
                                 double valY,
                                 double log_pr,
                                 int m) {

    return log(NU) - log((double) m) + log_pr + max_of_three(valM, valX, valY, DBL_EPSILON,true);
}

double pPIP::computeLK_M_local(double NU,
                               double valM,
                               double valX,
                               double valY,
                               bpp::Node *node,
                               MSAcolumn_t &sL,
                               MSAcolumn_t &sR,
                               int m,
                               std::map<MSAcolumn_t, double> &lkM,
                               std::vector< std::vector<double> > &lkM_pattern,
                               bool flag_map,
                               bool flag_pattern,
                               bool flag_RAM) {

    double log_pr;

    tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
    tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

    int sonLeftID = treemap_.right.at(vnode_left);
    bpp::Node *sonLeft = tree_->getNode(sonLeftID);

    int sonRightID = treemap_.right.at(vnode_right);
    bpp::Node *sonRight = tree_->getNode(sonRightID);

    int nodeID = node->getId();

    // create left + right column
    MSAcolumn_t s;
    s.append(sL);
    s.append(sR);

    // If the columns has been already seen in the alignment, then do not compute the value, but get it from the map
    bool is_found;
    std::map<MSAcolumn_t, double>::iterator it;

    if(flag_pattern){

    }else {
        if (flag_map) {
            it = lkM.find(s);
            is_found = (it == lkM.end());
        } else {
            is_found = true;
        }
    }

    if(flag_pattern){
        //log_pr=lkM_pattern[][];
    }else {
        if (is_found) {
            // is the first time that it computes the lk of this column

            int idx;

            // number of discrete gamma categories
            int num_gamma_categories = rDist_->getNumberOfCategories();

            double pr = 0.0;
            for (int catg = 0; catg < num_gamma_categories; catg++) {

                // computes the recursive Felsenstein's peeling weight on the left subtree
                idx = 0;
                bpp::ColMatrix<double> fvL = computeFVrec(sonLeft, sL, idx, catg);

                // computes the recursive Felsenstein's peeling weight on the right subtree
                idx = 0;
                bpp::ColMatrix<double> fvR = computeFVrec(sonRight, sR, idx, catg);

                // PrfvL = Pr_L * fv_L
                bpp::ColMatrix<double> PrfvL;
                bpp::MatrixTools::mult(prNode_[sonLeftID].at(catg), fvL, PrfvL);

                // PrfvR = Pr_R * fv_R
                bpp::ColMatrix<double> PrfvR;
                bpp::MatrixTools::mult(prNode_[sonRightID].at(catg), fvR, PrfvR);

                // fv = PrfvL * PrfvR
                bpp::ColMatrix<double> fv;
                bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

                // fv0 = pi * fv
                double fv0 = MatrixBppUtils::dotProd(fv, pi_);

                // match probability with gamma
                double p = rDist_->getProbability((size_t) catg) * \
                       iotasNode_[nodeID][catg] * \
                       betasNode_[nodeID][catg] * \
                       fv0;

                // marginal lk, marginalized over N gamma discrete classes
                pr += p;
            }

            log_pr = log((long double) pr);

            lkM[s] = log_pr;

        } else {
            // the lk of a column like this has been already computed,
            // the lk value can be retrieved from the map

            log_pr = it->second;
        }
    }

    //return  log_phi_gamma + log_pr + max_of_three(valM, valX, valY, DBL_EPSILON);
    return log(NU) - log((double) m) + log_pr + max_of_three(valM, valX, valY, DBL_EPSILON,flag_RAM);
}

double pPIP::computeLK_X_local(int nodeID,
                               int sonLeftID,
                               int sonRightID,
                               MSAcolumn_t &sL,
                               MSAcolumn_t &col_gap_R,
                               std::vector< bpp::ColMatrix<double> > &fvL,
                               std::vector< bpp::ColMatrix<double> > &fvR,
                               std::vector< bpp::ColMatrix<double> > &Fv_X_ij) {

    // number of discrete gamma categories
    int num_gamma_categories = rDist_->getNumberOfCategories();

    double pL;
    double pr = 0.0;
    for (int catg = 0; catg < num_gamma_categories; catg++) {

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(prNode_[sonLeftID].at(catg), fvL.at(catg), PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(prNode_[sonRightID].at(catg), fvR.at(catg), PrfvR);

        // fv = PrfvL * PrfvR
        bpp::ColMatrix<double> fv;
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

        Fv_X_ij[catg] = fv;

        // fv0 = pi * fv
        double fv0 = MatrixBppUtils::dotProd(fv, pi_);

        // gapX probability with gamma
        double p0 = rDist_->getProbability((size_t) catg) * \
                iotasNode_[nodeID][catg] * \
                betasNode_[nodeID][catg] * \
                fv0;

//        if (sonLeft->isLeaf()) {
//            pL += compute_lk_down(sonLeft, sL, catg);
//        } else {
//            pL = lk_down_.at(sonLeftID).at(idx);
//            pL = exp(pL);
//        }

        pr += p0;
    }

//    pL = lk_down_.at(sonLeftID).at(idx);
//    pL = exp(pL);
//    pr += pL;

    return log((long double) pr);

}

double pPIP::computeLK_X_local(double NU,
                               double valM,
                               double valX,
                               double valY,
                               bpp::Node *node,
                               MSAcolumn_t &sL,
                               MSAcolumn_t &col_gap_R,
                               int m,
                               std::map<MSAcolumn_t, double> &lkX,
                               std::vector< std::vector<double> > &lkX_pattern,
                               bool flag_map,
                               bool flag_RAM,
                               bool flag_pattern) {

    double log_pr;

    tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
    tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

    int sonLeftID = treemap_.right.at(vnode_left);
    bpp::Node *sonLeft = tree_->getNode(sonLeftID);

    int sonRightID = treemap_.right.at(vnode_right);
    bpp::Node *sonRight = tree_->getNode(sonRightID);

    int nodeID = node->getId();

    // create left + right column
    MSAcolumn_t s;
    bool is_found;
    std::map<MSAcolumn_t, double>::iterator it;

    if(flag_pattern){

    }else {
        if (flag_map) {
            s.append(sL);
            s.append(col_gap_R);

            it = lkX.find(s);
            is_found = (it != lkX.end());
        } else {
            is_found = false;
        }
    }

    if(flag_pattern) {
        //log_pr=lkX_pattern[][];
    }else{
        if (!is_found) {
            // is the first time that it computes the lk of this column

            int idx;

            // number of discrete gamma categories
            int num_gamma_categories = rDist_->getNumberOfCategories();

            double pL;
            double pr = 0.0;
            for (int catg = 0; catg < num_gamma_categories; catg++) {

                // computes the recursive Felsenstein's peeling weight on the left subtree
                idx = 0;
                bpp::ColMatrix<double> fvL = computeFVrec(sonLeft, sL, idx, catg);

                // computes the recursive Felsenstein's peeling weight on the right subtree
                idx = 0;
                bpp::ColMatrix<double> fvR = computeFVrec(sonRight, col_gap_R, idx, catg);

                // PrfvL = Pr_L * fv_L
                bpp::ColMatrix<double> PrfvL;
                bpp::MatrixTools::mult(prNode_[sonLeftID].at(catg), fvL, PrfvL);

                // PrfvR = Pr_R * fv_R
                bpp::ColMatrix<double> PrfvR;
                bpp::MatrixTools::mult(prNode_[sonRightID].at(catg), fvR, PrfvR);

                // fv = PrfvL * PrfvR
                bpp::ColMatrix<double> fv;
                bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

                // fv0 = pi * fv
                double fv0 = MatrixBppUtils::dotProd(fv, pi_);

                // gapX probability with gamma
                double p0 = rDist_->getProbability((size_t) catg) * \
                        iotasNode_[nodeID][catg] * \
                        betasNode_[nodeID][catg] * \
                        fv0;

                if (flag_RAM) {
                    if (sonLeft->isLeaf()) {
                        pL += compute_lk_down(sonLeft, sL, catg);
                    } else {
                        pL = lk_down_.at(sonLeftID).at(idx);
                        pL = exp(pL);
                    }
                } else {
                    pL = compute_lk_down(sonLeft, sL, catg);
                }

                pr += p0 + pL;
            }

            log_pr = log((long double) pr);

            if (flag_map) {
                lkX[s] = log_pr;
            }

        } else {
            // the lk of a column like this has been already computed,
            // the lk value can be retrieved from the map

            log_pr = it->second;
        }
    }
    //return log_phi_gamma + log_pr + max_of_three(valM, valX, valY, DBL_EPSILON);
    return log(NU) - log((double) m) + log_pr + max_of_three(valM, valX, valY, DBL_EPSILON,flag_RAM);
}

double pPIP::computeLK_Y_local(int nodeID,
                               int sonLeftID,
                               int sonRightID,
                               MSAcolumn_t &col_gap_L,
                               MSAcolumn_t &sR,
                               std::vector< bpp::ColMatrix<double> > &fvL,
                               std::vector< bpp::ColMatrix<double> > &fvR,
                               std::vector< bpp::ColMatrix<double> > &Fv_Y_ij) {

    // number of discrete gamma categories
    int num_gamma_categories = rDist_->getNumberOfCategories();

    int idx;

    double pR;
    double pr = 0.0;
    for (int catg = 0; catg < num_gamma_categories; catg++) {

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(prNode_[sonLeftID].at(catg), fvL.at(catg), PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(prNode_[sonRightID].at(catg), fvR.at(catg), PrfvR);

        // fv = PrfvL * PrfvR
        bpp::ColMatrix<double> fv;
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

        Fv_Y_ij[catg] = fv;

        // fv0 = pi * fv
        double fv0 = MatrixBppUtils::dotProd(fv, pi_);

        // gapY probability with gamma
        double p0 = rDist_->getProbability((size_t) catg) * \
                    iotasNode_[nodeID][catg] * \
                    betasNode_[nodeID][catg] * \
                    fv0;

        pr += p0;

    }

//    pR = lk_down_.at(sonRightID).at(idx);
//    pR = exp(pR);
//    pr += pR;

    return log((long double) pr);

}

double pPIP::computeLK_Y_local(double NU,
                               double valM,
                               double valX,
                               double valY,
                               bpp::Node *node,
                               MSAcolumn_t &col_gap_L,
                               MSAcolumn_t &sR,
                               int m,
                               std::map<MSAcolumn_t, double> &lkY,
                               std::vector< std::vector<double> > &lkY_pattern,
                               bool flag_map,
                               bool flag_RAM,
                               bool flag_pattern) {

    double log_pr;

    tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
    tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

    int sonLeftID = treemap_.right.at(vnode_left);
    bpp::Node *sonLeft = tree_->getNode(sonLeftID);

    int sonRightID = treemap_.right.at(vnode_right);
    bpp::Node *sonRight = tree_->getNode(sonRightID);

    int nodeID = node->getId();

    // create left + right column
    MSAcolumn_t s;
    bool is_found;
    std::map<MSAcolumn_t, double>::iterator it;

    if(flag_pattern){

    }else{
        if(flag_map){
            s.append(col_gap_L);
            s.append(sR);

            it = lkY.find(s);
            is_found= (it != lkY.end());
        }else{
            is_found=false;
        }
    }

    if(flag_pattern){
        //log_pr=lkY_pattern[][];
    }else{

        if (!is_found) {
            // is the first time that it computes the lk of this column

            // number of discrete gamma categories
            int num_gamma_categories = rDist_->getNumberOfCategories();

            int idx;

            double pR;
            double pr = 0.0;
            for (int catg = 0; catg < num_gamma_categories; catg++) {

                // computes the recursive Felsenstein's peeling weight on the left subtree
                idx = 0;
                bpp::ColMatrix<double> fvL = computeFVrec(sonLeft, col_gap_L, idx, catg);

                // computes the recursive Felsenstein's peeling weight on the right subtree
                idx = 0;
                bpp::ColMatrix<double> fvR = computeFVrec(sonRight, sR, idx, catg);

                // PrfvL = Pr_L * fv_L
                bpp::ColMatrix<double> PrfvL;
                bpp::MatrixTools::mult(prNode_[sonLeftID].at(catg), fvL, PrfvL);

                // PrfvR = Pr_R * fv_R
                bpp::ColMatrix<double> PrfvR;
                bpp::MatrixTools::mult(prNode_[sonRightID].at(catg), fvR, PrfvR);

                // fv = PrfvL * PrfvR
                bpp::ColMatrix<double> fv;
                bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

                // fv0 = pi * fv
                double fv0 = MatrixBppUtils::dotProd(fv, pi_);

                // gapY probability with gamma
                double p0 = rDist_->getProbability((size_t) catg) * \
                            iotasNode_[nodeID][catg] * \
                            betasNode_[nodeID][catg] * \
                            fv0;


                if (flag_RAM) {
                    if(sonRight->isLeaf()){
                        pR = compute_lk_down(sonRight, sR, catg);
                    }else{
                        pR = lk_down_.at(sonRightID).at(idx);
                        pR = exp(pR);
                    }
                } else {
                    pR = compute_lk_down(sonRight, sR, catg);
                }

                pr += p0 + pR;

            }

            log_pr = log((long double) pr);

            if(flag_map) {
                lkY[s] = log_pr;
            }

        } else {
            // the lk of a column like this has been already computed,
            // the lk value can be retrieved from the map

            log_pr = it->second;
        }

    }

    //return log_phi_gamma + log_pr + max_of_three(valM, valX, valY, DBL_EPSILON);
    return log(NU) - log((double) m) + log_pr + max_of_three(valM, valX, valY, DBL_EPSILON,flag_RAM);
}

void pPIP::DP3D_PIP_RAM(bpp::Node *node,
                        bool local,
                        bool flag_map,
                        bool flag_pattern) {

    // TODO: place as argument
    // used to select random when 2 or 3 lks (M,X,Y) have "exactly" the same value
    //bool randomSeed = true;

    // number of discrete gamma categories
    size_t num_gamma_categories = rDist_->getNumberOfCategories();

    if (local) {
        // recompute local tau, total tree length of a tree root at the given node
        _setTau(treemap_.left.at(node->getId()));

        // recompute the normalizing factor nu for the local tree
        _setNu();

        // recompute lambdas with the new normalizing factor (local tree), flag true = tree rooted here
        _setAllIotas(node, true);

        // recompute betas with the new normalizing factor (local tree), flag true = tree rooted here
        _setAllBetas(node, true);
    }

    int h, w;

    // Get the IDs of the sons nodes given the current node
    int nodeID = node->getId();

    tshlib::VirtualNode *vnode_left = treemap_.left.at(nodeID)->getNodeLeft(); // bpp::Node to tshlib::VirtualNode
    tshlib::VirtualNode *vnode_right = treemap_.left.at(nodeID)->getNodeRight(); // bpp::Node to tshlib::VirtualNode
    int nodeID_L = treemap_.right.at(vnode_left);
    int nodeID_R = treemap_.right.at(vnode_right);

    /*
    // Compute dimensions of the 3D block at current internal node.
    h = MSA_.at(nodeID_L).size() + 1; // dimension of the alignment on the left side
    w = MSA_.at(nodeID_R).size() + 1; // dimension of the alignment on the right side

    int d = (h - 1) + (w - 1) + 1; // third dimension of the DP matrix
    */

    h = MSA_.at(nodeID_L).size(); // dimension of the alignment on the left side
    w = MSA_.at(nodeID_R).size(); // dimension of the alignment on the right side

    int d = h + w + 1; // third dimension of the DP matrix

    // lk of a single empty column (full of gaps) with rate variation (gamma distribution)
    std::vector<double> pc0;

    // MSA columns
    MSAcolumn_t sLs; // left column
    MSAcolumn_t sRs; // right column
    MSAcolumn_t col_gap_Ls; // left column full of gaps
    MSAcolumn_t col_gap_Rs; //right column full of gaps

    int numLeavesLeft = seqNames_.at(nodeID_L).size(); // number of leaves in the left sub-tree
    int numLeavesRight = seqNames_.at(nodeID_R).size(); // number of leaves in the right sub-tree

    col_gap_Ls = createGapCol(numLeavesLeft); // create column of gaps for the left sub-tree
    col_gap_Rs = createGapCol(numLeavesRight); // create column of gaps for the right sub-tree

    signed long seed;
//    if (randomSeed) {
//        seed = std::chrono::system_clock::now().time_since_epoch().count(); // "random" seed
//    } else {
//        seed = 0; // fixed seed
//    }

    std::default_random_engine generator(seed_);              // jatiapp seed
    std::uniform_real_distribution<double> distribution(0.0,
                                                        1.0); // Uniform distribution for the selection of lks with the same value

    auto epsilon = DBL_EPSILON;

    std::vector< std::vector<double> > lkM_pattern;
    std::vector< std::vector<double> > lkX_pattern;
    std::vector< std::vector<double> > lkY_pattern;

    int pattern_size_L;
    int pattern_size_R;

    if(flag_pattern) {
        lkM_pattern.resize(pattern_size_L);
        lkX_pattern.resize(pattern_size_L);
        lkY_pattern.resize(pattern_size_L);

        for (int k = 0; k < pattern_size_L; k++) {
            lkM_pattern[k].resize(pattern_size_R, -std::numeric_limits<double>::infinity());
            lkX_pattern[k].resize(pattern_size_R, -std::numeric_limits<double>::infinity());
            lkY_pattern[k].resize(pattern_size_R, -std::numeric_limits<double>::infinity());
        }
    }



    //***************************************************************************************
    //***************************************************************************************
    if (local) {
        // compute the lk of a column full of gaps
        //pc0 = computeLK_GapColumn_local(node, col_gap_Ls, col_gap_Rs,true);
        pc0 = computeLK_GapColumn_local(nodeID,
                                        nodeID_L,
                                        nodeID_R,
                                        fv_empty_data_[nodeID_L],
                                        fv_empty_data_[nodeID_R],
                                        fv_empty_data_[nodeID]);
        lk_empty_down_[nodeID]=pc0;
    } else {
        //global
    }

    //***************************************************************************************
    // Initialisation of the datastructure
    //***************************************************************************************
    std::vector< vector< vector<double> > > LogM;   // DP sparse matrix for MATCH case (only 2 layer are needed)
    std::vector< vector< vector<double> > > LogX;   // DP sparse matrix for GAPX case (only 2 layer are needed)
    std::vector< vector< vector<double> > > LogY;   // DP sparse matrix for GAPY case (only 2 layer are needed)
    std::vector< vector< vector<int> > > TR;        // 3D traceback matrix
    std::vector< vector< vector<double> > > LK;     // 3D LK matrix, stores best lk at each position

    LogM.resize(2);
    LogX.resize(2);
    LogY.resize(2);
    TR.resize(d);
    LK.resize(d);

    /*
    auto ***LogM = new double **[2];
    auto ***LogX = new double **[2];
    auto ***LogY = new double **[2];
    auto ***TR = new int **[d];
    auto ***LK = new double **[d];
    */

    // allocate memory for the 2 layers
    for(int k = 0; k < 2; k++){
        //LogM[k] = new double *[h]();
        //LogX[k] = new double *[h]();
        //LogY[k] = new double *[h]();

        LogM[k].resize(h);
        LogX[k].resize(h);
        LogY[k].resize(h);

        for(int i = 0; i < h; i++){
            //LogM[k][i] = new double[w]();
            //LogX[k][i] = new double[w]();
            //LogY[k][i] = new double[w]();

            LogM[k][i].resize(w,0);
            LogX[k][i].resize(w,0);
            LogY[k][i].resize(w,0);
        }
    }

//    TR[0] = new int*[1]();
//    TR[0][0] = new int[1]();

    TR[0].resize(1);
    LK[0].resize(1);

    TR[0][0].resize(1,0);
    LK[0][0].resize(1,0);

//    LK[0] = new double*[1]();
//    LK[0][0] = new double[1]();
//    for(int k = 1; k < d; k++){
////        TR[k] = new int *[h]();
////        LK[k] = new double *[h]();
//        TR[k].resize(h);
//        LK[k].resize(h);
//        for(int i = 0; i < h; i++){
////            TR[k][i] = new int[w]();
////            LK[k][i] = new double[w]();
//            TR[k][i].resize(w,0);
//            LK[k][i].resize(w,0);
//        }
//    }

    //***************************************************************************************
    //***************************************************************************************
    // marginal likelihood for all empty columns with rate variation (gamma distribution)
    // phi(m,pc0,r) depends on the MSA length m
    // marginal phi marginalized over gamma categories
    double log_phi_gamma;

    auto **PHI = new double *[d];
    double PC0 = 0.0;
    double NU = 0.0;

    for (int i = 0; i < d; i++) {
        PHI[i] = new double[num_gamma_categories];
    }

    for (int catg = 0; catg < num_gamma_categories; catg++) {
        // log( P_gamma(r) * phi(0,pc0(r),r) ): marginal lk for all empty columns of an alignment of size 0
        // PHI[0][catg] = log(rDist_->getProbability((size_t)catg)) + (nu_.at(catg) * (pc0.at(catg) - 1.0));
        PC0 += rDist_->getProbability((size_t) catg) * pc0.at(catg);
        NU += rDist_->getProbability((size_t) catg) * nu_.at(catg);
    }

    // computes the marginal phi marginalized over all the gamma categories
    log_phi_gamma = NU * (PC0 - 1);
    //============================================================

    LogM[0][0][0] = log_phi_gamma;
    LogX[0][0][0] = log_phi_gamma;
    LogY[0][0][0] = log_phi_gamma;
    TR[0][0][0] = STOP_STATE;
    LK[0][0][0] = -std::numeric_limits<double>::infinity();

    double max_of_MXY = -std::numeric_limits<double>::infinity();                   // Max value found in the matrices M,X,Y


    int m_binary_this;            // Level Index during computation / current
    int m_binary_prev;            // Level Index during computation / old

    double valM_this,valM_prev;
    double valX_this,valX_prev;
    double valY_this,valY_prev;

    int tr; // traceback index
    double val_tr;

    max_val_str max_val_index;

    double score = -std::numeric_limits<double>::infinity();        // best likelihood value at this node

    signed long level_max_lk = INT_MIN;                             // Depth in M,X,Y with the highest lk value
    //int depth;

    //int last_d = d - 1;
    std::map<MSAcolumn_t, double> lkM;
    std::map<MSAcolumn_t, double> lkX;
    std::map<MSAcolumn_t, double> lkY;

    //============================================================
    // early stop condition flag
    bool flag_exit = false;
    int counter_to_early_stop = 0;
    int max_decrease_before_stop = 10;                              // hardcoded to prevent early-stops
    double prev_best_lk = -std::numeric_limits<double>::infinity(); // TODO: Check for the same bug in the other version
    // ============================================================
    // For each slice of the 3D cube, compute the values of each cell
    for (int m = 1; m < d; m++) {

        if (flag_exit) {
            break;
        }

        // alternate the two layers
        m_binary_this = m % 2;
        m_binary_prev = (m + 1) % 2;

        //===========================
        TR[m].resize(h);
        LK[m].resize(h);
        for(int z = 0; z < h; z++){
            TR[m][z].resize(w,0);
            LK[m][z].resize(w,0);
        }
        //===========================

        //***************************************************************************************
        //***************************************************************************************
        for (int catg = 0; catg < num_gamma_categories; catg++) {
            // computes the marginal phi(m,pc0(r),r) with gamma by multiplying the starting value
            // phi(0,pco(r),r) = log( P_gamma(r) * exp( nu(r) * (pc0(r)-1) ) ) with
            // 1/m * nu(r) at each new layer
            PHI[m][catg] = PHI[m - 1][catg] - log((long double) m) + log((long double) nu_.at(catg));
        }

        // computes the marginal phi marginalized over all the gamma categories
        log_phi_gamma = PHI[m][0];
        for (int catg = 1; catg < num_gamma_categories; catg++) {
            log_phi_gamma = pPIPUtils::add_lns(log_phi_gamma, PHI[m][catg]);
        }
        //***************************************************************************************
        //***************************************************************************************
        for (int i = 0; i < h; i++) {

            sLs = (MSA_.at(nodeID_L).at(i));

            for (int j = 0; j < w; j++) {

                sRs = (MSA_.at(nodeID_R).at(j));

                //=======================================
                // MATCH
                if( (i-1) >= 0 && (j-1) >= 0 ){

                    valM_prev = LogM[m_binary_prev][i-1][j-1];
                    valX_prev = LogX[m_binary_prev][i-1][j-1];
                    valY_prev = LogY[m_binary_prev][i-1][j-1];

                    if (std::isinf(valM_prev) && std::isinf(valX_prev) && std::isinf(valY_prev)) {

                        LogM[m_binary_this][i][j] = -std::numeric_limits<double>::infinity();

                    }else{

                        if (local) {

                            valM_this = computeLK_M_local(NU,
                                                          valM_prev,
                                                          valX_prev,
                                                          valY_prev,
                                                          node,
                                                          sLs,
                                                          sRs,
                                                          m,
                                                          lkM,
                                                          lkM_pattern,
                                                          flag_map,
                                                          true,
                                                          flag_pattern);

//                            if (std::isinf(valM_this)) {
//                                LOG(FATAL)
//                                        << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is infinite. Check call stack below.";
//                            }
//
                            if (std::isnan(valM_this)) {
                                LOG(FATAL)
                                        << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is nan. Check call stack below.";
                            }

                            LogM[m_binary_this][i][j] = valM_this;

                        } else {
                            //global
                        }

                    }

                }
                //=======================================
                // GAPX
                if( (i-1) >=0 ){

                    valM_prev = LogM[m_binary_prev][i-1][j];
                    valX_prev = LogX[m_binary_prev][i-1][j];
                    valY_prev = LogY[m_binary_prev][i-1][j];

                    if (std::isinf(valM_prev) && std::isinf(valX_prev) && std::isinf(valY_prev)) {

                        LogX[m_binary_this][i][j] = -std::numeric_limits<double>::infinity();

                    }else{

                        if (local) {

                            valX_this = computeLK_X_local(NU,
                                                          valM_prev,
                                                          valX_prev,
                                                          valY_prev,
                                                          node,
                                                          sLs,
                                                          col_gap_Rs,
                                                          m,
                                                          lkX,
                                                          lkX_pattern,
                                                          flag_map,
                                                          true,
                                                          flag_pattern);

//                            if (std::isinf(valX_this)) {
//                                LOG(FATAL)
//                                        << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is infinite. Check call stack below.";
//                            }
//
                            if (std::isnan(valX_this)) {
                                LOG(FATAL)
                                        << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is nan. Check call stack below.";
                            }

                            LogX[m_binary_this][i][j] = valX_this;

                        } else {
                            //global
                        }

                    }

                }
                //=======================================
                // GAPY
                if( (j-1)>=0 ){

                    valM_prev = LogM[m_binary_prev][i][j-1];
                    valX_prev = LogX[m_binary_prev][i][j-1];
                    valY_prev = LogY[m_binary_prev][i][j-1];

                    if (std::isinf(valM_prev) && std::isinf(valX_prev) && std::isinf(valY_prev)) {

                        LogY[m_binary_this][i][j] = -std::numeric_limits<double>::infinity();

                    }else{

                        if (local) {

                            valY_this = computeLK_Y_local(NU,
                                                          valM_prev,
                                                          valX_prev,
                                                          valY_prev,
                                                          node,
                                                          col_gap_Ls,
                                                          sRs,
                                                          m,
                                                          lkY,
                                                          lkY_pattern,
                                                          flag_map,
                                                          true,
                                                          flag_pattern);

//                            if (std::isinf(valY_this)) {
//                                LOG(FATAL)
//                                        << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is infinite. Check call stack below.";
//                            }
//
                            if (std::isnan(valY_this)) {
                                LOG(FATAL)
                                        << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is nan. Check call stack below.";
                            }

                            LogY[m_binary_this][i][j] = valY_this;

                        } else {
                            //global
                        }

                    }

                }
                //=======================================
                valM_this = LogM[m_binary_this][i][j];
                valX_this = LogX[m_binary_this][i][j];
                valY_this = LogY[m_binary_this][i][j];

                //valM_this = fabs((long double) valM_this) < epsilon ? -std::numeric_limits<double>::infinity() : valM_this;
                //valX_this = fabs((long double) valX_this) < epsilon ? -std::numeric_limits<double>::infinity() : valX_this;
                //valY_this = fabs((long double) valY_this) < epsilon ? -std::numeric_limits<double>::infinity() : valY_this;

                // Find which matrix contains the best value of LK found until this point.
                index_of_max(valM_this, valX_this, valY_this, epsilon, generator, distribution,max_val_index,true,tr,val_tr);

                switch (max_val_index.index) {
                    case MATCH_STATE:
                        LK[m][i][j]=valM_this;
                        break;
                    case GAP_X_STATE:
                        LK[m][i][j]=valX_this;
                        break;
                    case GAP_Y_STATE:
                        LK[m][i][j]=valY_this;
                        break;
                    default:
                        LOG(FATAL) <<"\nSomething went wrong in reading the TR value. TR is neither MATCH, nor GAPX, nor GAPY. ";
                }
                /*
                if (TR[m][i][j] != 0) {
                    LOG_IF(FATAL, TR[m][i][j] != 0 ) << "\nSomething went wrong in accessing TR at indices:[" << m << "][" << i << "][" << j << "] in function pPIP::DP3D_PIP. Check call stack below.";
                }
*/
                // Store the index for the traceback
                TR[m][i][j] = max_val_index.index;

                // If we reached the corner of the 3D cube, then:
                if ( (i == (h - 1)) & (j == (w - 1))) {
                    // the algorithm is filling the last column of 3D DP matrix where
                    // all the characters are in the MSA

                    //max_of_MXY = max_of_three(valM_this, valX_this, valY_this, epsilon,true);

                    max_of_MXY = max_val_index.val;


                    if (max_of_MXY > score) {
                        score = max_of_MXY;
                        level_max_lk = m;
                    }

                    //=====================================================================
                    // early stop condition
                    if (score < prev_best_lk) {
                        prev_best_lk = score;
                        counter_to_early_stop++;
                        if (counter_to_early_stop > max_decrease_before_stop) {
                            // if for max_decrease_before_stop consecutive times
                            // the lk decrease then exit, the maximum lk has been reached
                            flag_exit = true;
                        }
                    } else {
                        counter_to_early_stop = 0;
                    }
                    //=====================================================================

                }
            }

        }

    }

    // level (k position) in the DP matrix that contains the highest lk value
    //    depth     = level_max_lk;

    score_.at(nodeID) = score;

    lk_down_[nodeID].resize(level_max_lk);

    //==========================================================================================
    // start backtracing the 3 matrices (MATCH, GAPX, GAPY)
    TracebackPath_t traceback_path(level_max_lk, ' ');
    int id1 = h - 1;
    int id2 = w - 1;
    for (int lev = level_max_lk; lev > 0; lev--) {
        int state = TR[lev][id1][id2];
        switch (state) {
            case MATCH_STATE:

                lk_down_.at(nodeID).at(lev - 1)=LK[lev][id1][id2];

                id1 = id1 - 1;
                id2 = id2 - 1;

                traceback_path[lev - 1] = MATCH_CHAR;

                break;
            case GAP_X_STATE:

                lk_down_.at(nodeID).at(lev - 1)=LK[lev][id1][id2];

                id1 = id1 - 1;

                traceback_path[lev - 1] = GAP_X_CHAR;

                break;
            case GAP_Y_STATE:

                lk_down_.at(nodeID).at(lev - 1)=LK[lev][id1][id2];

                id2 = id2 - 1;

                traceback_path[lev - 1] = GAP_Y_CHAR;

                break;
            default:
                LOG(FATAL) << "\nSomething went wrong during the alignment reconstruction in function pPIP::DP3D_PIP. Check call stack below.";
        }
    }

    traceback_path_.at(nodeID) = traceback_path;

    // converts traceback path into an MSA
    build_MSA(node, traceback_path);


    //TODO: compress


    // assigns the sequence names of the new alligned sequences to the current MSA
    setMSAsequenceNames(node);
    //==========================================================================================

    //==========================================================================================
    // memory freeing
    //delete[] &LogM;
    //delete[] &LogX;
    //delete[] &LogY;
    //delete[] &LK;
    //delete[] &TR;
    //==========================================================================================
}
void pPIP::DP3D_PIP_RAM_FAST(bpp::Node *node) {

    std::cout<<"\n"<<node->getName()<<"\n\n";

    // four levels of optimization:
    // 1) to pre-compute DP 2D
    // 2) to use LK (lk_down)
    // 3) to re-use fv
    // 4) to use compressed data
    // 5) triangle


    // number of discrete gamma categories
    size_t num_gamma_categories = rDist_->getNumberOfCategories();

    // recompute local tau, total tree length of a tree root at the given node
    _setTau(treemap_.left.at(node->getId()));

    // recompute the normalizing factor nu for the local tree
    _setNu();

    // recompute lambdas with the new normalizing factor (local tree), flag true = tree rooted here
    _setAllIotas(node, true);

    // recompute betas with the new normalizing factor (local tree), flag true = tree rooted here
    _setAllBetas(node, true);

    // Get the IDs of the sons nodes given the current node
    int nodeID = node->getId();

    tshlib::VirtualNode *vnode_left = treemap_.left.at(nodeID)->getNodeLeft(); // bpp::Node to tshlib::VirtualNode
    tshlib::VirtualNode *vnode_right = treemap_.left.at(nodeID)->getNodeRight(); // bpp::Node to tshlib::VirtualNode
    int nodeID_L = treemap_.right.at(vnode_left);
    int nodeID_R = treemap_.right.at(vnode_right);

    //=======================================
    if(vnode_left->isTerminalNode()) {
        bpp::Node *sonLeft = tree_->getNode(nodeID_L);
        compressMSA(sonLeft);
        setFVleaf(sonLeft);
        set_lk_leaf(sonLeft);
        set_lk_empty_leaf(sonLeft);
    }
    if(vnode_right->isTerminalNode()) {
        bpp::Node *sonRight = tree_->getNode(nodeID_R);
        compressMSA(sonRight);
        setFVleaf(sonRight);
        set_lk_leaf(sonRight);
        set_lk_empty_leaf(sonRight);
    }
    //=======================================

    // Compute dimensions of the 3D block at current internal node.
    int h = MSA_.at(nodeID_L).size() + 1; // dimension of the alignment on the left side
    int w = MSA_.at(nodeID_R).size() + 1; // dimension of the alignment on the right side
    int d = (h - 1) + (w - 1) + 1; // third dimension of the DP matrix

    // lk of a single empty column (full of gaps) with rate variation (gamma distribution)
    std::vector<double> pc0;

    // MSA columns
    MSAcolumn_t sLs; // left column
    MSAcolumn_t sRs; // right column
    MSAcolumn_t col_gap_Ls; // left column full of gaps
    MSAcolumn_t col_gap_Rs; //right column full of gaps

    int numLeavesLeft = seqNames_.at(nodeID_L).size(); // number of leaves in the left sub-tree
    int numLeavesRight = seqNames_.at(nodeID_R).size(); // number of leaves in the right sub-tree

    col_gap_Ls = createGapCol(numLeavesLeft); // create column of gaps for the left sub-tree
    col_gap_Rs = createGapCol(numLeavesRight); // create column of gaps for the right sub-tree

    signed long seed;

    std::default_random_engine generator(seed_);              // jatiapp seed
    std::uniform_real_distribution<double> distribution(0.0,
                                                        1.0); // Uniform distribution for the selection of lks with the same value

    auto epsilon = DBL_EPSILON;

    //***************************************************************************************
    //***************************************************************************************
    // compute the lk of a column full of gaps
    pc0 = computeLK_GapColumn_local(nodeID,
                                    nodeID_L,
                                    nodeID_R,
                                    fv_empty_data_[nodeID_L],
                                    fv_empty_data_[nodeID_R],
                                    fv_empty_data_[nodeID]);
    lk_empty_down_[nodeID]=pc0;

    //***************************************************************************************
    // Initialisation of the data structure
    //***************************************************************************************
    std::vector< vector< vector<double> > > Log3DM;   // DP sparse matrix for MATCH case (only 2 layer are needed)
    std::vector< vector< vector<double> > > Log3DX;   // DP sparse matrix for GAPX case (only 2 layer are needed)
    std::vector< vector< vector<double> > > Log3DY;   // DP sparse matrix for GAPY case (only 2 layer are needed)
    std::vector< vector< vector<int> > > TR;        // 3D traceback matrix
    std::vector< vector< vector<double> > > LK3D;     // 3D LK matrix, stores best lk at each position

    std::vector< vector< vector<bool> > > bool_MXY;

    Log3DM.resize(2);
    Log3DX.resize(2);
    Log3DY.resize(2);

    bool_MXY.resize(2);

    TR.resize(d);
    LK3D.resize(d);

    // allocate memory for the 2 layers
    for(int k = 0; k < 2; k++){
        Log3DM[k].resize(h);
        Log3DX[k].resize(h);
        Log3DY[k].resize(h);

        bool_MXY[k].resize(h);

        for(int i = 0; i < h; i++){
            Log3DM[k][i].resize(w,-std::numeric_limits<double>::infinity());
            Log3DX[k][i].resize(w,-std::numeric_limits<double>::infinity());
            Log3DY[k][i].resize(w,-std::numeric_limits<double>::infinity());

            bool_MXY[k][i].resize(w,false);
        }
    }

    TR[0].resize(1);
    LK3D[0].resize(1);

    TR[0][0].resize(1,STOP_STATE);
    LK3D[0][0].resize(1,-std::numeric_limits<double>::infinity());

    //***************************************************************************************
    //***************************************************************************************
    // marginal likelihood for all empty columns with rate variation (gamma distribution)
    // phi(m,pc0,r) depends on the MSA length m
    // marginal phi marginalized over gamma categories
    double log_phi_gamma;

    auto **PHI = new double *[d];
    double PC0 = 0.0;
    double NU = 0.0;

    for (int i = 0; i < d; i++) {
        PHI[i] = new double[num_gamma_categories];
    }

    for (int catg = 0; catg < num_gamma_categories; catg++) {
        // log( P_gamma(r) * phi(0,pc0(r),r) ): marginal lk for all empty columns of an alignment of size 0
        // PHI[0][catg] = log(rDist_->getProbability((size_t)catg)) + (nu_.at(catg) * (pc0.at(catg) - 1.0));
        PC0 += rDist_->getProbability((size_t) catg) * pc0.at(catg);
        NU += rDist_->getProbability((size_t) catg) * nu_.at(catg);
    }

    // computes the marginal phi marginalized over all the gamma categories
    log_phi_gamma = NU * (PC0 - 1);
    //============================================================

    Log3DM[0][0][0] = log_phi_gamma;
    Log3DX[0][0][0] = log_phi_gamma;
    Log3DY[0][0][0] = log_phi_gamma;
    TR[0][0][0] = STOP_STATE;
    LK3D[0][0][0] = -std::numeric_limits<double>::infinity();

    //double max_of_MXY = -std::numeric_limits<double>::infinity();                   // Max value found in the matrices M,X,Y


    int m_binary_this;            // Level Index during computation / current
    int m_binary_prev;            // Level Index during computation / old

    double valM_this,valM_prev;
    double valX_this,valX_prev;
    double valY_this,valY_prev;

    max_val_str max_val_index;

    double score = -std::numeric_limits<double>::infinity();        // best likelihood value at this node

    signed long level_max_lk = INT_MIN;                             // Depth in M,X,Y with the highest lk value

    int tr;
    double max_of_3;

    int id1,id2;
    //============================================================
    // early stop condition flag
    bool flag_exit = false;
    int counter_to_early_stop = 0;
    int max_decrease_before_stop = 10;                              // hardcoded to prevent early-stops
    double prev_best_lk = -std::numeric_limits<double>::infinity(); // TODO: Check for the same bug in the other version
    //============================================================
    // first 2D-DP

    int h_compr = rev_map_compressed_seqs_.at(nodeID_L).size();
    int w_compr = rev_map_compressed_seqs_.at(nodeID_R).size();

    // memory allocation
    std::vector< vector<double> > Log2DM;
    std::vector<double> Log2DX;
    std::vector<double> Log2DY;

    std::vector< vector< vector< bpp::ColMatrix<double> > > > Fv_M;
    std::vector< vector< bpp::ColMatrix<double> > > Fv_X;
    std::vector< vector< bpp::ColMatrix<double> > > Fv_Y;

    Log2DM.resize(h_compr);
    Log2DX.resize(h_compr);
    Log2DY.resize(w_compr);

    Fv_M.resize(h_compr);
    Fv_X.resize(h_compr);
    Fv_Y.resize(w_compr);

    for(int i = 0; i < h_compr; i++){
        Log2DM[i].resize(w_compr);
        Fv_M[i].resize(w_compr);
        for(int j = 0; j < w_compr; j++){
            Fv_M[i][j].resize(num_gamma_categories);
        }
    }
    for(int i = 0; i < h_compr; i++){
        Fv_X[i].resize(num_gamma_categories);
    }
    for(int j = 0; j < w_compr; j++){
        Fv_Y[j].resize(num_gamma_categories);
    }
    //================================================================
    // MATCH
    for (int i = 0; i < h_compr; i++) {
        id1=rev_map_compressed_seqs_.at(nodeID_L).at(i);
        sLs = (MSA_.at(nodeID_L).at(id1));

        for (int j = 0; j < w_compr; j++) {
            id2=rev_map_compressed_seqs_.at(nodeID_R).at(j);
            sRs = (MSA_.at(nodeID_R).at(id2));

            Log2DM[i][j] = computeLK_M_local(nodeID,
                                             nodeID_L,
                                             nodeID_R,
                                             sLs,
                                             sRs,
                                             fv_data_.at(nodeID_L).at(i),
                                             fv_data_.at(nodeID_R).at(j),
                                             Fv_M[i][j]);

        }
    }
    //================================================================
    // GAPX
    for (int i = 0; i < h_compr; i++) {
        id1 = rev_map_compressed_seqs_.at(nodeID_L).at(i);
        sLs = (MSA_.at(nodeID_L).at(id1));

        Log2DX[i] = (computeLK_X_local(nodeID,
                                       nodeID_L,
                                       nodeID_R,
                                       sLs,
                                       col_gap_Rs,
                                       fv_data_[nodeID_L].at(i),
                                       fv_empty_data_[nodeID_R],
                                       Fv_X[i]) + \
                                      exp(lk_down_.at(nodeID_L).at(i)));

    }
    //================================================================
    // GAPY
    for (int j = 0; j < w_compr; j++) {
        id2 = rev_map_compressed_seqs_.at(nodeID_R).at(j);
        sRs = (MSA_.at(nodeID_R).at(id2));

        Log2DY[j] = (computeLK_Y_local(nodeID,
                                       nodeID_L,
                                       nodeID_R,
                                       col_gap_Ls,
                                       sRs,
                                       fv_empty_data_[nodeID_L],
                                       fv_data_[nodeID_R].at(j),
                                       Fv_Y[j]) + \
                                      exp(lk_down_.at(nodeID_R).at(j)));
    }
    //================================================================

    // For each slice of the 3D cube, compute the values of each cell

    bool_MXY[0][0][0]=true;
    for (int m = 1; m < d; m++) {

        //std::cout<<".....m="<<m<<std::endl;

        if (flag_exit) {
            break;
        }

        // alternate the two layers
        m_binary_this = m % 2;
        m_binary_prev = (m + 1) % 2;

        //===========================
        TR[m].resize(h);
        LK3D[m].resize(h);
        for(int z = 0; z < h; z++){
            TR[m][z].resize(w,0);
            LK3D[m][z].resize(w,0);
        }
        //===========================

        //***************************************************************************************
        //***************************************************************************************
        for (int catg = 0; catg < num_gamma_categories; catg++) {
            // computes the marginal phi(m,pc0(r),r) with gamma by multiplying the starting value
            // phi(0,pco(r),r) = log( P_gamma(r) * exp( nu(r) * (pc0(r)-1) ) ) with
            // 1/m * nu(r) at each new layer
            PHI[m][catg] = PHI[m - 1][catg] - log((long double) m) + log((long double) nu_.at(catg));
        }

        // computes the marginal phi marginalized over all the gamma categories
        log_phi_gamma = PHI[m][0];
        for (int catg = 1; catg < num_gamma_categories; catg++) {
            log_phi_gamma = pPIPUtils::add_lns(log_phi_gamma, PHI[m][catg]);
        }
        //***************************************************************************************
        //***************************************************************************************
        for (int i = 0; i < h; i++) {

            for (int j = 0; j < w; j++) {


                Log3DM[m_binary_this][i][j] = -std::numeric_limits<double>::infinity();
                Log3DX[m_binary_this][i][j] = -std::numeric_limits<double>::infinity();
                Log3DY[m_binary_this][i][j] = -std::numeric_limits<double>::infinity();

                bool_MXY[m_binary_this][i][j]=false;

                //=======================================
                // MATCH
                if( (i-1) >= 0 && (j-1) >= 0 ){

                    //id1 = rev_map_compressed_seqs_.at(nodeID_L).at(i-1);
                    //id2 = rev_map_compressed_seqs_.at(nodeID_R).at(j-1);
                    id1 = map_compressed_seqs_.at(nodeID_L).at(i-1);
                    id2 = map_compressed_seqs_.at(nodeID_R).at(j-1);

                    bool_MXY[m_binary_this][i][j]=true;

                    valM_prev = Log3DM[m_binary_prev][i-1][j-1];
                    valX_prev = Log3DX[m_binary_prev][i-1][j-1];
                    valY_prev = Log3DY[m_binary_prev][i-1][j-1];

                    if (std::isinf(valM_prev) && std::isinf(valX_prev) && std::isinf(valY_prev)) {
                    //if(!bool_MXY[m_binary_prev][i-1][j-1]){

                        Log3DM[m_binary_this][i][j] = -std::numeric_limits<double>::infinity();

                    }else{

                        valM_this = computeLK_MXY_local(NU,
                                                        valM_prev,
                                                        valX_prev,
                                                        valY_prev,
                                                        Log2DM[id1][id2],
                                                        m);

                        if (std::isnan(valM_this)) {
                            LOG(FATAL)
                                    << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is nan. Check call stack below.";
                        }

                        Log3DM[m_binary_this][i][j] = valM_this;

                    }

                }
                //=======================================
                // GAPX
                if( (i-1) >=0 ){

                    //id1 = rev_map_compressed_seqs_.at(nodeID_L).at(i-1);
                    id1 = map_compressed_seqs_.at(nodeID_L).at(i-1);

                    bool_MXY[m_binary_this][i][j]=true;

                    valM_prev = Log3DM[m_binary_prev][i-1][j];
                    valX_prev = Log3DX[m_binary_prev][i-1][j];
                    valY_prev = Log3DY[m_binary_prev][i-1][j];

                    if (std::isinf(valM_prev) && std::isinf(valX_prev) && std::isinf(valY_prev)) {
                    //if(!bool_MXY[m_binary_prev][i-1][j]){

                        Log3DX[m_binary_this][i][j] = -std::numeric_limits<double>::infinity();

                    }else{

                        valX_this = computeLK_MXY_local(NU,
                                                        valM_prev,
                                                        valX_prev,
                                                        valY_prev,
                                                        Log2DX[id1],
                                                        m);

                        if (std::isnan(valX_this)) {
                            LOG(FATAL)
                                    << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is nan. Check call stack below.";
                        }

                        Log3DX[m_binary_this][i][j] = valX_this;

                    }

                }
                //=======================================
                // GAPY
                if( (j-1)>=0 ){

                    //id2 = rev_map_compressed_seqs_.at(nodeID_R).at(j-1);
                    id2 = map_compressed_seqs_.at(nodeID_R).at(j-1);

                    bool_MXY[m_binary_this][i][j]=true;

                    valM_prev = Log3DM[m_binary_prev][i][j-1];
                    valX_prev = Log3DX[m_binary_prev][i][j-1];
                    valY_prev = Log3DY[m_binary_prev][i][j-1];

                    if (std::isinf(valM_prev) && std::isinf(valX_prev) && std::isinf(valY_prev)) {
                    //if(!bool_MXY[m_binary_prev][i][j-1]){

                        Log3DY[m_binary_this][i][j] = -std::numeric_limits<double>::infinity();

                    }else{

                        valY_this = computeLK_MXY_local(NU,
                                                        valM_prev,
                                                        valX_prev,
                                                        valY_prev,
                                                        Log2DY[id2],
                                                        m);

                        if (std::isnan(valY_this)) {
                            LOG(FATAL)
                                    << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is nan. Check call stack below.";
                        }

                        Log3DY[m_binary_this][i][j] = valY_this;

                    }

                }
                //=======================================
                valM_this = Log3DM[m_binary_this][i][j];
                valX_this = Log3DX[m_binary_this][i][j];
                valY_this = Log3DY[m_binary_this][i][j];

                // Find which matrix contains the best value of LK found until this point.
                index_of_max(valM_this, valX_this, valY_this, epsilon, generator, distribution,max_val_index,true,tr,max_of_3);

                //switch (max_val_index.index) {
                switch (tr) {
                    case MATCH_STATE:
                        LK3D[m][i][j]=valM_this;
                        break;
                    case GAP_X_STATE:
                        LK3D[m][i][j]=valX_this;
                        break;
                    case GAP_Y_STATE:
                        LK3D[m][i][j]=valY_this;
                        break;
                    default:
                        LK3D[m][i][j]=-std::numeric_limits<double>::infinity();
                        //LOG(FATAL) <<"\nSomething went wrong in reading the TR value. TR is neither MATCH, nor GAPX, nor GAPY. ";
                }

                //std::cout<<"TR["<<m<<"]["<<i<<"]["<<j<<"]\n";

                // Store the index for the traceback
                TR[m][i][j] = tr; //max_val_index.index;

                // If we reached the corner of the 3D cube, then:
                if ( (i == (h - 1)) && (j == (w - 1)) && (m>=h) && (m>=w) ) {
                //if ( (i == (h - 1)) && (j == (w - 1))) {
                //if(bool_MXY[m_binary_this][i][j]){
                    // the algorithm is filling the last column of 3D DP matrix where
                    // all the characters are in the MSA

                    //double max_of_3 = max_of_three(valM_this, valX_this, valY_this, epsilon, true);
                    //max_of_MXY = max_val_index.val;

                    if(tr==(int)STOP_STATE){
                        LOG(FATAL) <<"\nSomething went wrong in reading the TR value. TR is neither MATCH, nor GAPX, nor GAPY. ";
                    }

                    //if (max_val_index.val > score) {
                    if (max_of_3 > score) {
                        score = max_of_3;//max_val_index.val;
                        level_max_lk = m;
                    }

                    //=====================================================================
                    // early stop condition
//                    if (score < prev_best_lk) {
//                        prev_best_lk = score;
//                        counter_to_early_stop++;
//                        if (counter_to_early_stop > max_decrease_before_stop) {
//                            // if for max_decrease_before_stop consecutive times
//                            // the lk decrease then exit, the maximum lk has been reached
//                            flag_exit = true;
//                        }
//                    } else {
//                        counter_to_early_stop = 0;
//                    }
                    //=====================================================================

                }
            }

        }

    }

    // level (k position) in the DP matrix that contains the highest lk value
    score_.at(nodeID) = score;

    lk_down_[nodeID].resize(level_max_lk);
    fv_data_[nodeID].resize(level_max_lk);

    //==========================================================================================

    //TODO: compute fv_empty_data_[nodeID]=
    //fv_empty_data_[nodeID]=


    // start backtracing the 3 matrices (MATCH, GAPX, GAPY)
    TracebackPath_t traceback_path(level_max_lk, ' ');
    id1 = h - 1;
    id2 = w - 1;
    int idmL,idmR;
    for (int lev = level_max_lk; lev > 0; lev--) {
        int state = TR[lev][id1][id2];
        switch (state) {
            case MATCH_STATE:

                lk_down_.at(nodeID).at(lev - 1)=LK3D[lev][id1][id2];

                idmL = map_compressed_seqs_.at(nodeID_L).at(id1-1);
                idmR = map_compressed_seqs_.at(nodeID_R).at(id2-1);

                fv_data_.at(nodeID).at(lev - 1) = Fv_M[idmL][idmR];

                id1 = id1 - 1;
                id2 = id2 - 1;

                traceback_path[lev - 1] = MATCH_CHAR;

                break;
            case GAP_X_STATE:

                lk_down_.at(nodeID).at(lev - 1)=LK3D[lev][id1][id2];

                idmL = map_compressed_seqs_.at(nodeID_L).at(id1-1);

                fv_data_.at(nodeID).at(lev - 1) = Fv_X[idmL];

                id1 = id1 - 1;

                traceback_path[lev - 1] = GAP_X_CHAR;

                break;
            case GAP_Y_STATE:

                lk_down_.at(nodeID).at(lev - 1)=LK3D[lev][id1][id2];

                idmR = map_compressed_seqs_.at(nodeID_R).at(id2-1);

                fv_data_.at(nodeID).at(lev - 1) = Fv_Y[idmR];

                id2 = id2 - 1;

                traceback_path[lev - 1] = GAP_Y_CHAR;

                break;
            default:
                LOG(FATAL) << "\nSomething went wrong during the alignment reconstruction in function pPIP::DP3D_PIP. Check call stack below.";
        }
    }

    traceback_path_.at(nodeID) = traceback_path;

    // converts traceback path into an MSA
    build_MSA(node, traceback_path);

    // assigns the sequence names of the new alligned sequences to the current MSA
    setMSAsequenceNames(node);


    compressMSA(node);



    //==========================================================================================

    //==========================================================================================
    // memory freeing
    //delete[] &LogM;
    //delete[] &LogX;
    //delete[] &LogY;
    //delete[] &LK;
    //delete[] &TR;
    //==========================================================================================
}

//void compressMSA(bpp::Node *node){
//
//    auto MSA = new bpp::VectorSequenceContainer(alphabet_);
//
//    int nodeID = node->getId();
//
//    std::vector<std::string> vect_seqs = pPIPUtils::siteContainer_2_sequence_vector(MSA_.at(nodeID));
//
//    for(int i=0;i<vect_seqs.size();i++){
//        MSA->addSequence(*(new bpp::BasicSequence(seqNames_.at(nodeID).at(i),
//                                                         vect_seqs.at(i),
//                                                         alphabet_)), true);
//    }
//
//    auto siteContainerMSA = new bpp::VectorSiteContainer(*MSA);
//
//    auto siteContComprMSA = bpp::PatternTools::shrinkSiteSet(*siteContainerMSA);
//
//    auto map = bpp::PatternTools::getIndexes(*siteContainerMSA,*siteContComprMSA);
//
//    std::vector<int> rev_map = pPIPUtils::reverse_map(map);
//
//    setMapsComprSeqsleaf(sonRight);
//
//}
void pPIP::DP3D_PIP(bpp::Node *node, bool local,bool flag_map) {

    // TODO: place as argument
    // used to select random when 2 or 3 lks (M,X,Y) have "exactly" the same value
    //bool randomSeed = true;

    // number of discrete gamma categories
    size_t num_gamma_categories = rDist_->getNumberOfCategories();

    if (local) {
        // recompute local tau, total tree length of a tree root at the given node
        _setTau(treemap_.left.at(node->getId()));

        // recompute the normalizing factor nu for the local tree
        _setNu();

        // recompute lambdas with the new normalizing factor (local tree), flag true = tree rooted here
        _setAllIotas(node, true);

        // recompute betas with the new normalizing factor (local tree), flag true = tree rooted here
        _setAllBetas(node, true);
    }

    int up_corner_i;
    int up_corner_j;
    int bot_corner_i;
    int bot_corner_j;
    int lw;
    int h, w;

    int tr;

    // Get the IDs of the sons nodes given the current node
    int nodeID = node->getId();

    tshlib::VirtualNode *vnode_left = treemap_.left.at(nodeID)->getNodeLeft(); // bpp::Node to tshlib::VirtualNode
    tshlib::VirtualNode *vnode_right = treemap_.left.at(nodeID)->getNodeRight(); // bpp::Node to tshlib::VirtualNode
    int sequenceID_1 = treemap_.right.at(vnode_left);
    int sequenceID_2 = treemap_.right.at(vnode_right);

    // Compute dimensions of the 3D block at current internal node.
    h = MSA_.at(sequenceID_1).size() + 1; // dimension of the alignment on the left side
    w = MSA_.at(sequenceID_2).size() + 1; // dimension of the alignment on the riht side

    int d = (h - 1) + (w - 1) + 1; // third dimension of the DP matrix

    // lk of a single empty column (full of gaps) with rate variation (gamma distribution)
    std::vector<double> pc0;

    // MSA columns
    MSAcolumn_t sLs; // left column
    MSAcolumn_t sRs; // right column
    MSAcolumn_t col_gap_Ls; // left column full of gaps
    MSAcolumn_t col_gap_Rs; //right column full of gaps

    int numLeavesLeft = seqNames_.at(sequenceID_1).size(); // number of leaves in the left sub-tree
    int numLeavesRight = seqNames_.at(sequenceID_2).size(); // number of leaves in the right sub-tree

    col_gap_Ls = createGapCol(numLeavesLeft); // create column of gaps for the left sub-tree
    col_gap_Rs = createGapCol(numLeavesRight); // create column of gaps for the right sub-tree

//    signed long seed;
//    if (randomSeed) {
//        seed = std::chrono::system_clock::now().time_since_epoch().count(); // "random" seed
//    } else {
//        seed = 0; // fixed seed
//    }

    std::default_random_engine generator(seed_);                   // jatiapp seed
    std::uniform_real_distribution<double> distribution(0.0, 1.0); // Uniform distribution for the selection of lks with the same value

    auto epsilon = DBL_EPSILON;


    std::vector< std::vector<double> > lkM_pattern;
    std::vector< std::vector<double> > lkX_pattern;
    std::vector< std::vector<double> > lkY_pattern;


    //***************************************************************************************
    //***************************************************************************************
    if (local) {
        // compute the lk of a column full of gaps
        pc0 = computeLK_GapColumn_local(node, col_gap_Ls, col_gap_Rs,false);
    } else {
        /*
        pc0 = compute_pr_gap_all_edges_s(node,
                                         col_gap_Ls,
                                         col_gap_Rs,
                                         pi,
                                         originalAlphabetSize,
                                         alphabet);
        */
    }
    //***************************************************************************************
    //***************************************************************************************

    auto **LogM = new double *[2]; // DP sparse matrix for MATCH case (only 2 layer are needed)
    auto **LogX = new double *[2]; // DP sparse matrix for GAPX case (only 2 layer are needed)
    auto **LogY = new double *[2]; // DP sparse matrix for GAPY case (only 2 layer are needed)

    auto **TR = new int *[d]; // 3D traceback matrix

    // val num of cells occupied in a layer
    int numcells = int((w * (h + 1)) / 2);

    // allocate memory for the 2 layers
    LogM[0] = new double[numcells];
    LogX[0] = new double[numcells];
    LogY[0] = new double[numcells];
    LogM[1] = new double[numcells];
    LogX[1] = new double[numcells];
    LogY[1] = new double[numcells];

    //============================================================
    // marginal likelihood for all empty columns with rate variation (gamma distribution)
    // phi(m,pc0,r) depends on the MSA length m

    // marginal phi marginalized over gamma categories
    double log_phi_gamma;
    double prev_log_phi_gamma; // to store old value

    auto **PHI = new double *[d];
    double PC0 = 0.0;
    double NU = 0.0;

    for (int i = 0; i < d; i++) {
        PHI[i] = new double[num_gamma_categories];
    }

    for (int catg = 0; catg < num_gamma_categories; catg++) {
        // log( P_gamma(r) * phi(0,pc0(r),r) ): marginal lk for all empty columns of an alignment of size 0
        // PHI[0][catg] = log(rDist_->getProbability((size_t)catg)) + (nu_.at(catg) * (pc0.at(catg) - 1.0));
        PC0 += rDist_->getProbability((size_t) catg) * pc0.at(catg);
        NU += rDist_->getProbability((size_t) catg) * nu_.at(catg);
    }



    // computes the marginal phi marginalized over all the gamma categories
    log_phi_gamma = NU * (PC0 - 1);
    //log_phi_gamma = PHI[0][0];
    //for (int catg = 1; catg < num_gamma_categories; catg++) {
    //    log_phi_gamma=pPIPUtils::add_lns(log_phi_gamma,PHI[0][catg]);
    //}
    //============================================================

    LogM[0][0] = log_phi_gamma;
    LogX[0][0] = log_phi_gamma;
    LogY[0][0] = log_phi_gamma;

    TR[0] = new int[1]();
    TR[0][0] = STOP_STATE;

    double max_of_3 = -std::numeric_limits<double>::infinity();

    signed long level_max_lk = INT_MIN;
    double val;
    int m_binary_this;
    int m_binary_prev;

    double valM;
    double valX;
    double valY;

    signed long idx;
    max_val_str max_val_index;

    int coordSeq_1;
    int coordSeq_2;
    int coordTriangle_this_i;
    int coordTriangle_this_j;
    int coordTriangle_prev_i;
    int coordTriangle_prev_j;

    double score = -std::numeric_limits<double>::infinity();

    int depth;

    int last_d = d - 1;
    int size_tr, tr_up_i, tr_up_j, tr_down_i, tr_down_j;
    std::map<MSAcolumn_t, double> lkM;
    std::map<MSAcolumn_t, double> lkX;
    std::map<MSAcolumn_t, double> lkY;

    //============================================================
    // early stop condition flag
    bool flag_exit = false;
    int counter_to_early_stop;
    int max_decrease_before_stop = 10;
    double prev_lk = -std::numeric_limits<double>::infinity();

    // ============================================================
    // For each slice of the 3D cube, compute the values of each cell

    for (int m = 1; m < d; m++) {

        if (flag_exit) {
            break;
        }

        // alternate the two layers
        m_binary_this = m % 2;
        m_binary_prev = (m + 1) % 2;

        //***************************************************************************************
        //***************************************************************************************
        for (int catg = 0; catg < num_gamma_categories; catg++) {
            // computes the marginal phi(m,pc0(r),r) with gamma by multiplying the starting value
            // phi(0,pco(r),r) = log( P_gamma(r) * exp( nu(r) * (pc0(r)-1) ) ) with
            // 1/m * nu(r) at each new layer
            PHI[m][catg] = PHI[m - 1][catg] - log((long double) m) + log((long double) nu_.at(catg));
        }

        // store old value
        prev_log_phi_gamma = log_phi_gamma;

        // computes the marginal phi marginalized over all the gamma categories
        log_phi_gamma = PHI[m][0];
        for (int catg = 1; catg < num_gamma_categories; catg++) {
            log_phi_gamma = pPIPUtils::add_lns(log_phi_gamma, PHI[m][catg]);
        }
        //***************************************************************************************
        //***************************************************************************************

        //***************************************************************************************
        //***************************************************************************************
        // COMPUTES MATCH LK
        set_indeces_M(up_corner_i,
                      up_corner_j,
                      bot_corner_i,
                      bot_corner_j,
                      m, h, w);

        if (checkboundary(up_corner_i,
                          up_corner_j,
                          bot_corner_i,
                          bot_corner_j,
                          h, w)) {

            lw = 0;
            for (int i = up_corner_i; i <= bot_corner_i; i++) {

                coordTriangle_this_i = i;
                coordSeq_1 = coordTriangle_this_i - 1;
                coordTriangle_prev_i = coordTriangle_this_i - 1;

                // get left MSA column
                sLs = (MSA_.at(sequenceID_1).at(coordSeq_1));

                for (int j = 0; j <= lw; j++) {

                    coordTriangle_this_j = up_corner_j - j;
                    coordSeq_2 = coordTriangle_this_j - 1;
                    coordTriangle_prev_j = coordTriangle_this_j - 1;

                    // get right MSA column
                    sRs = (MSA_.at(sequenceID_2).at(coordSeq_2));

                    idx = get_indices_M(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valM = LogM[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valM = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_X(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valX = LogX[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valX = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_Y(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valY = LogY[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valY = -std::numeric_limits<double>::infinity();
                    }

                    if (std::isinf(valM) && std::isinf(valX) && std::isinf(valY)) {
                        LOG(FATAL) << "\nSomething went wrong during the comparison of valM, valX, valY in function pPIP::DP3D_PIP. Check call stack below.";
                    }

                    if (local) {
                        // compute MATCH lk
                        //valM -= prev_log_phi_gamma; // to avoid summing the marginal phi twice
                        //valX -= prev_log_phi_gamma;
                        //valY -= prev_log_phi_gamma;
                        //val = computeLK_M_local(log_phi_gamma,
                        val = computeLK_M_local(NU,
                                                valM,
                                                valX,
                                                valY,
                                                node,
                                                sLs,
                                                sRs,
                                                m,
                                                lkM,
                                                lkM_pattern,
                                                flag_map,
                                                false,
                                                false);
                    } else {
                        /*
                        val=computeLK_M_all_edges_s_opt(valM,
                                                        valX,
                                                        valY,
                                                        nu,
                                                        node,
                                                        sLs, sRs,
                                                        pi,
                                                        m,
                                                        lkM,
                                                        originalAlphabetSize, alphabet);
                        */
                    }

                    if (std::isinf(val)) {
                        LOG(FATAL) << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is infinite. Check call stack below.";
                    }

                    if (std::isnan(val)) {
                        LOG(FATAL) << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is nan. Check call stack below.";
                    }

                    idx = get_indices_M(coordTriangle_this_i,
                                        coordTriangle_this_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m, h, w);

                    LogM[m_binary_this][idx] = val;
                }
                lw++;
            }
        }
        //***************************************************************************************
        //***************************************************************************************
        // COMPUTES GAPX LK
        set_indeces_X(up_corner_i,
                      up_corner_j,
                      bot_corner_i,
                      bot_corner_j,
                      m, h, w);
        tr_down_i = bot_corner_i;
        tr_down_j = bot_corner_j;
        if (checkboundary(up_corner_i,
                          up_corner_j,
                          bot_corner_i,
                          bot_corner_j,
                          h, w)) {

            lw = 0;
            for (int i = up_corner_i; i <= bot_corner_i; i++) {

                coordTriangle_this_i = i;
                coordTriangle_prev_i = coordTriangle_this_i - 1;
                coordSeq_1 = coordTriangle_this_i - 1;

                // get left MSA column
                sLs = (MSA_.at(sequenceID_1).at(coordSeq_1));

                for (int j = 0; j <= lw; j++) {

                    coordTriangle_this_j = up_corner_j - j;
                    coordTriangle_prev_j = coordTriangle_this_j;

                    idx = get_indices_M(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valM = LogM[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valM = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_X(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valX = LogX[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valX = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_Y(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valY = LogY[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valY = -std::numeric_limits<double>::infinity();
                    }

                    if (std::isinf(valM) && std::isinf(valX) && std::isinf(valY)) {
                        LOG(FATAL) << "\nSomething went wrong during the comparison of valM, valX, valY in function pPIP::DP3D_PIP. Check call stack below.";
                    }

                    if (local) {
                        // compute GAPX lk
                        //valM -= prev_log_phi_gamma; // to avoid summing the marginal phi twice
                        //valX -= prev_log_phi_gamma;
                        //valY -= prev_log_phi_gamma;
                        //val= computeLK_X_local(log_phi_gamma,
                        val = computeLK_X_local(NU,
                                                valM,
                                                valX,
                                                valY,
                                                node,
                                                sLs,
                                                col_gap_Rs,
                                                m,
                                                lkX,
                                                lkX_pattern,
                                                flag_map,
                                                false,
                                                false);
                    } else {
                        /*
                        val=computeLK_X_all_edges_s_opt(valM,
                                                        valX,
                                                        valY,
                                                        nu,
                                                        node,
                                                        sLs, col_gap_Rs,
                                                        pi,
                                                        m,
                                                        lkX,
                                                        originalAlphabetSize, alphabet);
                        */
                    }

                    if (std::isinf(val)) {
                        LOG(FATAL) << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is infinite. Check call stack below.";

                    }

                    if (std::isnan(val)) {
                        LOG(FATAL) << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is nan. Check call stack below.";
                    }
                    idx = get_indices_X(coordTriangle_this_i,
                                        coordTriangle_this_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m, h, w);

                    LogX[m_binary_this][idx] = val;
                }
                lw++;
            }

        }
        //***************************************************************************************
        //***************************************************************************************
        // COMPUTES GAPY LK
        set_indeces_Y(up_corner_i,
                      up_corner_j,
                      bot_corner_i,
                      bot_corner_j,
                      m, h, w);
        tr_up_i = up_corner_i;
        tr_up_j = up_corner_j;
        if (checkboundary(up_corner_i,
                          up_corner_j,
                          bot_corner_i,
                          bot_corner_j,
                          h, w)) {

            lw = 0;
            for (int i = up_corner_i; i <= bot_corner_i; i++) {
                coordTriangle_this_i = i;
                coordTriangle_prev_i = coordTriangle_this_i;
                for (int j = 0; j <= lw; j++) {

                    coordTriangle_this_j = up_corner_j - j;
                    coordTriangle_prev_j = coordTriangle_this_j - 1;
                    coordSeq_2 = coordTriangle_this_j - 1;

                    // get right MSA column
                    sRs = (MSA_.at(sequenceID_2).at(coordSeq_2));

                    idx = get_indices_M(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valM = LogM[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valM = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_X(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valX = LogX[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valX = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_Y(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valY = LogY[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valY = -std::numeric_limits<double>::infinity();
                    }

                    if (std::isinf(valM) && std::isinf(valX) && std::isinf(valY)) {
                        LOG(FATAL) << "\nSomething went wrong during the comparison of valM, valX, valY in function pPIP::DP3D_PIP. Check call stack below.";
                    }

                    if (local) {
                        // compute GAPY lk
                        //valM -= prev_log_phi_gamma; // to avoid summing the marginal phi twice
                        //valX -= prev_log_phi_gamma;
                        //valY -= prev_log_phi_gamma;
                        //val= computeLK_Y_local(log_phi_gamma,
                        val = computeLK_Y_local(NU,
                                                valM,
                                                valX,
                                                valY,
                                                node,
                                                col_gap_Ls,
                                                sRs,
                                                m,
                                                lkY,
                                                lkY_pattern,
                                                flag_map,
                                                false,
                                                false);
                    } else {
                        /*
                        val=computeLK_Y_all_edges_s_opt(valM,
                                                        valX,
                                                        valY,
                                                        nu,
                                                        node,
                                                        col_gap_Ls, sRs,
                                                        pi,
                                                        m,
                                                        lkY,
                                                        originalAlphabetSize, alphabet);
                         */
                    }

                    if (std::isinf(val)) {
                        LOG(FATAL) << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is infinite. Check call stack below.";
                    }

                    if (std::isnan(val)) {
                        LOG(FATAL) << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is nan. Check call stack below.";

                    }

                    idx = get_indices_Y(coordTriangle_this_i,
                                        coordTriangle_this_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m, h, w);

                    LogY[m_binary_this][idx] = val;
                }
                lw++;
            }

        }

        size_tr = (int) ceil((tr_down_i - tr_up_i + 1) * (tr_up_j - tr_down_j + 1 + 1) / 2);

        /*TODO: optimize size TR*/
        TR[m] = new int[size_tr]();

        set_indeces_T(up_corner_i,
                      up_corner_j,
                      bot_corner_i,
                      bot_corner_j,
                      m, h, w);

        if (checkboundary(up_corner_i,
                          up_corner_j,
                          bot_corner_i,
                          bot_corner_j,
                          h, w)) {

            lw = 0;
            for (int i = up_corner_i; i <= bot_corner_i; i++) {
                coordTriangle_this_i = i;
                for (int j = 0; j <= lw; j++) {
                    coordTriangle_this_j = up_corner_j - j;

                    double mval;
                    double xval;
                    double yval;

                    idx = get_indices_M(coordTriangle_this_i,
                                        coordTriangle_this_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m, h, w);
                    if (idx >= 0) {
                        mval = LogM[m_binary_this][idx];
                    } else {
                        mval = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_X(coordTriangle_this_i,
                                        coordTriangle_this_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m, h, w);
                    if (idx >= 0) {
                        xval = LogX[m_binary_this][idx];
                    } else {
                        xval = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_Y(coordTriangle_this_i,
                                        coordTriangle_this_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m, h, w);
                    if (idx >= 0) {
                        yval = LogY[m_binary_this][idx];
                    } else {
                        yval = -std::numeric_limits<double>::infinity();
                    }


                    // TODO:: remove these 3 lines
                    mval = fabs((long double) mval) < epsilon ? -std::numeric_limits<double>::infinity() : mval;
                    xval = fabs((long double) xval) < epsilon ? -std::numeric_limits<double>::infinity() : xval;
                    yval = fabs((long double) yval) < epsilon ? -std::numeric_limits<double>::infinity() : yval;

                    //int ttrr;

                    index_of_max(mval, xval, yval, epsilon, generator, distribution,max_val_index,false,tr,max_of_3);

                    idx = get_indices_T(coordTriangle_this_i,
                                        coordTriangle_this_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m, h, w);

                    if (TR[m][idx] != 0) {
                        LOG(FATAL) << "\nSomething went wrong in accessing TR at indices:[" << m << "][" << idx << "] in function pPIP::DP3D_PIP. Check call stack below.";
                    }

                    TR[m][idx] = max_val_index.index;

                    if ((coordTriangle_this_i == (h - 1)) & (coordTriangle_this_j == (w - 1))) {
                        // the algorithm is filling the last column of 3D DP matrix where
                        // all the characters are in the MSA

                        //max_of_3 = max_of_three(mval, xval, yval, epsilon, false);


                        //max_of_3 = max_val_index.val;

                        if (max_val_index.val > score) {
                            score = max_val_index.val;
                            level_max_lk = m;
                        }

                        //=====================================================================
                        // early stop condition
                        if (score < prev_lk) {
                            prev_lk = score;
                            counter_to_early_stop++;
                            if (counter_to_early_stop > max_decrease_before_stop) {
                                // if for max_decrease_before_stop consecutive times
                                // the lk decrease then exit, the maximum lk has been reached
                                flag_exit = true;
                            }
                        } else {
                            counter_to_early_stop = 0;
                        }
                        //=====================================================================

                    }

                }
                lw++;
            }
        }
    }

    // level (k position) in the DP matrix that contains the highest lk value
    depth = level_max_lk;

    score_.at(nodeID) = score;

    //==========================================================================================
    // start backtracing the 3 matrices (MATCH, GAPX, GAPY)
    TracebackPath_t traceback_path(depth, ' ');
    int id1 = h - 1;
    int id2 = w - 1;
    for (int lev = depth; lev > 0; lev--) {
        set_indeces_T(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, lev, h, w);
        idx = get_indices_T(id1, id2, up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, lev, h, w);
        int state = TR[lev][idx];
        switch (TR[lev][idx]) {
            case MATCH_STATE:
                id1 = id1 - 1;
                id2 = id2 - 1;
                traceback_path[lev - 1] = MATCH_CHAR;
                break;
            case GAP_X_STATE:
                id1 = id1 - 1;
                traceback_path[lev - 1] = GAP_X_CHAR;
                break;
            case GAP_Y_STATE:
                id2 = id2 - 1;
                traceback_path[lev - 1] = GAP_Y_CHAR;
                break;
            default:
                LOG(FATAL) << "\nSomething went wrong during the alignment reconstruction in function pPIP::DP3D_PIP. Check call stack below.";
        }
    }

    traceback_path_.at(nodeID) = traceback_path;

    // converts traceback path into an MSA
    build_MSA(node, traceback_path);

    // assigns the sequence names of the new alligned sequences to the current MSA
    setMSAsequenceNames(node);
    //==========================================================================================

    //==========================================================================================
    // memory freeing
    delete[] LogM;
    delete[] LogX;
    delete[] LogY;
    delete[] TR;
    //==========================================================================================
}


//void pPIP::DP3D_PIP_SB(bpp::Node *node,UtreeBppUtils::treemap *tm,double gamma_rate, bool local,double temperature,int num_SB){
//
//    //TODO: place as argument
//    bool randomSeed = true;
//
//    //TODO: re-implement gamma distribution
//    double lambda_gamma = lambda_ * gamma_rate;
//    double mu_gamma = mu_ * gamma_rate;
//
//    if(local){
//        _setTau(node);
//        _setNu();
//        _setAllIotas(node,true);
//        _setAllBetas(node,true);
//    }else{
//    }
//
//    int lw;
//    int h,w;
//
//    auto sons = node->getSons();
//
//    int s1ID = sons.at(LEFT)->getId();
//    int s2ID = sons.at(RIGHT)->getId();
//
//    int nodeID = node->getId();
//
//    h = MSA_.at(s1ID).size()+1;
//    w = MSA_.at(s2ID).size()+1;
//
//    int d=(h-1)+(w-1)+1;
//
//    double pc0;
//
//    std::string sLs;
//    std::string sRs;
//    std::string col_gap_Ls;
//    std::string col_gap_Rs;
//
//    int numLeavesLeft = seqNames_.at(s1ID).size();
//    int numLeavesRight = seqNames_.at(s2ID).size();
//
//    col_gap_Ls=createGapCol(numLeavesLeft);
//    col_gap_Rs=createGapCol(numLeavesRight);
//
//    signed long seed;
//    if(randomSeed){
//        seed = std::chrono::system_clock::now().time_since_epoch().count();
//    }else{
//        seed = 0;
//    }
//
//    std::default_random_engine generator(seed);
//    std::uniform_real_distribution<double> distribution(0.0,1.0);
//
//    auto epsilon=DBL_EPSILON;
//
//    //***************************************************************************************
//    //***************************************************************************************
//    if(local){
//        pc0 = computeLK_GapColumn_local(node, col_gap_Ls, col_gap_Rs);
//    }
//
//    //else{
//    //    pc0 = compute_pr_gap_all_edges_s(node,
//    //                                     col_gap_Ls,
//    //                                     col_gap_Rs,
//    //                                     pi,
//    //                                     originalAlphabetSize,
//    //                                     alphabet);
//
//    //}
//
//
//    double ***LogM = new double**[d];
//    double ***LogX = new double**[d];
//    double ***LogY = new double**[d];
//    double ***Mp = new double**[d];
//    double ***Xp = new double**[d];
//    double ***Yp = new double**[d];
//    int ***TR = new int**[d];
//    for(int i =0; i<d; i++){
//        LogM[i] = new double*[h];
//        LogX[i] = new double*[h];
//        LogY[i] = new double*[h];
//        Mp[i] = new double*[h];
//        Xp[i] = new double*[h];
//        Yp[i] = new double*[h];
//        TR[i] = new int*[h];
//        for(int j =0; j<h; j++){
//            LogM[i][j] = new double[w];
//            LogX[i][j] = new double[w];
//            LogY[i][j] = new double[w];
//            Mp[i][j] = new double[w];
//            Xp[i][j] = new double[w];
//            Yp[i][j] = new double[w];
//            TR[i][j] = new int[w];
//            for(int k = 0; k<w;k++){
//                LogM[i][j][k] = -std::numeric_limits<double>::infinity();
//                LogX[i][j][k] = -std::numeric_limits<double>::infinity();
//                LogY[i][j][k] = -std::numeric_limits<double>::infinity();
//                Mp[i][j][k] = -std::numeric_limits<double>::infinity();
//                Xp[i][j][k] = -std::numeric_limits<double>::infinity();
//                Yp[i][j][k] = -std::numeric_limits<double>::infinity();
//                TR[i][j][k] = 0;
//            }
//        }
//    }
//
//    LogM[0][0][0]=nu_*(pc0-1.0);
//    LogX[0][0][0]=nu_*(pc0-1.0);
//    LogY[0][0][0]=nu_*(pc0-1.0);
//
//    //TR[0] = new int[1]();
//    TR[0][0][0]=STOP_STATE;
//
//    double max_of_3=-std::numeric_limits<double>::infinity();
//
//    signed long level_max_lk=INT_MIN;
//    double val;
//    int m_binary_this;
//    int m_binary_prev;
//
//    double valM;
//    double valX;
//    double valY;
//
//    signed long idx;
//
//    int coordSeq_1;
//    int coordSeq_2;
//    int coordTriangle_this_i;
//    int coordTriangle_this_j;
//    int coordTriangle_prev_i;
//    int coordTriangle_prev_j;
//
//    double score=-std::numeric_limits<double>::infinity();
//    int start_depth;
//    int depth;
//
//    bool flag_exit=false;
//    int last_d=d-1;
//    int size_tr,tr_up_i,tr_up_j,tr_down_i,tr_down_j;
//    std::map<std::string,double> lkM;
//    std::map<std::string,double> lkX;
//    std::map<std::string,double> lkY;
//
//    int m,i,j;
//
//    for(m=1;m<d;m++) {
//
//        if (flag_exit) {
//            break;
//        }
//
//        for (i = 0; i < h; i++) {
//
//            coordSeq_1 = i;
//
//            sLs = (MSA_.at(s1ID).at(coordSeq_1));
//
//            for (j = 0; j < w; j++) {
//
//                coordSeq_2 = j;
//
//                sRs = (MSA_.at(s2ID).at(coordSeq_2));
//
//                if (i - 1 > 0 && j - 1 > 0) {
//                    if (!(isinf(LogM[m - 1][i - 1][j - 1])) | !(isinf(LogX[m - 1][i - 1][j - 1])) |
//                        !(isinf(LogY[m - 1][i - 1][j - 1]))) {
//
//
//                        if (local) {
//                            val = computeLK_M_local(valM,
//                                                    valX,
//                                                    valY,
//                                                    node,
//                                                    sLs,
//                                                    sRs,
//                                                    m,
//                                                    lkM);
//                        } else {
//                            /*
//                            val=computeLK_M_all_edges_s_opt(valM,
//                                                            valX,
//                                                            valY,
//                                                            nu,
//                                                            node,
//                                                            sLs, sRs,
//                                                            pi,
//                                                            m,
//                                                            lkM,
//                                                            originalAlphabetSize, alphabet);
//                            */
//                        }
//                        double lk_c = val;
//                        //lk_c=log(iotaV0*betaV0*P(idx_i,idx_j));
//                        double lk = -log(m - 1) + log(nu_) + lk_c;
//                        double l1 = pPIPUtils::add_lns(LogM[m - 1][i - 1][j - 1], LogX[m - 1][i - 1][j - 1]);
//                        double l2 = pPIPUtils::add_lns(l1, LogY[m - 1][i - 1][j - 1]);
//                        LogM[m][i][j] = pPIPUtils::add_lns(lk, l2);
//                        Mp[m][i][j] = lk_c;
//                    }
//                }
//
//                if (i - 1 > 0) {
//                    if (!(isinf(LogM[m - 1][i - 1][j])) || !(isinf(LogX[m - 1][i - 1][j])) ||
//                        !(isinf(LogY[m - 1][i - 1][j]))) {
//
//                        if (local) {
//                            val = computeLK_X_local(valM,
//                                                    valX,
//                                                    valY,
//                                                    node,
//                                                    sLs, col_gap_Rs,
//                                                    m,
//                                                    lkX);
//                        } else {
//                            /*
//                            val=computeLK_X_all_edges_s_opt(valM,
//                                                            valX,
//                                                            valY,
//                                                            nu,
//                                                            node,
//                                                            sLs, col_gap_Rs,
//                                                            pi,
//                                                            m,
//                                                            lkX,
//                                                            originalAlphabetSize, alphabet);
//                            */
//                        }
//
//                        double lk_c = val;
//                        //lk_c = log(iotaV0 * betaV0 * P(idx_i, idx_j) + iotaV1 + betaV1 * FV(idx_i));
//                        double lk = -log(m - 1) + log(nu_) + lk_c;
//                        double l1 = pPIPUtils::add_lns(LogM[m - 1][i - 1][j], LogX[m - 1][i - 1][j]);
//                        double l2 = pPIPUtils::add_lns(l1, LogY[m - 1][i - 1][j]);
//                        LogX[m][i][j] = pPIPUtils::add_lns(lk, l2);
//                        Xp[m][i][j] = lk_c;
//                    }
//                }
//
//                if (j - 1 > 0) {
//                    if (!(isinf(LogM[m - 1][i][j - 1])) || !(isinf(LogX[m - 1][i][j - 1])) ||
//                        !(isinf(LogY[m - 1][i][j - 1]))) {
//
//                        if (local) {
//                            val = computeLK_Y_local(valM,
//                                                    valX,
//                                                    valY,
//                                                    node,
//                                                    col_gap_Ls, sRs,
//                                                    m,
//                                                    lkY);
//                        } else {
//                            /*
//                            val=computeLK_Y_all_edges_s_opt(valM,
//                                                            valX,
//                                                            valY,
//                                                            nu,
//                                                            node,
//                                                            col_gap_Ls, sRs,
//                                                            pi,
//                                                            m,
//                                                            lkY,
//                                                            originalAlphabetSize, alphabet);
//                             */
//                        }
//                        double lk_c = val;
//                        //lk_c = log(iotaV0 * betaV0 * P(idx_i, idx_j) + iotaV2 + betaV2 * FV(idx_j));
//                        double lk = -log(m - 1) + log(nu_) + lk_c;
//                        double l1 = pPIPUtils::add_lns(LogM[m - 1][i][j - 1], LogX[m - 1][i][j - 1]);
//                        double l2 = pPIPUtils::add_lns(l1, LogY[m - 1][i][j - 1]);
//                        LogY[m][i][j] = pPIPUtils::add_lns(lk, l2);
//                        Yp[m][i][j] = lk_c;
//                    }
//                }
//            }
//        }
//    }
//
//    double pm;
//    double pmn;
//    double log_pm;
//    double px;
//    double pxn;
//    double log_px;
//    double py;
//    double pyn;
//    double log_py;
//    double z;
//    double lk;
//    double p0;
//    double random_number;
//    double log_P;
//    char T;
//    double max_M,max_X,max_Y;
//    int idx_M,idx_X,idx_Y;
//    int idxMax;
//    for(int sb=0;sb<num_SB;sb++) {
//
//        pPIPUtils::max_val_in_column(LogM,d,h,w,max_M,idx_M);
//        pPIPUtils::max_val_in_column(LogX,d,h,w,max_X,idx_X);
//        pPIPUtils::max_val_in_column(LogY,d,h,w,max_Y,idx_Y);
//
//        score=max_of_three(max_M,max_X,max_Y,epsilon);
//
//        idxMax = index_of_max(max_M,max_X,max_Y,epsilon,generator,distribution);
//        switch(idxMax){
//            case MATCH_STATE:
//                T = MATCH_CHAR;
//                score = max_M;
//                m = idx_M;
//                break;
//            case GAP_X_STATE:
//                T = GAP_X_STATE;
//                score = max_X;
//                m = idx_X;
//                break;
//            case GAP_Y_STATE:
//                T = GAP_Y_CHAR;
//                score = max_Y;
//                m = idx_Y;
//                break;
//            default:
//                perror("state not recognized");
//        }
//
//        i = h;
//        j = w;
//
//
//        double log_Zm = LogM[m][i][j];
//        double log_Zx = LogX[m][i][j];
//        double log_Zy = LogY[m][i][j];
//
//        if(isinf(log_Zm) && isinf(log_Zx) && isinf(log_Zy)){
//            perror("ERROR 1: Zm, Zx and Zy are inf");
//        }
//
//        double log_Zmx = pPIPUtils::add_lns(log_Zm, log_Zx);
//        double log_Z = pPIPUtils::add_lns(log_Zmx, log_Zy);
//
//        if(isinf(log_Z)){
//            perror("ERROR 2 Z: is inf");
//        }
//
//        if(isinf(log_Zm)){
//            pm = 0;
//            pmn = 0;
//        }else{
//            log_pm = log_Zm - log_Z;
//            pm = exp(log_pm);
//            pmn = exp(-(1 - pm) / temperature);
//        }
//
//        if(isinf(log_Zx)){
//            px = 0;
//            pxn = 0;
//        }else{
//            log_px = log_Zx - log_Z;
//            px = exp(log_px);
//            pxn = exp(-(1 - px) / temperature);
//        }
//
//        if(isinf(log_Zy)){
//            py = 0;
//            pyn = 0;
//        }else{
//            log_py = log_Zy - log_Z;
//            py = exp(log_py);
//            pyn = exp(-(1 - py) / temperature);
//        }
//
//        z = pmn + pxn + pyn;
//        pm = pmn/z;
//        px = pxn/z;
//        py = pyn/z;
//
//        TracebackPath_t traceback;
//
//        m = 1;
//        lk = -log(m) + log(nu_) + p0;
//
//        while (i > 1 || j > 1 || m > 1) {
//
//            random_number  = distribution(generator);
//
//            if (random_number < pm) {
//                log_P = Mp[m][i][j];
//                i = i - 1;
//                j = j - 1;
//                m = m - 1;
//                T = MATCH_CHAR;
//            }else if(random_number < (pm + px)){
//                log_P = Xp[m][i][j];
//                i = i - 1;
//                m = m - 1;
//                T = GAP_X_CHAR;
//            }else{
//                log_P = Yp[m][i][j];
//                j = j - 1;
//                m = m - 1;
//                T = GAP_Y_CHAR;
//            }
//
//            if(isinf(log_P)){
//                perror("ERROR 3: P inf");
//            }
//
//            lk = lk + log_P;
//
//            traceback.append(&T);
//
//            log_Zm = LogM[m][i][j];
//            log_Zx = LogX[m][i][j];
//            log_Zy = LogY[m][i][j];
//
//            if (isinf(log_Zm) && isinf(log_Zx) && isinf(log_Zy)) {
//                perror("ERROR 1: Zm, Zx and Zy are inf");
//            }
//
//            log_Zmx = pPIPUtils::add_lns(log_Zm, log_Zx);
//            log_Z = pPIPUtils::add_lns(log_Zmx, log_Zy);
//
//            if(isinf(log_Z)){
//                perror("ERROR 2 Z: is inf");
//            }
//
//            if(isinf(log_Zm)) {
//                pm = 0;
//                pmn = 0;
//            }else {
//                log_pm = log_Zm - log_Z;
//                pm = exp(log_pm);
//                pmn = exp(-(1 - pm) / temperature);
//            }
//
//            if(isinf(log_Zx)){
//                px = 0;
//                pxn = 0;
//            }else {
//                log_px = log_Zx - log_Z;
//                px = exp(log_px);
//                pxn = exp(-(1 - px) / temperature);
//            }
//
//            if(isinf(log_Zy)) {
//                py = 0;
//                pyn = 0;
//            }else {
//                log_py = log_Zy - log_Z;
//                py = exp(log_py);
//                pyn = exp(-(1 - py) / temperature);
//            }
//
//            z = pmn + pxn + pyn;
//            pm = pmn/z;
//            px = pxn/z;
//            py = pyn/z;
//
//        }
//
//        reverse(traceback.begin(),traceback.end());
//
//        traceback_path_ensemble_.at(nodeID).push_back(traceback);
//
//        score_ensemble_.at(nodeID).push_back(score);
//
//    }
//
//    //TODO:free memory
//
//}

void pPIP::setSubstModel(bpp::SubstitutionModel *smodel) {
    substModel_ = smodel;
}

void pPIP::setTree(const Tree *tree) {
    tree_ = new TreeTemplate<Node>(*tree);
}

void pPIP::setFVleaf(bpp::Node *node) {

    int idx;

    int nodeID = node->getId();

    MSA_t MSA = MSA_.at(nodeID);

    size_t num_gamma_categories = rDist_->getNumberOfCategories();

    int lenComprSeqs = rev_map_compressed_seqs_.at(nodeID).size();

    fv_data_[nodeID].resize(lenComprSeqs);

    for (int i = 0; i < lenComprSeqs; i++) {
        fv_data_[nodeID][i].resize(num_gamma_categories);
    }

    for (int i = 0; i < lenComprSeqs; i++) {

        idx = rev_map_compressed_seqs_.at(nodeID).at(i);
        MSAcolumn_t s = MSA.at(idx);

        bpp::ColMatrix<double> fv;
        fv.resize(extendedAlphabetSize_, 1);
        bpp::MatrixTools::fill(fv, 0.0);

        idx = alphabet_->charToInt(&s[0]);
        idx = idx < 0 ? alphabetSize_ : idx;

        fv(idx, 0) = 1.0;
        for (int catg = 0; catg < num_gamma_categories; catg++) {
            fv_data_[nodeID].at(i).at(catg) = fv;
        }

    }

    bpp::ColMatrix<double> fv;
    fv.resize(extendedAlphabetSize_, 1);
    bpp::MatrixTools::fill(fv, 0.0);
    fv(alphabetSize_, 0) = 1.0;

    fv_empty_data_[nodeID].resize(num_gamma_categories);
    for (int catg = 0; catg < num_gamma_categories; catg++) {
        fv_empty_data_[nodeID].at(catg) = fv;
    }

}

void pPIP::set_lk_leaf(bpp::Node *node) {

    int nodeID = node->getId();

    size_t num_gamma_categories = rDist_->getNumberOfCategories();

    int len_seq_comp = rev_map_compressed_seqs_.at(nodeID).size();

    lk_down_.at(nodeID).resize(len_seq_comp);

    for(int i = 0; i < len_seq_comp; i++) {

        double p = 0.0;
        for (int catg = 0; catg < num_gamma_categories; catg++) {

            double fv0 = MatrixBppUtils::dotProd(fv_data_[nodeID].at(i).at(catg), pi_);

            p += rDist_->getProbability((size_t) catg) * \
                       iotasNode_[nodeID][catg] * \
                       betasNode_[nodeID][catg] * \
                       fv0;
        }

        lk_down_.at(nodeID).at(i) = p;

    }

}

void pPIP::set_lk_empty_leaf(bpp::Node *node) {

    int nodeID = node->getId();

    size_t num_gamma_categories = rDist_->getNumberOfCategories();

    lk_empty_down_.at(nodeID).resize(num_gamma_categories);

    double p;
    for (int catg = 0; catg < num_gamma_categories; catg++) {

        double fv0 = MatrixBppUtils::dotProd(fv_empty_data_[nodeID].at(catg), pi_);

        p = rDist_->getProbability((size_t) catg) * (iotasNode_[nodeID][catg] - \
                   iotasNode_[nodeID][catg] * betasNode_[nodeID][catg] + \
                   iotasNode_[nodeID][catg] * betasNode_[nodeID][catg] * fv0);

        lk_empty_down_.at(nodeID).at(catg) = p;
    }



}

void pPIP::compressMSA(bpp::Node *node) {

    int nodeID = node->getId();

    MSA_t MSA = MSA_.at(nodeID);

    auto sequences = new bpp::VectorSequenceContainer(alphabet_);

    std::vector<std::string> seqs = pPIPUtils::siteContainer_2_sequence_vector(MSA);

    for(int i = 0; i < seqs.size(); i++) {
        sequences->addSequence(*(new bpp::BasicSequence(seqNames_.at(nodeID).at(i),
                                                        seqs.at(i),
                                                        alphabet_)), true);
    }

    auto siteContainer = new bpp::VectorSiteContainer(*sequences);
    auto siteContCompr = bpp::PatternTools::shrinkSiteSet(*siteContainer);
    auto map_seqs = bpp::PatternTools::getIndexes(*siteContainer, *siteContCompr);

    map_compressed_seqs_.at(nodeID) = map_seqs;

    std::vector<int> rev_map_seqs = pPIPUtils::reverse_map(map_seqs);

    rev_map_compressed_seqs_.at(nodeID) = rev_map_seqs;

}
void pPIP::PIPAligner(std::vector<tshlib::VirtualNode *> &list_vnode_to_root,
                      bool local,
                      bool flag_RAM,
                      bool flag_map,
                      bool flag_pattern,
                      bool flag_fv) {

    // progressive PIP aligner
    // local: local subtree, rooted at the current node

    // resize vectors and maps
    _reserve(list_vnode_to_root);

    // set lambdas with rate variation (gamma distribution)
    _setLambda(substModel_->getParameter("lambda").getValue());

    // set mus with rate variation (gamma distribution)
    _setMu(substModel_->getParameter("mu").getValue());

    // copy pi
    _setPi(substModel_->getFrequencies());

    // set substitution/deletion probabilities with rate variation (gamma distribution)
    _getPrFromSubstutionModel(list_vnode_to_root);

    /*
    if(!local){
        utree_->addVirtualRootNode();
        _setTau(utree_->rootnode);
        utree_->removeVirtualRootNode();
        _setNu();
        _setAllIotas(list_vnode_to_root);
        _setAllBetas(list_vnode_to_root);
    }
     */

    size_t i = 1;

    for (auto &vnode:list_vnode_to_root) {
        // traverses the list of nodes and aligns the MSAs on the left and right side
        // if node is a leaf the resulting MSA is the sequence itself

        ApplicationTools::displayGauge(i, list_vnode_to_root.size());
        i++;

        auto node = tree_->getNode(treemap_.right.at(vnode), false);

        VLOG(1) << "[pPIP] Processing node " << node->getId();

        if (node->isLeaf()) {

            std::string seqname = sequences_->getSequencesNames().at((int) vnode->vnode_seqid);

            // associate the sequence name to the leaf node
            setMSAsequenceNames(node, seqname);

            // create a column containing the sequence associated to the leaf node
            setMSAleaves(node, sequences_->getSequence(seqname).toString());

        } else {

            // align using progressive 3D DP PIP

            if(flag_fv){
                DP3D_PIP_RAM_FAST(node);
            }else {
                if (flag_RAM) {
                    DP3D_PIP_RAM(node, local, flag_map, flag_pattern); // local: tree rooted at the given node
                } else {
                    DP3D_PIP(node, local, flag_map); // local: tree rooted at the given node
                }
            }
        }
    }


}

bpp::SiteContainer *pPIPUtils::pPIPmsa2Sites(bpp::pPIP *progressivePIP) {
    auto MSA = progressivePIP->getMSA(progressivePIP->getRootNode());

    auto sequences = new bpp::VectorSequenceContainer(progressivePIP->getAlphabet());

    auto seqNames = progressivePIP->getSeqnames(progressivePIP->getRootNode());

    int msaLen = MSA.size();

    int numLeaves = seqNames.size();
    for (int j = 0; j < numLeaves; j++) {
        std::string seqname = seqNames.at(j);
        std::string seqdata;
        seqdata.resize(msaLen);
        for (int i = 0; i < msaLen; i++) {
            seqdata.at(i) = MSA.at(i).at(j);
        }
        sequences->addSequence(*(new bpp::BasicSequence(seqname, seqdata, progressivePIP->getAlphabet())), true);
    }

    return new bpp::VectorSiteContainer(*sequences);
}

double pPIPUtils::add_lns(double a_ln, double b_ln) {
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

void pPIPUtils::max_val_in_column(double ***M, int depth, int height, int width, double &val, int &level) {

    val = -std::numeric_limits<double>::infinity();
    level = 0;

    for (int k = 0; k < depth; k++) {
        if (M[k][height - 1][width - 1] > val) {
            val = M[k][height - 1][width - 1];
            level = k;
        }
    }


}

std::vector<std::string> pPIPUtils::siteContainer_2_sequence_vector(std::vector<bpp::pPIP::MSAcolumn_t> &MSA){

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

std::vector<int> pPIPUtils::reverse_map(std::vector<int> &m){

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
