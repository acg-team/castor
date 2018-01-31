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
 * @file RHomogeneousTreeLikelihood_PIP.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 11 01 2018
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

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/PatternTools.h>
#include <glog/logging.h>
#include "Utilities.hpp"

using namespace bpp;

#include "RHomogeneousTreeLikelihood_PIP.hpp"

RHomogeneousTreeLikelihood_PIP::RHomogeneousTreeLikelihood_PIP(
        const Tree &tree,
        TransitionModel *model,
        DiscreteDistribution *rDist,
        UtreeBppUtils::treemap *tm,
        bool checkRooted,
        bool verbose,
        bool usePatterns)
throw(Exception) :
        AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
        likelihoodData_(0),
        minusLogLik_(-1.),
        treemap_(*tm){


    init_(usePatterns);

}

RHomogeneousTreeLikelihood_PIP::RHomogeneousTreeLikelihood_PIP(
        const Tree &tree,
        const SiteContainer &data,
        TransitionModel *model,
        DiscreteDistribution *rDist,
        UtreeBppUtils::treemap *tm,
        bool checkRooted,
        bool verbose,
        bool usePatterns)
throw(Exception) :
        AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
        likelihoodData_(0),
        minusLogLik_(-1.),
        treemap_(*tm) {

    init_(usePatterns);
    setData(data);


}

RHomogeneousTreeLikelihood_PIP::RHomogeneousTreeLikelihood_PIP(
        const RHomogeneousTreeLikelihood_PIP &lik) :
        AbstractHomogeneousTreeLikelihood(lik),
        likelihoodData_(0),
        minusLogLik_(lik.minusLogLik_) {
    likelihoodData_ = dynamic_cast<DRASRTreeLikelihoodData *>(lik.likelihoodData_->clone());
    likelihoodData_->setTree(tree_);
}


RHomogeneousTreeLikelihood_PIP &RHomogeneousTreeLikelihood_PIP::operator=(const RHomogeneousTreeLikelihood_PIP &lik) {
    AbstractHomogeneousTreeLikelihood::operator=(lik);
    if (likelihoodData_) delete likelihoodData_;
    likelihoodData_ = dynamic_cast<DRASRTreeLikelihoodData *>(lik.likelihoodData_->clone());
    likelihoodData_->setTree(tree_);
    minusLogLik_ = lik.minusLogLik_;
    return *this;
}


RHomogeneousTreeLikelihood_PIP::~RHomogeneousTreeLikelihood_PIP() {
    delete likelihoodData_;
}


void RHomogeneousTreeLikelihood_PIP::init_(bool usePatterns) throw(Exception) {
    // This call initialises the data structure to compute the partial likelihoods of the nodes (it allows for ASRV distributions to be added on top of the s.m.)
    likelihoodData_ = new DRASRTreeLikelihoodData(tree_, rateDistribution_->getNumberOfCategories(), usePatterns);
    likelihoodEmptyData_ = new DRASRTreeLikelihoodData(tree_, rateDistribution_->getNumberOfCategories(), usePatterns);


}


void RHomogeneousTreeLikelihood_PIP::setData(const SiteContainer &sites) throw(Exception) {

    if (data_)
        delete data_;

    data_ = PatternTools::getSequenceSubset(sites, *tree_->getRootNode());

    tau_ = tree_->getTotalLength();
    computeNu();

    if (verbose_)
        ApplicationTools::displayTask("Initializing data structure");

    likelihoodData_->initLikelihoods(*data_, *model_);

    // Initialise gap column
    auto gaps = new std::vector<int>(sites.getNumberOfSequences(), sites.getAlphabet()->getGapCharacterCode());
    auto emptySite = new Site(*gaps, sites.getAlphabet());

    auto emptyContainer = new VectorSiteContainer(sites.getSequencesNames(), sites.getAlphabet());
    emptyContainer->addSite(*emptySite, false);

    // Prepare FV empty
    likelihoodEmptyData_->initLikelihoods(*emptyContainer, *model_);


    if (verbose_)
        ApplicationTools::displayTaskDone();

    nbSites_ = likelihoodData_->getNumberOfSites();
    nbDistinctSites_ = likelihoodData_->getNumberOfDistinctSites();
    nbStates_ = likelihoodData_->getNumberOfStates();


    // Initialise iotas and betas maps
    for(bpp::Node *node:tree_->getNodes()) {
        iotasData_.insert(std::make_pair(node, 0));
        betasData_.insert(std::make_pair(node, 0));
        indicatorFun_[node].resize(nbSites_);
    }

    // Initialise vector of likelihood nodes
    //likelihoodNodes_.clear();
    //computePostOrderNodeList(tree_->getRootNode());

    // Set all iotas
    setAllIotas();

    // Set all betas
    setAllBetas();

    // Initialise vectors for storing insertion histories values
    InitialiseInsertionHistories();


    if (verbose_)
        ApplicationTools::displayResult("Number of distinct sites", TextTools::toString(nbDistinctSites_));

    initialized_ = false;

}


double RHomogeneousTreeLikelihood_PIP::getLikelihood() const {
    double l = 1.;
    for (size_t i = 0; i < nbSites_; i++) {
        l *= getLikelihoodForASite(i);
    }
    return l;
}


double RHomogeneousTreeLikelihood_PIP::getLogLikelihood() const {

    return computeLikelihoodWholeAlignment();

}


double RHomogeneousTreeLikelihood_PIP::getLikelihoodForASite(size_t site) const {
    double l = 0;
    for (size_t i = 0; i < nbClasses_; i++) {
        l += getLikelihoodForASiteForARateClass(site, i) * rateDistribution_->getProbability(i);
    }
    return l;
}


double RHomogeneousTreeLikelihood_PIP::getLogLikelihoodForASite(size_t site) const {
    double l = 0;
    for (size_t i = 0; i < nbClasses_; i++) {
        double li = getLikelihoodForASiteForARateClass(site, i) * rateDistribution_->getProbability(i);
        if (li > 0) l += li; //Corrects for numerical instabilities leading to slightly negative likelihoods
    }
    return log(l);
}


double RHomogeneousTreeLikelihood_PIP::getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const {
    double l = 0;
    Vdouble *la = &likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass];
    for (size_t i = 0; i < nbStates_; i++) {
        //cout << (*la)[i] << "\t" << rootFreqs_[i] << endl;
        double li = (*la)[i] * rootFreqs_[i];
        if (li > 0) l += li; //Corrects for numerical instabilities leading to slightly negative likelihoods
    }
    return l;
}


double RHomogeneousTreeLikelihood_PIP::getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const {
    double l = 0;
    Vdouble *la = &likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass];
    for (size_t i = 0; i < nbStates_; i++) {
        l += (*la)[i] * rootFreqs_[i];
    }
    //if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
    return log(l);
}


double RHomogeneousTreeLikelihood_PIP::getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const {
    return likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass][static_cast<size_t>(state)];
}


double RHomogeneousTreeLikelihood_PIP::getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const {
    return log(likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass][static_cast<size_t>(state)]);
}


void RHomogeneousTreeLikelihood_PIP::setParameters(const ParameterList &parameters)
throw(ParameterNotFoundException, ConstraintException) {
    setParametersValues(parameters);
}


double RHomogeneousTreeLikelihood_PIP::getValue() const
throw(Exception) {
    if (!isInitialized()) throw Exception("RHomogeneousTreeLikelihood::getValue(). Instance is not initialized.");
    return minusLogLik_;
}


void RHomogeneousTreeLikelihood_PIP::computeTreeLikelihood() {

    computePostOrderNodeList(tree_->getRootNode());

    // Set insertion histories
    setInsertionHistories(*data_);

    computeSubtreeLikelihood();

}


void RHomogeneousTreeLikelihood_PIP::computeTreeLikelihood(std::vector<Node *> nodeList) {

    likelihoodNodes_ = nodeList;

    // Set insertion histories
    setInsertionHistories(*data_);

    computeSubtreeLikelihood();

}


void RHomogeneousTreeLikelihood_PIP::initializeLikelihoodMatrix_(VVVdouble *_likelihoods_node) {

    for (size_t i = 0; i < nbSites_; i++) { //For each site in the sequence
        for (size_t c = 0; c < nbClasses_; c++) { //For each rate classe
            for (size_t x = 0; x < nbStates_; x++) { //For each initial state
                (*_likelihoods_node)[i][c][x] = 1.;
            }
        }
    }

}


void RHomogeneousTreeLikelihood_PIP::initializeLikelihoodEmptyMatrix_(VVVdouble *_likelihoods_empty_node) {

    for (size_t c = 0; c < nbClasses_; c++) { //For each rate classe
        for (size_t x = 0; x < nbStates_; x++) { //For each initial state
            (*_likelihoods_empty_node)[0][c][x] = 1.;
        }
    }

}

/*
void RHomogeneousTreeLikelihood_PIP::hadamardMultFvSons_(Node *node) const{

    VVVdouble *_likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());

    size_t nbNodes = node->getNumberOfSons();

    double val;
    for (size_t i = 0; i < nbSites_; i++) {
        for (size_t c = 0; c < nbClasses_; c++) {
            for (size_t x = 0; x < nbStates_; x++) {
                val=1.0;
                for (size_t l = 0; l < nbNodes; l++) {

                    //VLOG(3) << "["<<i<<","<<c<<","<<x<<","<<l<<"]"<<std::endl;


                    const Node *son = node->getSon(l);
                    VVVdouble *_likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());
                    val *= (*_likelihoods_son)[i][c][x];
                    //val *= (*LikelihoodArraySonsPtr.at(l))[i][c][x];
                }
                (*_likelihoods_node)[i][c][x] = val;
            }
        }
    }

}*/

void RHomogeneousTreeLikelihood_PIP::hadamardMultFvSons_(Node *node) const{

    VVVdouble *_likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());

    //size_t nbNodes = node->getNumberOfSons();

    // Get mapped node on Utree
    std::vector<int> sonsIDs;
    tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
    tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

    sonsIDs.push_back(treemap_.right.at(vnode_left));
    sonsIDs.push_back(treemap_.right.at(vnode_right));
    size_t nbNodes = sonsIDs.size();

    double val;
    for (size_t i = 0; i < nbSites_; i++) {
        for (size_t c = 0; c < nbClasses_; c++) {
            for (size_t x = 0; x < nbStates_; x++) {
                val=1.0;
                for (size_t l = 0; l < nbNodes; l++) {

                    //VLOG(3) << "["<<i<<","<<c<<","<<x<<","<<l<<"]"<<std::endl;


                    //const Node *son = node->getSon(l);
                    //VVVdouble *_likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

                    VVVdouble *_likelihoods_son = &likelihoodData_->getLikelihoodArray(sonsIDs.at(l));
                    val *= (*_likelihoods_son)[i][c][x];
                    //val *= (*LikelihoodArraySonsPtr.at(l))[i][c][x];
                }
                (*_likelihoods_node)[i][c][x] = val;
            }
        }
    }

}

/*
void RHomogeneousTreeLikelihood_PIP::hadamardMultFvEmptySons_(Node *node) const {

    VVVdouble *_likelihoods_empty_node = &likelihoodEmptyData_->getLikelihoodArray(node->getId());

    size_t nbNodes = node->getNumberOfSons();

    double val;
    for (size_t c = 0; c < nbClasses_; c++) {
        for (size_t x = 0; x < nbStates_; x++) {
            val=1.0;
            for (size_t l = 0; l < nbNodes; l++) {
                const Node *son = node->getSon(l);
                VVVdouble *_likelihoods_empty_son = &likelihoodEmptyData_->getLikelihoodArray(son->getId());
                val *= (*_likelihoods_empty_son)[0][c][x];

                //val *= LikelihoodArrayEmptySonsPtr.at(l)->[0][c][x];
            }
            (*_likelihoods_empty_node)[0][c][x] = val;
        }
    }

}
*/
void RHomogeneousTreeLikelihood_PIP::hadamardMultFvEmptySons_(Node *node) const {

    VVVdouble *_likelihoods_empty_node = &likelihoodEmptyData_->getLikelihoodArray(node->getId());

    //size_t nbNodes = node->getNumberOfSons();

    // Get mapped node on Utree
    std::vector<int> sonsIDs;
    tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
    tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

    sonsIDs.push_back(treemap_.right.at(vnode_left));
    sonsIDs.push_back(treemap_.right.at(vnode_right));
    size_t nbNodes = sonsIDs.size();

    double val;
    for (size_t c = 0; c < nbClasses_; c++) {
        for (size_t x = 0; x < nbStates_; x++) {
            val=1.0;
            for (size_t l = 0; l < nbNodes; l++) {
                //const Node *son = node->getSon(l);
                VVVdouble *_likelihoods_empty_son = &likelihoodEmptyData_->getLikelihoodArray(sonsIDs.at(l));

                //VVVdouble *_likelihoods_empty_son = &likelihoodEmptyData_->getLikelihoodArray(son->getId());
                val *= (*_likelihoods_empty_son)[0][c][x];

                //val *= LikelihoodArrayEmptySonsPtr.at(l)->[0][c][x];
            }
            (*_likelihoods_empty_node)[0][c][x] = val;
        }
    }

}


/*
void RHomogeneousTreeLikelihood_PIP::SingleRateCategoryHadamardMultFvSons_(Node *node,unsigned long site,unsigned long rate,Vdouble *fv_out) const {

    size_t nbSons = node->getNumberOfSons();

    double val;
    for (size_t x = 0; x < nbStates_; x++) {
        val=1.0;
        for (size_t l = 0; l < nbSons; l++) {
            const Node *son = node->getSon(l);
            VVVdouble *_likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());
            val *= (*_likelihoods_son)[site][rate][x];
        }
        (*fv_out)[x] = val;
    }

}
*/

void RHomogeneousTreeLikelihood_PIP::SingleRateCategoryHadamardMultFvSons_(Node *node,unsigned long site,unsigned long rate,Vdouble *fv_out) const {

    // Get mapped node on Utree
    std::vector<int> sonsIDs;
    tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
    tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

    sonsIDs.push_back(treemap_.right.at(vnode_left));
    sonsIDs.push_back(treemap_.right.at(vnode_right));
    size_t nbSons = sonsIDs.size();

    double val;
    for (size_t x = 0; x < nbStates_; x++) {
        val=1.0;
        for (size_t l = 0; l < nbSons; l++) {
            VVVdouble *_likelihoods_son = &likelihoodData_->getLikelihoodArray(sonsIDs.at(l));
            val *= (*_likelihoods_son)[site][rate][x];
        }

        (*fv_out)[x] = val;
    }


}


/*
void RHomogeneousTreeLikelihood_PIP::SingleRateCategoryHadamardMultFvEmptySons_(Node *node,unsigned long rate,Vdouble *fv_out) const {

    size_t nbSons = node->getNumberOfSons();

    double val;
    for (size_t x = 0; x < nbStates_; x++) {
        val=1.0;
        for (size_t l = 0; l < nbSons; l++) {
            const Node *son = node->getSon(l);
            VVVdouble *_likelihoods_empty_son = &likelihoodEmptyData_->getLikelihoodArray(son->getId());
            val *= (*_likelihoods_empty_son)[0][rate][x];
        }
        (*fv_out)[x] = val;
    }

}
*/

void RHomogeneousTreeLikelihood_PIP::computePrTimesFv_(Node *node) const{

    VVVdouble *pxy__node = &pxy_[node->getId()];

    VVVdouble *_likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());

    double val;
    for (size_t c = 0; c < nbClasses_; c++) {
        VVdouble *pxy__node_c = &(*pxy__node)[c];
        for (size_t i = 0; i < nbSites_; i++) {
            std::vector<double> tmp = (*_likelihoods_node)[i][c];
            for (size_t x = 0; x < nbStates_; x++) {

                val=0.0;
                for (size_t y = 0; y < nbStates_; y++) {
                    val += (*pxy__node_c)[x][y] * tmp[y];
                }
                (*_likelihoods_node)[i][c][x]=val;
            }
        }
    }

    //alternative: first sites and then classes
    /*
    double val;
    for (size_t i = 0; i < nbSites_; i++) {
        for (size_t c = 0; c < nbClasses_; c++) {

            VVdouble *pxy__node_c = &(*pxy__node)[c];

            std::vector<double> tmp = (*_likelihoods_node)[i][c];

            for (size_t x = 0; x < nbStates_; x++) {

                val=0.0;
                for (size_t y = 0; y < nbStates_; y++) {
                    //val += (*pxy__node_c)[x][y] * (*_likelihoods_node)[i][c][y];
                    val += (*pxy__node_c)[x][y] * tmp[y];
                }
                (*_likelihoods_node)[i][c][x]=val;
            }
        }
    }
    */

}


void RHomogeneousTreeLikelihood_PIP::computePrTimesFvEmpty_(Node *node) const {

    VVVdouble *pxy__node = &pxy_[node->getId()];

    VVVdouble *_likelihoods_empty_node = &likelihoodEmptyData_->getLikelihoodArray(node->getId());

    double val;
    for (size_t c = 0; c < nbClasses_; c++) {
        VVdouble *pxy__node_c = &(*pxy__node)[c];
        std::vector<double> tmp = (*_likelihoods_empty_node)[0][c];
        for (size_t x = 0; x < nbStates_; x++) {

            val=0.0;
            for (size_t y = 0; y < nbStates_; y++) {
                     val += (*pxy__node_c)[x][y] * tmp[y];
            }
            (*_likelihoods_empty_node)[0][c][x]=val;
        }
    }

}


void RHomogeneousTreeLikelihood_PIP::computePrTimesIndicator_(Node *node) const {

    VVVdouble *pxy__node = &pxy_[node->getId()];

    VVVdouble *_likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());

    double val;
    for (size_t c = 0; c < nbClasses_; c++) {
        VVdouble *pxy__node_c = &(*pxy__node)[c];
        for (size_t i = 0; i < nbSites_; i++) {
            for (size_t x = 0; x < nbStates_; x++) {

                /*
                if ((*_likelihoods_node)[i][c][x] == 1) {
                    for (size_t y = 0; y < nbStates_; y++) {
                        (*_likelihoods_node)[i][c][y] = (*pxy__node_c)[y][x];
                    }
                    break;
                }
                */

                val=0.0;
                for (size_t y = 0; y < nbStates_; y++) {
                    val += (*pxy__node_c)[x][y] * indicatorFun_[node][i][y];

                }
                (*_likelihoods_node)[i][c][x]=val;
            }
        }
    }

    //alternative: first sites and then classes
    /*
    for (size_t i = 0; i < nbSites_; i++) {
        for (size_t c = 0; c < nbClasses_; c++) {

            VVdouble *pxy__node_c = &(*pxy__node)[c];

            for (size_t row = 0; row < nbStates_; row++) {
                if ((*_likelihoods_node)[i][c][row] == 1) {
                    for (size_t col = 0; col < nbStates_; col++) {
                        (*_likelihoods_node)[i][c][col] = (*pxy__node_c)[col][row];
                    }
                    break;
                }
            }
        }
    }
    */

}


void RHomogeneousTreeLikelihood_PIP::computePrTimesIndicatorEmpty_(Node *node) const {

    VVVdouble *pxy__node = &pxy_[node->getId()];

    VVVdouble *_likelihoods_empty_node = &likelihoodEmptyData_->getLikelihoodArray(node->getId());

    //double val;
    for (size_t c = 0; c < nbClasses_; c++) {
        VVdouble *pxy__node_c = &(*pxy__node)[c];

        for (size_t x = 0; x < nbStates_; x++) {
            (*_likelihoods_empty_node)[0][c][x] = (*pxy__node_c)[x][nbStates_ - 1];
        }

        /*
        for (size_t row = 0; row < nbStates_; row++) {
            if ((*_likelihoods_empty_node)[0][c][row] == 1) {
                for (size_t col = 0; col < nbStates_; col++) {
                    (*_likelihoods_empty_node)[0][c][col] = (*pxy__node_c)[col][nbStates_ - 1];
                }
                break;
            }
        }
        */
    }


}


void RHomogeneousTreeLikelihood_PIP::recombineFvAfterMove() const{

    for (auto &node:likelihoodNodes_) {
        //std::cout << "[RecombineFV]" << vnode->printNeighbours() << std::endl;
        recombineFvAtNode(node);
        //vnode->recombineFv();
        //vnode->_printFV();
    }
}


void RHomogeneousTreeLikelihood_PIP::recombineFvAtNode(Node *node) const{

    if (!node->isLeaf()) {
        // Internal node
        hadamardMultFvSons_(node);
        hadamardMultFvEmptySons_(node);

        if (node->hasFather()) {
            // Root --> special case for internal node
            computePrTimesFv_(node);
            computePrTimesFvEmpty_(node);
        }

    } else {
        // Leaf node
        computePrTimesIndicator_(node);
        computePrTimesIndicatorEmpty_(node);
    }

}


void RHomogeneousTreeLikelihood_PIP::computeSubtreeLikelihood() {

    for(auto &node:likelihoodNodes_) {

        if (!node->isLeaf()) {
            // Internal node
            hadamardMultFvSons_(node);
            hadamardMultFvEmptySons_(node);

            if (node->hasFather()) {
                // Root --> special case for internal node
                computePrTimesFv_(node);
                computePrTimesFvEmpty_(node);
            }

        } else {
            // Leaf node
            computePrTimesIndicator_(node);
            computePrTimesIndicatorEmpty_(node);
        }

        //--------------------------------------------------------------------------------------------
        // Debug ** array fetching + printing **
        VVVdouble *_likelihoods_empty_node = &likelihoodEmptyData_->getLikelihoodArray(node->getId());
        VVVdouble *_likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());
        printFV(node, _likelihoods_node);
        printFV(node, _likelihoods_empty_node);
        //--------------------------------------------------------------------------------------------

    }

}


void RHomogeneousTreeLikelihood_PIP::displayLikelihood(const Node *node) {
    VLOG(2) << "Likelihoods at node " << node->getName() << ": ";
    displayLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId()));
    VLOG(2) << "                                         ***";
}


void RHomogeneousTreeLikelihood_PIP::fireParameterChanged(const ParameterList &params) {
    applyParameters();

    if (rateDistribution_->getParameters().getCommonParametersWith(params).size() > 0
        || model_->getParameters().getCommonParametersWith(params).size() > 0) {
        //Rate parameter changed, need to recompute all probs:
        computeAllTransitionProbabilities();
    } else if (params.size() > 0) {
        //We may save some computations:
        for (size_t i = 0; i < params.size(); i++) {
            std::string s = params[i].getName();
            if (s.substr(0, 5) == "BrLen") {
                //Branch length parameter:
                computeTransitionProbabilitiesForNode(nodes_[TextTools::to<size_t>(s.substr(5))]);
            }
        }
        rootFreqs_ = model_->getFrequencies();
    }

    // Calls the routine to compute the FV values
    computeTreeLikelihood();

    minusLogLik_ = -getLogLikelihood();
}


void RHomogeneousTreeLikelihood_PIP::InitialiseInsertionHistories() const{

    for (int i = 0;i < nbSites_;i++) {
        for (auto &node:tree_->getNodes()) {

            // Initialize vectors descCount_ and setA_ and indicatorFunctionVector
            std::vector<int> descCount_;
            std::vector<bool> setA_;
            descCount_.resize(nbSites_);
            setA_.resize(nbSites_);
            indicatorFun_[node].at(i).resize(nbStates_);

            descCountData_.insert(std::make_pair(node->getId(), std::make_pair(descCount_, node)));
            setAData_.insert(std::make_pair(node->getId(), std::make_pair(setA_, node)));
        }
    }
}


void RHomogeneousTreeLikelihood_PIP::setInsertionHistories(const SiteContainer &sites) const {

    for(int i=0;i<nbSites_;i++){

        // Compute the total number of characters (exc. gap) for the current site
        int nonGaps_ = 0;

        for(int s=0;s<sites.getNumberOfSequences();s++ ){
            //siteValue = sites.getSite(i).getValue(s);
            if (sites.getAlphabet()->getGapCharacterCode() != sites.getSite(i).getValue(s)){
                nonGaps_++;
            }
        }

        for (auto &node:likelihoodNodes_) {

            // Computing descCount + setA
            if (node->isLeaf()){
                descCountData_[node->getId()].first.at(i) = (sites.getSequence(node->getName()).getValue(i) == sites.getAlphabet()->getGapCharacterCode() ? 0:1);

                indicatorFun_[node].at(i).at(sites.getSequence(node->getName()).getValue(i))=1;

            }else{

                std::vector<int> sonsIDs;
                tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
                tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

                sonsIDs.push_back(treemap_.right.at(vnode_left));
                sonsIDs.push_back(treemap_.right.at(vnode_right));
                size_t nbNodes = sonsIDs.size();

                //for(auto &son: node->getSons()){
                for (size_t l = 0; l < nbNodes; l++) {

                    descCountData_[node->getId()].first.at(i) += getNodeDescCountForASite(tree_->getNode(sonsIDs.at(l)),i); //descCountData_[son->getId()].first.at(i);
                }

            }

            setAData_[node->getId()].first.at(i) = (getNodeDescCountForASite(node,i) == nonGaps_);

        }
    }

}

/*
void RHomogeneousTreeLikelihood_PIP::setAllSetAData(const SiteContainer &sites) const {

    for(int i=0;i<nbSites_;i++){
        for(bpp::Node *node:tree_->getNodes()){

            // Compute the total number of characters (exc. gap) for the current site
            int nonGaps_ = 0;

            for(int s=0;s<sites.getNumberOfSequences();s++ ){
                //siteValue = sites.getSite(i).getValue(s);
                if (sites.getAlphabet()->getGapCharacterCode() != sites.getSite(i).getValue(s)){
                    nonGaps_++;
                }
            }

            setAData_[node->getId()].first.at(i) = (getNodeDescCountForASite(node,i) == nonGaps_);

        }
    }

}
*/

void RHomogeneousTreeLikelihood_PIP::setAllIotas() {

    Parameter mu_ = model_->getParameter("mu");
    double tau_ = tree_->getTotalLength();
    double T_;

    if (fabs(mu_.getValue()) < 1e-8) {
        perror("ERROR in set_iota: mu too small");
    }

    // Compute tau value;
    T_ = tau_ + 1 / mu_.getValue();

    if (fabs(T_) < 1e-8) {
        perror("ERROR in set_iota: T too small");
    }

    for (auto &node:tree_->getNodes()) {

        if (!node->hasFather()){

            iotasData_[node] = (1 / mu_.getValue()) / T_;

        } else {

            iotasData_[node] = node->getDistanceToFather() / T_;

        }

    }

}


void RHomogeneousTreeLikelihood_PIP::setAllBetas() {

    Parameter mu_ = model_->getParameter("mu");

    if (fabs(mu_.getValue()) < 1e-8) {
        perror("ERROR in set_iota: mu too small");
    }


    for (auto &node:tree_->getNodes()) {

        if (!node->hasFather()){

            betasData_[node] = 1.0;

        } else {
            betasData_[node] = (1.0 - exp(-mu_.getValue() * node->getDistanceToFather())) / (mu_.getValue() * node->getDistanceToFather());
            //node->vnode_beta = (1.0 - exp(-this->mu * vnode->vnode_branchlength)) / (this->mu * vnode->vnode_branchlength);
        }

    }

}


void RHomogeneousTreeLikelihood_PIP::resetNodeLikelihoodArrays(const Node *node) {

    //if (node->isLeaf()) return;
    size_t nbSites = likelihoodData_->getLikelihoodArray(node->getId()).size();

    // Must reset the likelihood array first (i.e. set all of them to 1):
    VVVdouble *_likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());
    VVVdouble *_likelihoods_node_empty = &likelihoodEmptyData_->getLikelihoodArray(node->getId());

    for (size_t i = 0; i < nbSites; i++) {
        //For each site in the sequence,
        VVdouble *_likelihoods_node_i = &(*_likelihoods_node)[i];
        VVdouble *_likelihoods_node_i_empty = &(*_likelihoods_node_empty)[0];

        for (size_t c = 0; c < nbClasses_; c++) {

            //For each rate classe,
            Vdouble *_likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
            Vdouble *_likelihoods_node_i_c_empty = &(*_likelihoods_node_i_empty)[c];

            for (size_t x = 0; x < nbStates_; x++) {
                //For each initial state,
                double initValue = 1.0;

                if(setAData_[node->getId()].first.at(i)){
                    if(x==0) VLOG(3) << "[FV   ] on site "<<i<<" -- setA active @ node: " << node->getName();
                    initValue = rootFreqs_[x] * iotasData_[node] * betasData_[node];
                    //if(node->isLeaf()) initValue *= rootFreqs_[x];
                }else{
                    if(node->isLeaf()) initValue = (*_likelihoods_node_i_c)[x] * rootFreqs_[x];
                }

                (*_likelihoods_node_i_c)[x] = initValue;

                // do it for empty column vector (only once)
                if(i==0){
                    if(setAData_[node->getId()].first.at(i)) {
                        if(x==0) VLOG(3) << "[EMPTY] on site "<<i<<" -- setA active @ node: " << node->getName();
                        initValue = rootFreqs_[x] * iotasData_[node] * (1-betasData_[node]) + betasData_[node];
                    }else{
                        if(node->isLeaf()) initValue = (*_likelihoods_node_i_c_empty)[x];
                    }
                    (*_likelihoods_node_i_c_empty)[x] = initValue;
                }

            }
        }
    }

}


double RHomogeneousTreeLikelihood_PIP::getLogLikelihoodSubtree(const Node *node) const {

    double ll = 0;
    std::vector<double> la(nbSites_);
    for (size_t i = 0; i < nbSites_; i++) {
        la[i] = getLogLikelihoodSubtreeForASite(i);
    }

    sort(la.begin(), la.end());
    for (size_t i = nbSites_; i > 0; i--) {
        ll += la[i - 1];
    }
    return ll;

}


double RHomogeneousTreeLikelihood_PIP::getLogLikelihoodSubtreeForASite(size_t site) const {
    double l = 0;
    for (size_t i = 0; i < nbClasses_; i++) {
        double li = getLogLikelihoodSubtreeForASiteForARateClass(site, i) * rateDistribution_->getProbability(i);
        if (li > 0) l += li; //Corrects for numerical instabilities leading to slightly negative likelihoods
    }
    return log(l);
}


double RHomogeneousTreeLikelihood_PIP::getLogLikelihoodSubtreeForASiteForARateClass(size_t site, size_t rateClass) const {

    double l = 0;
    Vdouble *la = &likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass];
    for (size_t i = 0; i < nbStates_; i++) {
        l += (*la)[i] * rootFreqs_[i];
    }
    //if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
    return log(l);
}


double RHomogeneousTreeLikelihood_PIP::computePhi(double lkEmptyColumn) const {

    double p;
    double log_factorial_m;
    size_t m = nbSites_;

    log_factorial_m = 0;
    for (int i = 1; i <= m; i++) {
        log_factorial_m += log(i);
    }

    p = -log_factorial_m + m * log(nu_) + (nu_ * (lkEmptyColumn - 1));

    return p;
}


void RHomogeneousTreeLikelihood_PIP::computeNu() {

    double lambda_ = model_->getParameter("lambda").getValue();
    double mu_ = model_->getParameter("mu").getValue();

    if (fabs(mu_) < 1e-8) {
        perror("ERROR in setNu: mu too small");
    }

    nu_ = lambda_ * (tau_ + 1 / mu_);

}


void RHomogeneousTreeLikelihood_PIP::computePostOrderNodeList(Node *startNode) const{

    if (!startNode->isLeaf()) {

        for(auto &son:startNode->getSons()){
            computePostOrderNodeList(son);

        }

        likelihoodNodes_.push_back(startNode);

    } else {

        likelihoodNodes_.push_back(startNode);

    }

}


void RHomogeneousTreeLikelihood_PIP::printFV(Node *node, VVVdouble *likelihoodvector){

    // Print FV if required ------------------------------------------------------------------------------------------------
    std::ostringstream sout;
    std::ostringstream sout_siteline;

    unsigned long sites = (*likelihoodvector).size();
    unsigned long classes = (*likelihoodvector)[0].size();
    unsigned long states = (*likelihoodvector)[0][0].size();

    VVVdouble lknode_i_c = (*likelihoodvector);

    for(int c=0;c<classes; c++) {

        VLOG(3) << "Node: " << node->getName() << " - ASVR class [" << c << "]";

        for (int x = 0; x < states; x++) {

            for (int i=0; i < sites; i++) {

                //unsigned long idxSite = likelihoodData_->getRootArrayPosition(i);

                sout_siteline << std::setfill('0') << std::setw(3) << i << "      | ";
                sout <<  std::fixed << std::setw( 8 ) << std::setprecision( 15 ) << (*likelihoodvector)[i][c][x] << " | ";
            }
            if(x==0){
                VLOG(3) << sout_siteline.str();
                sout_siteline.str("");
                sout_siteline.clear();
            }
            VLOG(3) << sout.str();
            sout.str("");
            sout.clear();
        }
    }
    VLOG(3) << "-----------------------------------------------------";


}


void RHomogeneousTreeLikelihood_PIP::printPrMatrix(Node *node, VVdouble *pr){
    std::ostringstream sout;
    VLOG(3) << "Node: " << node->getName() << " - distance [" << node->getDistanceToFather() << "]";

    for(size_t x=0; x<nbStates_; x++) {
        for (size_t y = 0; y < nbStates_; y++) {

            sout << std::fixed << std::setw( 8 ) << std::setprecision( 6 ) << (*pr)[x][y] << "\t";

        }
        VLOG(3) <<  sout.str() ;
        sout.str("");
        sout.clear();
        //VLOG(3) << "\n";
    }
    VLOG(3) << "---------------------------------";


}


double RHomogeneousTreeLikelihood_PIP::computeLikelihoodWholeAlignment()const {

    // compute loglk for the whole alignment
    double logLK = computeLikelihoodWholeSites();

    // compute likelihood empty column
    double lk_site_empty = computeLikelihoodWholeAlignmentEmptyColumn();

    // compute PHi
    double log_phi_value = computePhi(lk_site_empty);
    logLK += log_phi_value;

    return logLK;

}


double RHomogeneousTreeLikelihood_PIP::computeLikelihoodWholeSites() const {
    double logLK;// Initialise vector of likelihood per each site
    std::vector<double> lk_sites(nbSites_);

    // For each site in the alignment
    for (size_t i = 0; i < nbSites_; i++) {

        double lk_site = computeLikelihoodSite(i);

        lk_sites[i] = log(lk_site);
        VLOG(3) << "log_lk[" << i << "]=" << lk_sites[i] << std::endl;
    }

    // Sum all the values stored in the lk vector
    logLK= MatrixBppUtils::sumVector(&lk_sites);
    VLOG(3) << "LK Sites [BPP] "<< logLK;

    return logLK;
}


double RHomogeneousTreeLikelihood_PIP::computeLikelihoodWholeAlignmentEmptyColumn() const {
    // Add iota and beta quantities on nodes with actived SetA for empty column
    double lk_site_empty = 0;

    for (auto &node:likelihoodNodes_) {

        double fv_site = 0;

        for (size_t c = 0; c < nbClasses_; c++) {

            if (!node->isLeaf()) {

                auto fv_sons = new Vdouble(nbStates_);
                SingleRateCategoryHadamardMultFvEmptySons_(node, c, fv_sons);

                double partialLK = MatrixBppUtils::dotProd(fv_sons, &rootFreqs_);

                if (node->hasFather()) {
                    fv_site += iotasData_[node] * (1 - betasData_[node] + betasData_[node] * partialLK);
                } else {
                    fv_site += iotasData_[node] * betasData_[node] * partialLK;
                }


            } else {
                fv_site += iotasData_[node] * (1 - betasData_[node]);
            }


            VLOG(3) << "lk empty at node[" << node->getName() << "] -> " << fv_site << std::endl;


            double li = fv_site * rateDistribution_->getProbability(c);
            lk_site_empty += li;
        }

    }


    // Empty column
    VLOG(3) << "LK Empty [BPP] "<< lk_site_empty;
    return lk_site_empty;
}

double RHomogeneousTreeLikelihood_PIP::computeLikelihoodSite(size_t i) const {
    double lk_site = 0;

    for (auto &node:likelihoodNodes_) {

            double fv_site = 0;

            for (size_t c = 0; c < nbClasses_; c++) {
                //setAData_[treemap_.left.at(node->getId())].first.at(site)
                if (setAData_[node->getId()].first.at(i)) {

                    if (!node->isLeaf()) {

                        auto fv_sons = new Vdouble(nbStates_);
                        SingleRateCategoryHadamardMultFvSons_(node, i, c, fv_sons);

                        double partialLK = MatrixBppUtils::dotProd(fv_sons, &rootFreqs_);
                        fv_site += iotasData_[node] * betasData_[node] * partialLK;

                    } else {
                        std::vector<double> *indFun = &indicatorFun_[node][i];
                        double partialLK = MatrixBppUtils::dotProd(indFun, &rootFreqs_);
                        fv_site += (iotasData_[node] * betasData_[node]) * partialLK;
                    }

                }

                // Multiply the likelihood value of the site for the ASVR category probability
                double li = fv_site * rateDistribution_->getProbability(c);
                lk_site += li;
            }

        }
    return lk_site;
}


double RHomogeneousTreeLikelihood_PIP::computeLikelihoodOnTreeRearrangment(std::vector<tshlib::VirtualNode *> &listNodes) const {
    double logLK;

    // 0. convert the list of tshlib::VirtualNodes into bpp::Node
    likelihoodNodes_ = remapVirtualNodeLists(listNodes);

    // 1. Recombine FV arrays after move
    recombineFvAfterMove();

    // 2. Compute the lk of the empty column
    likelihoodNodes_.clear();
    computePostOrderNodeList(tree_->getRootNode());
    double lk_site_empty = computeLikelihoodWholeAlignmentEmptyColumn();

    // Set ancestral histories
    setInsertionHistories(*data_);

    // 3. Compute the likelihood of each site
    std::vector<double> lk_sites(nbSites_);
    std::vector<Node *> tempExtendedNodeList;

    for(unsigned long i = 0; i<nbSites_;i++) {

        // Extend it
        ExtendNodeListOnSetA(listNodes.back(), tempExtendedNodeList, i);

        // Overwrite the list of nodes on which computing the likelihood
        //likelihoodNodes_ = tempExtendedNodeList;

        // call to function which retrives the lk value for each site
        lk_sites[i] =  log(computeLikelihoodSite(i));

        tempExtendedNodeList.clear();

    }

    // Sum all the values stored in the lk vector
    logLK= MatrixBppUtils::sumVector(&lk_sites);
    VLOG(3) << "LK Sites [BPP] "<< logLK;


    // compute PHi
    double log_phi_value = computePhi(lk_site_empty);
    logLK += log_phi_value;


    return logLK;
}


void RHomogeneousTreeLikelihood_PIP::ExtendNodeListOnSetA(tshlib::VirtualNode *qnode, std::vector<Node *> &listNodes, unsigned long site) const {

    tshlib::VirtualNode *temp = qnode;
    // Get the corresponding node on the BPP tree
    Node *node = tree_->getNode(treemap_.right.at(temp));
    listNodes.push_back(node);

    do {

        if (node->isLeaf()) {
            break;
        }

        if (setAData_[treemap_.right.at(temp->getNodeLeft())].first.at(site)) {

            temp = temp->getNodeLeft();

        } else if (setAData_[treemap_.right.at(temp->getNodeRight())].first.at(site)) {

            temp = temp->getNodeRight();

        } else {

            break;
        }

        listNodes.push_back(tree_->getNode(treemap_.right.at(temp)));


    } while (setAData_[treemap_.right.at(temp)].first.at(site));

}


double RHomogeneousTreeLikelihood_PIP::getLogLikelihood(std::vector<tshlib::VirtualNode *> &listNodes) const {
    return computeLikelihoodOnTreeRearrangment(listNodes);
}


std::vector<Node *> RHomogeneousTreeLikelihood_PIP::remapVirtualNodeLists(std::vector<tshlib::VirtualNode *> &inputList) const {

    std::vector<Node *> newList;

    for(auto &vnode:inputList){

        newList.push_back(tree_->getNode(treemap_.right.at(vnode)));
    }

    return newList;
}


void RHomogeneousTreeLikelihood_PIP::SingleRateCategoryHadamardMultFvEmptySons_(Node *node, unsigned long rate, Vdouble *fv_out) const {

    //size_t nbSons = node->getNumberOfSons();
    std::vector<int> sonsIDs;
    tshlib::VirtualNode *vnode_left = treemap_.left.at(node->getId())->getNodeLeft();
    tshlib::VirtualNode *vnode_right = treemap_.left.at(node->getId())->getNodeRight();

    sonsIDs.push_back(treemap_.right.at(vnode_left));
    sonsIDs.push_back(treemap_.right.at(vnode_right));


    double val;
    for (size_t x = 0; x < nbStates_; x++) {
        val=1.0;
        for (size_t l = 0; l < sonsIDs.size(); l++) {
        //    const Node *son = node->getSon(l);
            VVVdouble *_likelihoods_empty_son = &likelihoodEmptyData_->getLikelihoodArray(sonsIDs.at(l));
            val *= (*_likelihoods_empty_son)[0][rate][x];
        }
        //VVVdouble *_likelihoods_empty_son_left = &likelihoodEmptyData_->getLikelihoodArray(son_leftID);
        //VVVdouble *_likelihoods_empty_son_right = &likelihoodEmptyData_->getLikelihoodArray(son_rightID);

        //val *= (*_likelihoods_empty_son_left)[0][rate][x];
        //val *= (*_likelihoods_empty_son_right)[0][rate][x];

        (*fv_out)[x] = val;
    }

}

