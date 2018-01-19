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

using namespace bpp;

#include "RHomogeneousTreeLikelihood_PIP.hpp"

RHomogeneousTreeLikelihood_PIP::RHomogeneousTreeLikelihood_PIP(
        const Tree &tree,
        TransitionModel *model,
        DiscreteDistribution *rDist,
        bool checkRooted,
        bool verbose,
        bool usePatterns)
throw(Exception) :
        AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
        likelihoodData_(0),
        minusLogLik_(-1.) {


    init_(usePatterns);

}

RHomogeneousTreeLikelihood_PIP::RHomogeneousTreeLikelihood_PIP(
        const Tree &tree,
        const SiteContainer &data,
        TransitionModel *model,
        DiscreteDistribution *rDist,
        bool checkRooted,
        bool verbose,
        bool usePatterns)
throw(Exception) :
        AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
        likelihoodData_(0),
        minusLogLik_(-1.) {

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

    // Add vectors for storing SetA array
    for(int i=0;i<nbSites_;i++){

        unsigned long idxSite = likelihoodData_->getRootArrayPositions().at(i);

        for(bpp::Node *node:tree_->getNodes()){

            // Initialize vectors descCount_ and setA_ and indicatorFunctionVector
            std::vector<int> descCount_;
            std::vector<bool> setA_;
            descCount_.resize(nbSites_);
            setA_.resize(nbSites_);
            indicatorFun_[node].at(i).resize(nbStates_);


            descCountData_.insert(std::make_pair(node->getId(),std::make_pair(descCount_,node)));
            setAData_.insert(std::make_pair(node->getId(),std::make_pair(setA_,node)));

            // Computing descCount + setA
            if (node->isLeaf()){

                descCountData_[node->getId()].first.at(i) = (sites.getSequence(node->getName()).getValue(i) == sites.getAlphabet()->getGapCharacterCode() ? 0:1);

                indicatorFun_[node].at(i).at(sites.getSequence(node->getName()).getValue(i))=1;

            }else{

                for(auto &son: node->getSons()){

                    descCountData_[node->getId()].first.at(i) += getNodeDescCountForASite(son,i); //descCountData_[son->getId()].first.at(i);

                }

            }

            // Compute the total number of characters (exc. gap) for the current site
            int nonGaps_ = 0;

            for(int s=0;s<sites.getNumberOfSequences();s++ ){
                //siteValue = sites.getSite(i).getValue(s);
                if (sites.getAlphabet()->getGapCharacterCode() != sites.getSite(i).getValue(s)) nonGaps_++;
            }

            setAData_[node->getId()].first.at(i) = (getNodeDescCountForASite(node,i) == nonGaps_);

        }
    }

    // Set all iotas
    setAllIotas();

    // Set all betas
    setAllBetas();

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

    //std::vector<Node *> postOrderList_Complete;

    computePostOrderNodeList(tree_->getRootNode());

    computeSubtreeLikelihood();

    //computeSubtreeLikelihood(tree_->getRootNode());

}

void RHomogeneousTreeLikelihood_PIP::computeTreeLikelihood(std::vector<Node *> nodeList) {


    likelihoodNodes_ = nodeList;

    //computePostOrderNodeList(tree_->getRootNode());

    computeSubtreeLikelihood();


}



void RHomogeneousTreeLikelihood_PIP::computeSubtreeLikelihood() {

    for(auto &node:likelihoodNodes_) {
        size_t nbSites = likelihoodData_->getLikelihoodArray(node->getId()).size();

        // Retrieve the likelihood arrays for the current node (fv + empty-column)
        VVVdouble *_likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());
        VVVdouble *_likelihoods_empty_node = &likelihoodEmptyData_->getLikelihoodArray(node->getId());

        if (!node->isLeaf()) {

            for (size_t i = 0; i < nbSites; i++) {
                //For each site in the sequence,
                VVdouble *_likelihoods_node_i = &(*_likelihoods_node)[i];

                for (size_t c = 0; c < nbClasses_; c++) {
                    //For each rate classe,
                    Vdouble *_likelihoods_node_i_c = &(*_likelihoods_node_i)[c];

                    for (size_t x = 0; x < nbStates_; x++) {
                        //For each initial state,
                        (*_likelihoods_node_i_c)[x] = 1.;

                        if (i == 0) {
                            (*_likelihoods_empty_node)[i][c][x] = 1.;
                        }
                    }
                }
            }


            size_t nbNodes = node->getNumberOfSons();


            for (size_t l = 0; l < nbNodes; l++) {
                //For each son node,

                const Node *son = node->getSon(l);
                std::vector<size_t> *_patternLinks_node_son = &likelihoodData_->getArrayPositions(node->getId(), son->getId());
                VVVdouble *_likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

                // For empty column
                VVVdouble *_likelihoods_empty_son = &likelihoodEmptyData_->getLikelihoodArray(son->getId());

                for (size_t i = 0; i < nbSites; i++) {
                    //For each site in the sequence,
                    size_t patternid = (*_patternLinks_node_son)[i];
                    VVdouble *_likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_node_son)[i]];
                    VVdouble *_likelihoods_node_i = &(*_likelihoods_node)[i];

                    // For empty column (only one site)
                    VVdouble *_likelihoods_empty_son_i = &(*_likelihoods_empty_son)[0];
                    VVdouble *_likelihoods_empty_node_i = &(*_likelihoods_empty_node)[0];
                    Vdouble *_likelihoods_empty_son_i_c;
                    Vdouble *_likelihoods_empty_node_i_c;

                    for (size_t c = 0; c < nbClasses_; c++) {
                        //For each rate classe,
                        Vdouble *_likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
                        Vdouble *_likelihoods_node_i_c = &(*_likelihoods_node_i)[c];

                        // For empty column (to be opened only at the first site of the alignment)
                        if (i == 0) {
                            //For each rate class,
                            _likelihoods_empty_son_i_c = &(*_likelihoods_empty_son_i)[c];
                            _likelihoods_empty_node_i_c = &(*_likelihoods_empty_node_i)[c];
                        }

                        // All-against-all states
                        for (size_t x = 0; x < nbStates_; x++) {

                            // store the likelihood value state-wise
                            (*_likelihoods_node_i_c)[x] *= (*_likelihoods_son_i_c)[x];

                            // For empty column (to be opened only at the first site of the alignment)
                            if (i == 0) (*_likelihoods_empty_node_i_c)[x] *= (*_likelihoods_empty_son_i_c)[x];
                        }
                    }
                }

            }

            if (node->hasFather()) {
                VVVdouble *pxy__node = &pxy_[node->getId()];
                Vdouble temp_fv_empty;
                for (size_t i = 0; i < nbSites; i++) {
                    //For each site in the sequence,
                    VVdouble *_likelihoods_node_i = &(*_likelihoods_node)[i];

                    for (size_t c = 0; c < nbClasses_; c++) {
                        //For each rate classe,
                        Vdouble *_likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
                        VVdouble *pxy__node_c = &(*pxy__node)[c];

                        Vdouble temp_fv = (*_likelihoods_node_i_c);

                        if (i == 0) {
                            temp_fv_empty = (*_likelihoods_empty_node)[i][c];
                        }

                        for (size_t x = 0; x < nbStates_; x++) {
                            //For each initial state,
                            //double temp_x = 0;
                            (*_likelihoods_node_i_c)[x] = 0;

                            if (i == 0) {
                                (*_likelihoods_empty_node)[i][c][x] = 0;
                            }

                            for (size_t y = 0; y < nbStates_; y++) {

                                (*_likelihoods_node_i_c)[x] += (*pxy__node_c)[x][y] * temp_fv[y];

                                if (i == 0) {
                                    (*_likelihoods_empty_node)[i][c][x] += (*pxy__node_c)[x][y] * temp_fv_empty[y];
                                }

                            }

                        }

                    }

                }

            }

        } else {


            VVVdouble *pxy__node = &pxy_[node->getId()];

            VVVdouble *_likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());

            // For empty column
            VVVdouble *_likelihoods_empty_node = &likelihoodEmptyData_->getLikelihoodArray(node->getId());


            for (size_t i = 0; i < nbSites; i++) {
                //For each site in the sequence,
                VVdouble *_likelihoods_node_i = &(*_likelihoods_node)[i];

                // For empty column (only one site)
                VVdouble *_likelihoods_empty_node_i = &(*_likelihoods_empty_node)[0];
                Vdouble *_likelihoods_empty_node_i_c;

                for (size_t c = 0; c < nbClasses_; c++) {
                    //For each rate classe,
                    Vdouble *_likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
                    VVdouble *pxy__node_c = &(*pxy__node)[c];

                    // For empty column (to be opened only at the first site of the alignment)
                    if (i == 0) {
                        //For each rate class,
                        _likelihoods_empty_node_i_c = &(*_likelihoods_empty_node_i)[c];
                    }

                    for (size_t col = 0; col < nbStates_; col++) {

                        if ((*_likelihoods_node_i_c)[col] == 1) {

                            for (size_t row = 0; row < nbStates_; row++) {
                                (*_likelihoods_node_i_c)[row] = (*pxy__node_c)[row][col];

                                // For empty column (to be opened only at the first site of the alignment)
                                if (i == 0) {
                                    (*_likelihoods_empty_node_i_c)[row] = (*pxy__node_c)[row][col];
                                }
                            }
                            break;
                        }
                    }
                }
            }

        }
        printFV(node, _likelihoods_node);
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

    computeTreeLikelihood();

    minusLogLik_ = -getLogLikelihood();
}


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


void RHomogeneousTreeLikelihood_PIP::recomputeFVarrays(const std::vector<Node *> nodelist) {

    size_t nbSites = likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId()).size();

    for (auto &node:nodelist) {

        size_t nbNodes = node->getNumberOfSons();
        VVVdouble *_likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());
        resetNodeLikelihoodArrays(node);

        for (size_t l = 0; l < nbNodes; l++) {
            //For each son node,
            const Node *son = node->getSon(l);

            VVVdouble *pxy__son = &pxy_[son->getId()];
            std::vector<size_t> *_patternLinks_node_son = &likelihoodData_->getArrayPositions(node->getId(), son->getId());
            VVVdouble *_likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

            for (size_t i = 0; i < nbSites; i++) {

                //For each site in the sequence,
                VVdouble *_likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_node_son)[i]];
                VVdouble *_likelihoods_node_i = &(*_likelihoods_node)[i];


                for (size_t c = 0; c < nbClasses_; c++) {
                    //For each rate classe,
                    Vdouble *_likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
                    Vdouble *_likelihoods_node_i_c = &(*_likelihoods_node_i)[c];

                    // cwise of each son to recompute the parent node (for site i and class c)
                    for (size_t e = 0; e < _likelihoods_node_i_c->size(); e++) {

                        (*_likelihoods_node_i_c)[e] *= (*_likelihoods_son_i_c)[e];

                    } // end for loop on states

                } // end for loop on classes

            } // end for loop on sites

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


double RHomogeneousTreeLikelihood_PIP::recomputeLikelihood(std::vector<Node *> nodeList) {

    recomputeFVarrays(nodeList);

    double ll = 0;
    double ll_empty = 0;
    // Initialise vector of likelihood per each site
    std::vector<double> la(nbSites_);
    std::vector<double> la_empty(nbSites_);

    // Per each site of the alignment compress the likelihood of each node and each asvr category
    for (size_t i = 0; i < nbSites_; i++) {

        double lic = 0;
        double lic_empty = 0;

        for (size_t c = 0; c < nbClasses_; c++) {

            double l = 0;
            double l_empty = 0;

            Vdouble *la = &likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(i)][c];
            Vdouble *la_empty = &likelihoodEmptyData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodEmptyData_->getRootArrayPosition(i)][c];

            for (size_t i = 0; i < nbStates_; i++) {

                l += (*la)[i] * rootFreqs_[i];
                l_empty += (*la_empty)[i] * rootFreqs_[i];

            }

            double li = log(l);
            double li_empty = log(l_empty);

            if (li > 0) lic += li;
            if (li_empty > 0) lic_empty += li_empty;
        }
        la[i] = log(lic);
        la_empty[i] = log(lic_empty);
    }

    // Sum the likelihood value for each site
    sort(la.begin(), la.end());
    for (size_t i = nbSites_; i > 0; i--) {
        ll += la[i - 1];
    }

    // Not necessary probably
    sort(la_empty.begin(), la_empty.end());
    for (size_t i = nbSites_; i > 0; i--) {
        ll_empty += la_empty[i - 1];
    }



    return ll;


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


void RHomogeneousTreeLikelihood_PIP::computePostOrderNodeList(Node *startNode) {

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



    for(int c=0;c<classes; c++) {

        VLOG(3) << "Node: " << node->getName() << " - ASVR class [" << c << "]";

        for (int x = 0; x < states; x++) {

            for (int i = 0; i < sites; i++) {

                VVdouble *_likelihoods_node_i = &(*likelihoodvector)[i];
                Vdouble *_likelihoods_node_i_c = &(*_likelihoods_node_i)[c];


                sout_siteline << std::setfill('0') << std::setw(3) << i << "      | ";


                sout <<  std::fixed << std::setw( 8 ) << std::setprecision( 6 ) << (*_likelihoods_node_i_c)[x] << " | ";
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

double RHomogeneousTreeLikelihood_PIP::computeLikelihoodWholeAlignment()const {

    double ll = 0;
    double ll_empty = 0;

    // Initialise vector of likelihood per each site
    std::vector<double> la(nbSites_);
    std::vector<double> la_empty(1);

    for (size_t i = 0; i < nbSites_; i++) {

        //unsigned long idxSite = likelihoodData_->getRootArrayPosition(i);

        double lk_site = 0;
        double lk_site_empty = 0;

        for (auto &node:likelihoodNodes_) {

            for (size_t c = 0; c < nbClasses_; c++) {

                double fv_site = 0;
                double fv_empty = 0;

                Vdouble *lknode_i_c = &likelihoodData_->getLikelihoodArray(node->getId())[likelihoodData_->getRootArrayPosition(i)][c];
                Vdouble *lknode_c_empty = &likelihoodEmptyData_->getLikelihoodArray(node->getId())[likelihoodEmptyData_->getRootArrayPosition(0)][c];


                if (setAData_[node->getId()].first.at(i)) {
                    VLOG(3) << "Likelihood for setA of " << node->getName() << "@"<<i;
                    if (!node->isLeaf()) {
                        // this compresses the vector into a scalar (dot product)
                        double partialLK = 0;

                        for (size_t s = 0; s < nbStates_; s++) {

                            partialLK += (*lknode_i_c)[s] * rootFreqs_[s];                        // this produces a scalar
                        }
                        fv_site += iotasData_[node] * betasData_[node] * partialLK;

                    } else {
                        double partialLK = 0;
                        //TODO: The rootFreq should be equal to the true site/node  @MAX  LG: DONE!
                        for (size_t s = 0; s < nbStates_; s++) {

                            partialLK += indicatorFun_[node][i][s] * rootFreqs_[s];                        // this produces a scalar
                        }
                        fv_site += (iotasData_[node] * betasData_[node]) * partialLK;
                    }


                    // If the setA is active for the node and the site, then multiply each entry of the fv vector for its frequency

                }


                // empty
                if(node->isLeaf()) {
                    double partialLK = 0;
                    for (size_t s = 0; s < nbStates_; s++) {

                        partialLK += indicatorFun_[node][i][s] * rootFreqs_[s];                        // this produces a scalar
                    }

                    fv_empty += iotasData_[node] * (1 - betasData_[node] + betasData_[node] * partialLK);


                }else{

                    double partialLK = 0;
                    for (size_t s = 0; s < nbStates_; s++) {
                        //TODO: The rootFreq should be equal to the true site/node  @MAX
                        partialLK += (*lknode_c_empty)[s] * rootFreqs_[s];                        // this produces a scalar
                    }

                    fv_empty += iotasData_[node] * (1 - betasData_[node] + betasData_[node] * partialLK);

                }


                // Multiply the likelihood value of the site for the ASVR category probability
                double li = fv_site * rateDistribution_->getProbability(c);
                double li_empty = fv_empty * rateDistribution_->getProbability(c);
                lk_site = li;
                lk_site_empty = li_empty;
            }

        }

        la[i] = log(lk_site);
        la_empty[0] = (lk_site_empty);

    }


    // Corrects for negative likelihood values
    //if (li > 0) lk_site += li;
    //if (li_empty > 0) lk_site_empty += li_empty;

    /*                  } else {
                          fv_site += (*lknode_i_c)[s];            //* rootFreqs_[i];          // this produces a scalar
                      }

                      fv_empty += (*lknode_c_empty)[s] * rootFreqs_[s];       //* rootFreqs_[i];          // this produces a scalar

                  }
              }else{

                  // If the setA is active for the node and the site, then add iota and beta too
                  if (setAData_[node->getId()].first.at(i)) {
                      fv_site *= iotasData_[node] * betasData_[node];                         // this produces a scalar
                  }

              }



              fv_empty *= betasData_[node];
              fv_empty += (iotasData_[node] * (1-betasData_[node]));



      la[i] = log(lk_site);
      la_empty[0] = (lk_site_empty);
  }

*/
  // Sum the likelihood value for each site
  sort(la.begin(), la.end());
  for (size_t i = nbSites_; i > 0; i--) {
      ll += la[i - 1];
  }


  VLOG(3) << "LK Sites [BPP] "<< ll;



  // Empty column
  ll_empty = la_empty[0];

  VLOG(3) << "LK Empty [BPP] "<< ll_empty;

  // compute PHi
  double log_phi_value = computePhi(ll_empty);
  ll += log_phi_value;

    return ll;

}

