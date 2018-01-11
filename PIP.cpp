/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti & Massimo Maiolo
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
 * @file PIP.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 20 12 2017
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
#include "PIP.hpp"
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/PatternTools.h>
#include <glog/logging.h>

PIP_Nuc::PIP_Nuc(const NucleicAlphabet *alpha, double lambda, double mu, SubstitutionModel *basemodel) :
        AbstractParameterAliasable("PIP."),
        AbstractReversibleNucleotideSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "PIP."),
        lambda_(lambda), mu_(mu), r_(), l_(), k_(), exp1_(), exp2_(), p_(size_, size_) {


    // Inheriting basemodel parameters
    ParameterList parlist = basemodel->getParameters();

    for (int i = 0; i < parlist.size(); i++) {
        addParameter_(new Parameter("PIP." + parlist[i].getName(), parlist[i].getValue(), parlist[i].getConstraint()));
    }

    name_ = basemodel->getName() + "+PIP";

    generator_.resize(alpha->getSize() + 1, alpha->getSize() + 1);
    size_ = alpha->getSize() + 1;
    freq_.resize(alpha->getSize() + 1);


    // Copy the generator from substitution model + extend it
    const bpp::Matrix<double> &qmatrix = basemodel->getGenerator();

    int cols = qmatrix.getNumberOfColumns();
    int rows = qmatrix.getNumberOfRows();

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {

            if (i == j) {
                generator_(i, j) = qmatrix(i, j) - mu;
            } else {

                generator_(i, j) = qmatrix(i, j);
            }

        }

        generator_(i, cols) = mu;
    }

    // Add frequency for gap character
    freq_.at(alpha->getSize()) = 0;


    addParameter_(new Parameter("PIP.lambda", lambda, &Parameter::R_PLUS_STAR));
    addParameter_(new Parameter("PIP.mu", mu, &Parameter::R_PLUS_STAR));

    //updateMatrices();
}

void PIP_Nuc::updateMatrices() {


}

double PIP_Nuc::Pij_t(size_t i, size_t j, double d) const {

    return getPij_t(d)(i, j);

}

double PIP_Nuc::dPij_dt(size_t i, size_t j, double d) const {

    std::cerr << "PIP_Nuc::dPij_dt has been called";
    return 0;
}

double PIP_Nuc::d2Pij_dt2(size_t i, size_t j, double d) const {

    std::cerr << "PIP_Nuc::d2Pij_dt2 has been called";
    return 0;

}

const bpp::Matrix<double> &PIP_Nuc::getPij_t(double d) const {

    RowMatrix<double> md;

    MatrixTools::copy(generator_, md);
    MatrixTools::scale(md, d);
    MatrixTools::exp(md, md);

    return md;


}

const bpp::Matrix<double> &PIP_Nuc::getdPij_dt(double d) const {

    std::cerr << "PIP_Nuc::getdPij_dt has been called";


}

const bpp::Matrix<double> &PIP_Nuc::getd2Pij_dt2(double d) const {

    std::cerr << "PIP_Nuc::getd2dPij_dt2 has been called";


}

void PIP_Nuc::setFreq(std::map<int, double> &freqs) {}


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
    likelihoodData_ = new DRASRTreeLikelihoodData(tree_,
                                                  rateDistribution_->getNumberOfCategories(),
                                                  usePatterns);
}


void RHomogeneousTreeLikelihood_PIP::setData(const SiteContainer &sites) throw(Exception) {

    if (data_)
        delete data_;

    data_ = PatternTools::getSequenceSubset(sites, *tree_->getRootNode());

    if (verbose_)
        ApplicationTools::displayTask("Initializing data structure");

    likelihoodData_->initLikelihoods(*data_, *model_);

    if (verbose_)
        ApplicationTools::displayTaskDone();

    nbSites_ = likelihoodData_->getNumberOfSites();
    nbDistinctSites_ = likelihoodData_->getNumberOfDistinctSites();
    nbStates_ = likelihoodData_->getNumberOfStates();

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
    double ll = 0;
    std::vector<double> la(nbSites_);
    for (size_t i = 0; i < nbSites_; i++) {
        la[i] = getLogLikelihoodForASite(i);
    }
    std::sort(la.begin(), la.end());
    for (size_t i = nbSites_; i > 0; i--) {
        ll += la[i - 1];
    }
    return ll;
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


void RHomogeneousTreeLikelihood_PIP::computeTreeLikelihood()
{
    //computeSubtreeLikelihood(tree_->getRootNode());
}

void RHomogeneousTreeLikelihood_PIP::displayLikelihood(const Node* node)
{
    VLOG(2) << "Likelihoods at node " << node->getName() << ": ";
    displayLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId()));
    VLOG(2) << "                                         ***";
}


void RHomogeneousTreeLikelihood_PIP::fireParameterChanged(const ParameterList& params)
{
    applyParameters();

    if (rateDistribution_->getParameters().getCommonParametersWith(params).size() > 0
        || model_->getParameters().getCommonParametersWith(params).size() > 0)
    {
        //Rate parameter changed, need to recompute all probs:
        computeAllTransitionProbabilities();
    }
    else if (params.size() > 0)
    {
        //We may save some computations:
        for (size_t i = 0; i < params.size(); i++)
        {
            std::string s = params[i].getName();
            if (s.substr(0, 5) == "BrLen")
            {
                //Branch length parameter:
                computeTransitionProbabilitiesForNode(nodes_[TextTools::to<size_t>(s.substr(5))]);
            }
        }
        rootFreqs_ = model_->getFrequencies();
    }

    computeTreeLikelihood();

    minusLogLik_ = -getLogLikelihood();
}
