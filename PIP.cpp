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

/******************************************************************************/

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


RHomogeneousTreeLikelihood_PIP& RHomogeneousTreeLikelihood_PIP::operator=(const RHomogeneousTreeLikelihood_PIP& lik) {
    AbstractHomogeneousTreeLikelihood::operator=(lik);
    if (likelihoodData_) delete likelihoodData_;
    likelihoodData_ = dynamic_cast<DRASRTreeLikelihoodData*>(lik.likelihoodData_->clone());
    likelihoodData_->setTree(tree_);
    minusLogLik_ = lik.minusLogLik_;
    return *this;
}

RHomogeneousTreeLikelihood_PIP::~RHomogeneousTreeLikelihood_PIP()
{
    delete likelihoodData_;
}


void RHomogeneousTreeLikelihood_PIP::init_(bool usePatterns) throw (Exception)
{
    likelihoodData_ = new DRASRTreeLikelihoodData(
            tree_,
            rateDistribution_->getNumberOfCategories(),
            usePatterns);
}