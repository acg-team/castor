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
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Phyl/PatternTools.h>

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

    //¨generator_.resize(alpha->getSize(), alpha->getSize());
    //size_ = alpha->getSize();
    //freq_.resize(alpha->getSize());



    addParameter_(new Parameter("PIP.lambda", lambda, &Parameter::R_PLUS_STAR));
    addParameter_(new Parameter("PIP.mu", mu, &Parameter::R_PLUS_STAR));

    updateMatrices(basemodel);
}

void PIP_Nuc::updateMatrices(SubstitutionModel *basemodel) {

    lambda_ = getParameterValue("lambda");
    mu_ = getParameterValue("mu");


    //AbstractReversibleSubstitutionModel::updateMatrices();

    // Add frequency for gap character

    freq_ = basemodel->getFrequencies();
    freq_[alphabet_->getGapCharacterCode()] = 0; // hack for updateMatrices()


    // Copy the generator from substitution model + extend it
    const bpp::Matrix<double> &qmatrix = basemodel->getGenerator();

    int cols = qmatrix.getNumberOfColumns();
    int rows = qmatrix.getNumberOfRows();


    for (int i = 0; i < rows - 1; i++) {
        for (int j = 0; j < cols; j++) {

            if (i == j) {
                generator_(i, j) = qmatrix(i, j) - mu_;
            } else {

                generator_(i, j) = qmatrix(i, j);
            }

        }
        generator_(i, cols - 1) = mu_;
    }



    // Copy the exchangeability from substitution model + extend it
    const bpp::Matrix<double> &exMatrix = basemodel->getExchangeabilityMatrix();

    // Exchangeability
    for (int i = 0; i < rows - 1; i++) {
        for (int j = 0; j < cols - 1; j++) {

            if (i == j) {
                exchangeability_(i, j) = exMatrix(i, j) - mu_;
            } else {
                exchangeability_(i, j) = exMatrix(i, j);
            }

        }

        exchangeability_(i, cols - 1) = mu_;
    }

    // Normalization:
    setDiagonal();
    //normalize();






    // Compute eigen values and vectors:
    if (enableEigenDecomposition()) {
        EigenValue<double> ev(generator_);
        rightEigenVectors_ = ev.getV();
        eigenValues_ = ev.getRealEigenValues();
        iEigenValues_ = ev.getImagEigenValues();
        try {
            MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);
            isNonSingular_ = true;
            isDiagonalizable_ = true;
            for (size_t i = 0; i < size_ && isDiagonalizable_; i++) {
                if (abs(iEigenValues_[i]) > NumConstants::TINY())
                    isDiagonalizable_ = false;
            }
        }
        catch (ZeroDivisionException &e) {
            ApplicationTools::displayMessage("Singularity during diagonalization. Taylor series used instead.");

            isNonSingular_ = false;
            isDiagonalizable_ = false;
            MatrixTools::Taylor(generator_, 30, vPowGen_);
        }
    }

    //freq_[alphabet_->getGapCharacterCode()] = 0;

}

double PIP_Nuc::Pij_t(size_t i, size_t j, double d) const {

    return getPij_t(d)(i, j);

}

double PIP_Nuc::dPij_dt(size_t i, size_t j, double d) const {

    std::cout << "PIP_Nuc::dPij_dt has been called";
    return 0;
}

double PIP_Nuc::d2Pij_dt2(size_t i, size_t j, double d) const {

    std::cout << "PIP_Nuc::d2Pij_dt2 has been called";
    return 0;

}

/*
const bpp::Matrix<double> &PIP_Nuc::getPij_t(double d) const {

    MatrixTools::getId(size_, pijt_);

    MatrixTools::copy(generator_, pijt_);
    MatrixTools::scale(pijt_, d);
    MatrixTools::exp(pijt_, pijt_);

    return pijt_;


}
*/
//const bpp::Matrix<double> &PIP_Nuc::getdPij_dt(double d) const {

//std::cerr << "PIP_Nuc::getdPij_dt has been called";


//}

//const bpp::Matrix<double> &PIP_Nuc::getd2Pij_dt2(double d) const {

//std::cerr << "PIP_Nuc::getd2dPij_dt2 has been called";


//}

void PIP_Nuc::setFreq(std::map<int, double> &freqs) {
    std::vector<double> values;
    for (auto const &element : freqs) {
        values.push_back(element.second);
    }
    freq_ = values;

}

void PIP_Nuc::setFreqFromData(const SequenceContainer &data, double pseudoCount) {
    std::map<int, int> counts;
    SequenceContainerTools::getCounts(data, counts);
    std::map<int, double> freqs;

    int gapkey = data.getAlphabet()->getGapCharacterCode();
    std::map<int, int>::iterator iter = counts.find(gapkey);
    if (iter != counts.end())
        counts.erase(iter);
    else puts("not found");


    std::vector<int> retval;
    for (auto const &element : counts) {
        retval.push_back(element.first);
    }

    double t = 0;
    for (auto &ci : counts)
        t += ci.second;

    t += pseudoCount * (double) counts.size();

    //for (int i = 0; i < static_cast<int>(size_)-1; i++)
    for (auto &key:retval) {
        freqs[key] = (static_cast<double>(counts[key]) + pseudoCount) / t;
    }

    freqs[data.getAlphabet()->getGapCharacterCode()] = 0;
    // Re-compute generator and eigen values:
    setFreq(freqs);
}

