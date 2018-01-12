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
 * @file PIP.hpp
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
#ifndef MINIJATI_PIP_HPP
#define MINIJATI_PIP_HPP

#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Phyl/Model/AbstractSubstitutionModel.h>
#include <Bpp/Phyl/Model/Nucleotide/NucleotideSubstitutionModel.h>
#include <Bpp/Phyl/Model/Protein/ProteinSubstitutionModel.h>


using namespace bpp;
namespace bpp {
    class PIP_Nuc : public AbstractReversibleNucleotideSubstitutionModel {
private:
    double lambda_, mu_, tau_, nu_;
    std::string name_;

    double kappa_, r_;
    mutable double l_, k_, exp1_, exp2_;
    mutable RowMatrix<double> p_;

public:

    explicit PIP_Nuc(const NucleicAlphabet *alpha, double lambda = 0.1, double mu = 0.1, SubstitutionModel *basemodel = nullptr);

    virtual ~PIP_Nuc() = default;

    PIP_Nuc *clone() const { return new PIP_Nuc(*this); }

public:
    double Pij_t(size_t i, size_t j, double d) const;

    double dPij_dt(size_t i, size_t j, double d) const;

    double d2Pij_dt2(size_t i, size_t j, double d) const;

    //const Matrix<double> &getPij_t(double d) const;

    //const Matrix<double> &getdPij_dt(double d) const;

    //const Matrix<double> &getd2Pij_dt2(double d) const;

    void setFreqFromData(const SequenceContainer& data, double pseudoCount);

    size_t getNumberOfStates() const { return stateMap_.get()->getNumberOfModelStates(); }

    std::string getName() const { return name_; }

protected:

    void setFreq(std::map<int, double> &freqs);

    void updateMatrices();

};


}
#endif //MINIJATI_PIP_HPP

