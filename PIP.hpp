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
#include <Bpp/Phyl/Model/FrequenciesSet/ProteinFrequenciesSet.h>

using namespace bpp;
namespace bpp {

    class PIP_Nuc : public AbstractReversibleNucleotideSubstitutionModel {
    private:
        double lambda_, mu_, tau_, nu_;
        std::string name_;
        std::string modelname_;

        double kappa_, r_;
        mutable double l_, k_, exp1_, exp2_;
        mutable RowMatrix<double> p_;
        mutable SubstitutionModel *submodel_;

    public:

        explicit PIP_Nuc(const NucleicAlphabet *alpha, double lambda = 0.1, double mu = 0.1, SubstitutionModel *basemodel = nullptr);

        virtual ~PIP_Nuc() = default;

        PIP_Nuc *clone() const { return new PIP_Nuc(*this); }

    public:
        //double Pij_t(size_t i, size_t j, double d) const;

        //double dPij_dt(size_t i, size_t j, double d) const;

        //double d2Pij_dt2(size_t i, size_t j, double d) const;

        //const Matrix<double> &getPij_t(double d) const;

        //const Matrix<double> &getdPij_dt(double d) const;

        //const Matrix<double> &getd2Pij_dt2(double d) const;

        void setFreqFromData(const SequenceContainer &data, double pseudoCount);

        size_t getNumberOfStates() const { return stateMap_.get()->getNumberOfModelStates(); }

        std::string getName() const { return name_; }

    protected:

        void setFreq(std::map<int, double> &freqs);

        void updateMatrices();

    };


    class PIP_AA : public AbstractReversibleProteinSubstitutionModel {
    private:
        double lambda_, mu_, tau_, nu_;
        std::string name_;
        std::string modelname_;
        ProteinFrequenciesSet *freqSet_;
        mutable SubstitutionModel *submodel_;

    public:
        /**
         * @brief Build a simple PIP_AA model, with original equilibrium frequencies.
         *
         * @param alpha A proteic alphabet.
         */
        PIP_AA(const ProteicAlphabet *alpha, double lambda = 0.1, double mu = 0.1, SubstitutionModel *basemodel = nullptr);

        /**
         * @brief Build a PIP_AA model with special equilibrium frequencies.
         *
         * @param alpha A proteic alphabet.
         * @param freqSet A pointer toward a protein frequencies set, which will be owned by this instance.
         * @param initFreqs Tell if the frequency set should be initialized with the original PIP_AA values.
         * Otherwise, the values of the set will be used.
         */
        //PIP_AA(const ProteicAlphabet *alpha, ProteinFrequenciesSet *freqSet, bool initFreqs = false,);

        PIP_AA(const PIP_AA &model) :
                AbstractParameterAliasable(model),
                AbstractReversibleProteinSubstitutionModel(model),
                freqSet_(dynamic_cast<ProteinFrequenciesSet *>(model.freqSet_->clone())) {}

        PIP_AA &operator=(const PIP_AA &model) {
            AbstractParameterAliasable::operator=(model);
            AbstractReversibleProteinSubstitutionModel::operator=(model);
            if (freqSet_) delete freqSet_;
            freqSet_ = dynamic_cast<ProteinFrequenciesSet *>(model.freqSet_->clone());
            return *this;
        }

        virtual ~PIP_AA() { delete freqSet_; }

        PIP_AA *clone() const { return new PIP_AA(*this); }

    public:
        std::string getName() const {
            if (freqSet_->getNamespace().find("PIP_AA+F.") != std::string::npos)
                return name_ + "+F";
            else
                return name_;
        }

        void fireParameterChanged(const ParameterList &parameters) {
            freqSet_->matchParametersValues(parameters);
            freq_ = freqSet_->getFrequencies();
            AbstractReversibleSubstitutionModel::fireParameterChanged(parameters);
        }

        void setNamespace(const std::string &prefix) {
            AbstractParameterAliasable::setNamespace(prefix);
            freqSet_->setNamespace(prefix + freqSet_->getName() + ".");
        }

        void setFrequenciesSet(const ProteinFrequenciesSet &freqSet) {
            delete freqSet_;
            freqSet_ = dynamic_cast<ProteinFrequenciesSet *>(freqSet.clone());
            resetParameters_();
            addParameters_(freqSet_->getParameters());
        }

        const FrequenciesSet *getFrequenciesSet() const { return freqSet_; }

        void setFreqFromData(const SequenceContainer &data, double pseudoCount = 0);

    protected:

        size_t getNumberOfStates() const { return 21; };

        void updateMatrices();

    };


}
#endif //MINIJATI_PIP_HPP

