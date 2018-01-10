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
#include <Bpp/Phyl/Likelihood/AbstractHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRASRTreeLikelihoodData.h>




using namespace bpp;

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

    const Matrix<double> &getPij_t(double d) const;

    const Matrix<double> &getdPij_dt(double d) const;

    const Matrix<double> &getd2Pij_dt2(double d) const;


    std::string getName() const { return name_; }

protected:

    void setFreq(std::map<int, double> &freqs);

    void updateMatrices();

};

namespace bpp {
    class RHomogeneousTreeLikelihood_PIP : public AbstractHomogeneousTreeLikelihood {
    private:

        mutable DRASRTreeLikelihoodData *likelihoodData_;

    protected:
        double minusLogLik_;

    public:
        /**
         * @brief Build a new RHomogeneousTreeLikelihood object without data.
         *
         * This constructor only initialize the parameters.
         * To compute a likelihood, you will need to call the setData() and the computeTreeLikelihood() methods.
         *
         * @param tree The tree to use.
         * @param model The substitution model to use.
         * @param rDist The rate across sites distribution to use.
         * @param checkRooted Tell if we have to check for the tree to be unrooted.
         * If true, any rooted tree will be unrooted before likelihood computation.
         * @param verbose Should I display some info?
         * @param usePatterns Tell if recursive site compression should be performed.
         * @throw Exception in an error occured.
         */
        RHomogeneousTreeLikelihood_PIP(
                const Tree &tree,
                TransitionModel *model,
                DiscreteDistribution *rDist,
                bool checkRooted = true,
                bool verbose = true,
                bool usePatterns = true)
        throw(Exception);

        /**
         * @brief Build a new RHomogeneousTreeLikelihood object with data.
         *
         * This constructor initializes all parameters, data, and likelihood arrays.
         *
         * @param tree The tree to use.
         * @param data Sequences to use.
         * @param model The substitution model to use.
         * @param rDist The rate across sites distribution to use.
         * @param checkRooted Tell if we have to check for the tree to be unrooted.
         * If true, any rooted tree will be unrooted before likelihood computation.
         * @param verbose Should I display some info?
         * @param usePatterns Tell if recursive site compression should be performed.
         * @throw Exception in an error occured.
         */
        RHomogeneousTreeLikelihood_PIP(
                const Tree &tree,
                const SiteContainer &data,
                TransitionModel *model,
                DiscreteDistribution *rDist,
                bool checkRooted = true,
                bool verbose = true,
                bool usePatterns = true)
        throw(Exception);

        RHomogeneousTreeLikelihood_PIP(const RHomogeneousTreeLikelihood_PIP &lik);

        RHomogeneousTreeLikelihood_PIP &operator=(const RHomogeneousTreeLikelihood_PIP &lik);

        virtual ~RHomogeneousTreeLikelihood_PIP();

        //RHomogeneousTreeLikelihood_PIP *clone() const { return new RHomogeneousTreeLikelihood_PIP(*this); }


    private:

        /**
         * @brief Method called by constructors.
         */
        void init_(bool usePatterns) throw(Exception);

    public:

        /**
         * @name The TreeLikelihood interface.
         *
         * Other methods are implemented in the AbstractHomogeneousTreeLikelihood class.
         *
         * @{
         */
        //void setData(const SiteContainer &sites) throw(Exception);

        //double getLikelihood() const;

        //double getLogLikelihood() const;

        //double getLikelihoodForASite(size_t site) const;

        //double getLogLikelihoodForASite(size_t site) const;
        /** @} */


        /**
         * @name The DiscreteRatesAcrossSites interface implementation:
         *
         * @{
         */
        //double getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

        //double getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

        //double getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const;

        //double getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const;
        /** @} */

        /**
         * @brief Implements the Function interface.
         *
         * Update the parameter list and call the applyParameters() method.
         * Then compute the likelihoods at each node (computeLikelihood() method)
         * and call the getLogLikelihood() method.
         *
         * If a subset of the whole parameter list is passed to the function,
         * only these parameters are updated and the other remain constant (i.e.
         * equal to their last value).
         *
         * @param parameters The parameter list to pass to the function.
         */
        //void setParameters(const ParameterList &parameters) throw(ParameterNotFoundException, ConstraintException);

        //double getValue() const throw(Exception);

        //size_t getSiteIndex(size_t site) const throw(IndexOutOfBoundsException) { return likelihoodData_->getRootArrayPosition(site); }

        /**
         * @name DerivableFirstOrder interface.
         *
         * @{
         */
        //double getFirstOrderDerivative(const std::string &variable) const throw(Exception);
        /** @} */

        /**
         * @name DerivableSecondOrder interface.
         *
         * @{
         */
        //double getSecondOrderDerivative(const std::string &variable) const throw(Exception);

        //double getSecondOrderDerivative(const std::string &variable1, const std::string &variable2) const throw(Exception) { return 0; } // Not implemented for now.
        /** @} */

    public:    // Specific methods:

        //DRASRTreeLikelihoodData *getLikelihoodData() { return likelihoodData_; }

        //const DRASRTreeLikelihoodData *getLikelihoodData() const { return likelihoodData_; }

        //void computeTreeLikelihood();

        //virtual double getDLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

        //virtual double getDLikelihoodForASite(size_t site) const;

        //virtual double getDLogLikelihoodForASite(size_t site) const;

        //virtual double getDLogLikelihood() const;

        //virtual void computeTreeDLikelihood(const std::string &variable);

        //virtual double getD2LikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

        //virtual double getD2LikelihoodForASite(size_t site) const;

        //virtual double getD2LogLikelihoodForASite(size_t site) const;

        //virtual double getD2LogLikelihood() const;

        //virtual void computeTreeD2Likelihood(const std::string &variable);


    protected:

        /**
         * @brief Compute the likelihood for a subtree defined by the Tree::Node <i>node</i>.
         *
         * @param node The root of the subtree.
         */
        //virtual void computeSubtreeLikelihood(const Node *node); //Recursive method.
        //virtual void computeDownSubtreeDLikelihood(const Node *);

        //virtual void computeDownSubtreeD2Likelihood(const Node *);

        //void fireParameterChanged(const ParameterList &params);

        /**
         * @brief This method is mainly for debugging purpose.
         *
         * @param node The node at which likelihood values must be displayed.
         */
        //virtual void displayLikelihood(const Node *node);

        //friend class RHomogeneousMixedTreeLikelihood;
    };
}
#endif //MINIJATI_PIP_HPP

