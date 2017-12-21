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

namespace bpp {

    class PIP : public bpp::AbstractReversibleSubstitutionModel {
    private:

        double lambda;
        double mu;
        double tau;
        double nu;
        std::string name;

    public:

        PIP(bpp::SubstitutionModel *model, double valLambda, double valMu) :
                AbstractParameterAliasable(model->getName()),
                AbstractReversibleSubstitutionModel(model->getAlphabet(), &model->getStateMap(), model->getName()) {}

        virtual ~PIP() {}

        PIP* clone() const { return new  PIP(*this); }

        std::string getName() const { return name; }


    };
}
#endif //MINIJATI_PIP_HPP
