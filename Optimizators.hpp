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
 * @file Optimizators.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 24 02 2018
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
#ifndef MINIJATI_OPTIMIZATORS_HPP
#define MINIJATI_OPTIMIZATORS_HPP

namespace bpp {

    class Optimizators {

    public:

        Optimizators();

        virtual ~Optimizators();

        /**
       * @brief Optimize parameters according to options.
       *
       * @param tl               The TreeLikelihood function to optimize.
       * @param parameters       The initial list of parameters to optimize.
       *                         Use tl->getIndependentParameters() in order to estimate all parameters.
       * @param params           The attribute map where options may be found.
       * @param suffix           A suffix to be applied to each attribute name.
       * @param suffixIsOptional Tell if the suffix is absolutely required.
       * @param verbose          Print some info to the 'message' output stream.
       * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
       * @throw Exception        Any exception that may happen during the optimization process.
       * @return A pointer toward the final likelihood object.
       * This pointer may be the same as passed in argument (tl), but in some cases the algorithm
       * clone this object. We may change this bahavior in the future...
       * You hence should write something like
       * @code
       * tl = PhylogeneticsApplicationTools::optimizeParameters(tl, ...);
       * @endcode
       */
        static TreeLikelihood *optimizeParameters(
                TreeLikelihood *inTL,
                const ParameterList &parameters,
                std::map<std::string, std::string> &params,
                const std::string &suffix = "",
                bool suffixIsOptional = true,
                bool verbose = true,
                int warn = 1) throw(Exception);

    };

}
#endif //MINIJATI_OPTIMIZATORS_HPP
