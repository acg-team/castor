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
 * @file TSHLIBSearchable.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 05 02 2018
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
#ifndef MINIJATI_TSHLIBSEARCHABLE_HPP
#define MINIJATI_TSHLIBSEARCHABLE_HPP

#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TopologySearch.h>

namespace bpp {

    class TSHSearchable : public TopologyListener, public virtual Clonable {
    public:
        TSHSearchable() {}

        virtual ~TSHSearchable() {}

        virtual TSHSearchable *clone() const = 0;

    public:

        /**
         * @brief Send the score of a TSH movement, without performing it.
         *
         * This methods sends the score variation.
         * This variation must be negative if the new point is better,
         * i.e. the object is to be used with a minimizing optimization
         * (for consistence with Optimizer objects).
         *
         * @param nodeId The id of the node defining the TSH movement.
         * @return The score variation of the TSH.
         * @throw NodeException If the node does not define a valid TSH.
         */
        virtual double testTSHRearrangement(int nodeId) const throw(NodeException) = 0;

        /**
         * @brief Perform a TSH movement.
         *
         * @param nodeId The id of the node defining the TSH movement.
         * @throw NodeException If the node does not define a valid TSH.
         */
        virtual void doTSHRearrangment(int nodeId) throw(NodeException) = 0;

        /**
         * @brief Get the tree associated to this TSHSearchable object.
         *
         * @return The tree associated to this instance.
         */
        virtual const Tree &getTopology() const = 0;

        /**
         * @brief Get the current score of this TSHSearchable object.
         *
         * @return The current score of this instance.
         */
        virtual double getTopologyValue() const throw(Exception) = 0;

    };

} //end of namespace bpp.
#endif //MINIJATI_TSHLIBSEARCHABLE_HPP
