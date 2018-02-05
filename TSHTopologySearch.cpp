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
 * @file TSHTopologySearch.cpp
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
#include "TSHTopologySearch.hpp"


using namespace bpp;


void TSHTopologySearch::notifyAllPerformed(const TopologyChangeEvent &event) {
    searchableTree_->topologyChangePerformed(event);
    for (size_t i = 0; i < topoListeners_.size(); i++) {
        topoListeners_[i]->topologyChangePerformed(event);
    }
}

void TSHTopologySearch::notifyAllTested(const TopologyChangeEvent &event) {
    searchableTree_->topologyChangeTested(event);
    for (size_t i = 0; i < topoListeners_.size(); i++) {
        topoListeners_[i]->topologyChangeTested(event);
    }
}

void TSHTopologySearch::notifyAllSuccessful(const TopologyChangeEvent &event) {
    searchableTree_->topologyChangeSuccessful(event);
    for (size_t i = 0; i < topoListeners_.size(); i++) {
        topoListeners_[i]->topologyChangeSuccessful(event);
    }
}

void TSHTopologySearch::search() throw(Exception) {

}