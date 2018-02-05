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
 * @file TSHTopologySearch.hpp
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
#ifndef MINIJATI_TSHTOPOLOGYSEARCH_HPP
#define MINIJATI_TSHTOPOLOGYSEARCH_HPP

#include <Bpp/Phyl/TopologySearch.h>
#include "TSHSearchable.hpp"

namespace bpp {

/**
 * @brief TSH topology search method.
 *
 */
    class TSHTopologySearch : public virtual TopologySearch {
    public:
        const static std::string FAST;
        const static std::string BETTER;
        const static std::string PHYML;

    private:
        TSHSearchable *searchableTree_;
        std::string algorithm_;
        unsigned int verbose_;
        std::vector<TopologyListener *> topoListeners_;

    public:
        TSHTopologySearch(
                TSHSearchable &tree,
                const std::string &algorithm = FAST,
                unsigned int verbose = 2) :
                searchableTree_(&tree), algorithm_(algorithm), verbose_(verbose), topoListeners_() {}

        /*
        TSHTopologySearch(const TSHSearchable& ts) :
                searchableTree_(ts.searchableTree_),
                algorithm_(ts.algorithm_),
                verbose_(ts.verbose_),
                topoListeners_(ts.topoListeners_)
        {
            //Hard-copy all listeners:
            for (unsigned int i = 0; i < topoListeners_.size(); i++)
                topoListeners_[i] = dynamic_cast<TopologyListener*>(ts.topoListeners_[i]->clone());
        }
        */
        TSHTopologySearch &operator=(const TSHTopologySearch &ts) {
            searchableTree_ = ts.searchableTree_;
            algorithm_ = ts.algorithm_;
            verbose_ = ts.verbose_;
            topoListeners_ = ts.topoListeners_;
            //Hard-copy all listeners:
            for (unsigned int i = 0; i < topoListeners_.size(); i++)
                topoListeners_[i] = dynamic_cast<TopologyListener *>(ts.topoListeners_[i]->clone());
            return *this;
        }


        virtual ~TSHTopologySearch() {
            for (std::vector<TopologyListener *>::iterator it = topoListeners_.begin();
                 it != topoListeners_.end();
                 it++)
                delete *it;
        }

    public:
        void search() throw(Exception);

        /**
         * @brief Add a listener to the list.
         *
         * All listeners will be notified in the order of the list.
         * The first listener to be notified is the NNISearchable object itself.
         *
         * The listener will be owned by this instance, and copied when needed.
         */
        void addTopologyListener(TopologyListener *listener) {
            if (listener)
                topoListeners_.push_back(listener);
        }

    public:
        /**
         * @brief Retrieve the tree.
         *
         * @return The tree associated to this instance.
         */
        const Tree &getTopology() const { return searchableTree_->getTopology(); }

        /**
         * @return The NNISearchable object associated to this instance.
         */
        TSHSearchable *getSearchableObject() { return searchableTree_; }

        /**
         * @return The NNISearchable object associated to this instance.
         */
        const TSHSearchable *getSearchableObject() const { return searchableTree_; }

    protected:
        //void searchFast()   throw (Exception);
        //void searchBetter() throw (Exception);
        //void searchPhyML()  throw (Exception);

        /**
         * @brief Process a TopologyChangeEvent to all listeners.
         */
        void notifyAllPerformed(const TopologyChangeEvent &event);

        /**
         * @brief Process a TopologyChangeEvent to all listeners.
         */
        void notifyAllTested(const TopologyChangeEvent &event);

        /**
         * @brief Process a TopologyChangeEvent to all listeners.
         */
        void notifyAllSuccessful(const TopologyChangeEvent &event);

    };

} //end of namespace bpp.


#endif //MINIJATI_TSHTOPOLOGYSEARCH_HPP
