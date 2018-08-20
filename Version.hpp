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
 * @file Version.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 05 01 2018
 * @version 1.0.10
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
 * @see For more information visit: https://bitbucket.org/acg-team/minijati/wiki/Home
 */
#ifndef MINIJATI_VERSION_HPP
#define MINIJATI_VERSION_HPP



#include <string>

namespace software{

    //version 1.0.0 (master 1a8e6107, 16 Jun 2013, 17:50:26)
    std::string version(PRJ_VERSION);
    std::string releasegitbranch(PRJ_GITBRANCH);
    std::string releasegitref(PRJ_GITREF);
    std::string releasedate(PRJ_DATE);
    std::string releasetime(PRJ_TIME);

    std::string build = version +  " (" +releasegitbranch + " " + releasegitref + ", "+ releasedate + ", " + releasetime + ")";
    std::string name(PRJ_NAME);
    std::string name_extended(PRJ_DESC);

    std::string desc = name_extended + " ("+ name + ") " + build;
}
#endif //MINIJATI_VERSION_HPP


