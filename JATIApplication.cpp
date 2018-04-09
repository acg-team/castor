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
 * @file JATIApplication.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 06 02 2018
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
#include "JATIApplication.hpp"
#include <iostream>
#include <Bpp/Utils/AttributesTools.h>
#include <glog/logging.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>

using namespace bpp;

JATIApplication::JATIApplication(int argc, char *argv[], const std::string &name, const std::string &strVersion, const std::string &build_date) :
        appName_(name), appBuild_(build_date), appVersion_(strVersion), params_(), timerStarted_(false) {
    //LOG(INFO) << "Parsing options:";
    params_ = bpp::AttributesTools::parseOptions(argc, argv);
    bool showversion = bpp::ApplicationTools::getBooleanParameter("version", params_, false, "", true, 3);
    bpp::ApplicationTools::warningLevel = bpp::ApplicationTools::getIntParameter("--warning", params_, 0, "", true, 3);
    bool noint = bpp::ApplicationTools::getBooleanParameter("--noninteractive", params_, false, "", true, 3);
    bpp::ApplicationTools::interactive = !noint;
    long seed = bpp::ApplicationTools::getParameter<long>("--seed", params_, -1, "", true, 3);
    if (seed >= 0) {
        bpp::RandomTools::setSeed(seed);
        bpp::ApplicationTools::displayResult("Random seed set to", seed);
    }
    if (showversion) {
        this->version();
        exit(0);
    }
}

void JATIApplication::startTimer() {
    ApplicationTools::startTimer();
    timerStarted_ = true;
}

void JATIApplication::done() {
    LOG(INFO) << appName_ << "'s done. Bye.";
    if (timerStarted_)
        bpp::ApplicationTools::displayTime("Total execution time:");
}
