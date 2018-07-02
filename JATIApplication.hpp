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
 * @file JATIApplication.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 06 02 2018
 * @version 1.0.7
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
#ifndef MINIJATI_JATIAPPLICATION_HPP
#define MINIJATI_JATIAPPLICATION_HPP
// From the STL:
#include <string>
#include <map>
#include <Bpp/Exceptions.h>
#include <glog/logging.h>
#include <iostream>
#include <Bpp/App/ApplicationTools.h>

#include <boost/asio/ip/host_name.hpp>


namespace bpp {
    class JATIApplication {
    private:
        std::string appName_;
        std::string appBuild_;
        std::string appVersion_;
        mutable std::map<std::string, std::string> params_;
        bool timerStarted_;
        long seed_;

    public:
        JATIApplication(int argc, char *argv[], const std::string &name, const std::string &strVersion, const std::string &build_date);

    public:
        void startTimer();

        void done();

        std::map<std::string, std::string> &getParams() { return params_; }

        const std::string &getParam(const std::string &name) const {
            if (params_.find(name) == params_.end()) throw bpp::Exception("BppApplication::getParam(). Parameter '" + name + "' not found.");
            return params_[name];
        }

        std::string &getParam(const std::string &name) { return params_[name]; }

        long getSeed() {return seed_;}

        void help() {
            std::cout << appName_ << std::endl << std::endl;
            std::cout << "Usage: miniJATI [arguments] or [params=file.txt]" << std::endl;
            std::cout << "Documentation can be found at https://bitbucket.org/acg-team/minijati/" << std::endl;
        }


        void banner() {

            auto host_name = boost::asio::ip::host_name();

            bpp::ApplicationTools::displayMessage("------------------------------------------------------------------------------");
            bpp::ApplicationTools::displayMessage(appName_);
            bpp::ApplicationTools::displayMessage("Phylogenetic Tree Inference and Multiple Sequence Alignment under Indel models");
            bpp::ApplicationTools::displayMessage("Authors: Lorenzo Gatti & Massimo Maiolo");
            bpp::ApplicationTools::displayMessage("Build on commit: " + appVersion_);
            bpp::ApplicationTools::displayMessage("On date: "+ appBuild_);
            bpp::ApplicationTools::displayMessage("------------------------------------------------------------------------------");
            bpp::ApplicationTools::displayResult("Execution started on:", host_name);




        }

        void version() {
            std::cout << appName_ << std::endl;
            std::cout << appVersion_ << std::endl;
            std::cout << appBuild_ << std::endl;
        }

    };


} //end of namespace bpp;



#endif //MINIJATI_JATIAPPLICATION_HPP
