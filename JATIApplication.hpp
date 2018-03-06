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
#ifndef MINIJATI_JATIAPPLICATION_HPP
#define MINIJATI_JATIAPPLICATION_HPP
// From the STL:
#include <string>
#include <map>
#include <Bpp/Exceptions.h>
#include <glog/logging.h>
#include <iostream>

namespace bpp {
    class JATIApplication {
    private:
        std::string appName_;
        mutable std::map<std::string, std::string> params_;
        bool timerStarted_;

    public:
        JATIApplication(int argc, char *argv[], const std::string &name);

    public:
        void startTimer();

        void done();

        std::map<std::string, std::string> &getParams() { return params_; }

        const std::string &getParam(const std::string &name) const {
            if (params_.find(name) == params_.end()) throw bpp::Exception("BppApplication::getParam(). Parameter '" + name + "' not found.");
            return params_[name];
        }

        std::string &getParam(const std::string &name) { return params_[name]; }


        void help() {
            std::cout << appName_ << std::endl << std::endl;
            std::cout << "Usage: miniJATI [options]" << std::endl;
            std::cout << std::endl << "**** Input options ****" << std::endl << std::endl;
            std::cout << "\talphabet={DNA,RNA,Protein}      Dataset alphabet     [requested]\t(inherited from bpp documentation)" << std::endl;
            std::cout << "\tinput_sequences=<filepath>                           [requested]" << std::endl;
            std::cout << "\tinput_tree=<filepath>                                [not-requested]" << std::endl;

            std::cout << std::endl << "**** Alignment options ****" << std::endl << std::endl;
            std::cout << "\talignment=<bool>                                     [requested]" << std::endl;

            std::cout << std::endl << "**** Substitution model options ****" << std::endl << std::endl;
            std::cout << "\tmodel_substitution=<string>                          [requested]\t(inherited from bpp documentation)" << std::endl;
            std::cout << "\tmodel_setfreqsfromdata=<bool>                        [requested]\t(inherited from bpp documentation)" << std::endl;
            std::cout << "\tmodel_pip_lambda=<float>                             (if models_indels=true)" << std::endl;
            std::cout << "\tmodel_pip_mu=<float>                                 (if models_indels=true)" << std::endl;

            std::cout << std::endl << "**** Likelihood computation options ****" << std::endl << std::endl;

            std::cout << std::endl << "**** Parameter optimisation options ****" << std::endl << std::endl;
            std::cout << "\toptimization=<string>                                        [requested]\t(inherited from bpp documentation)" << std::endl;
            std::cout << "\toptimization.message_handler=<filepath>                      [not-requested]\t(inherited from bpp documentation)" << std::endl;
            std::cout << "\toptimization.profiler=<filepath>                             [not-requested]\t(inherited from bpp documentation)" << std::endl;
            std::cout << "\toptimization.reparametrization=<bool>                        [not-requested]\t(inherited from bpp documentation)" << std::endl;
            std::cout << "\toptimization.max_number_f_eval=<int>                         [not-requested]\t(inherited from bpp documentation)" << std::endl;
            std::cout << "\toptimization.tolerance=<float>                               [not-requested]\t(inherited from bpp documentation)" << std::endl;
            std::cout << "\toptimization.final ={powell|simplex}                         [not-requested]\t(inherited from bpp documentation)" << std::endl;
            std::cout << "\toptimization.topology=<bool>                                 [not-requested]\t(inherited from bpp documentation)" << std::endl;
            std::cout << "\toptimization.topology.algorithm={greedy|hillclimbing|swarm}                                  [not-requested]" << std::endl;
            std::cout << "\toptimization.topology.algorithm.hillclimbing.startnodes=<int>                                [not-requested] number of random nodes to use during the"
                    " hillclimbing" << std::endl;
            std::cout << "\toptimization.topology.algorithm.operations={nni-search|spr-search|tbr-search|best-search}    [not-requested]" << std::endl;
            std::cout << "\toptimization.topology.algorithm.maxcycles=<int>                                              [not-requested]" << std::endl;
            std::cout << "\toptimization.topology.likelihood={single,double}                                             [not-requested]" << std::endl;

            std::cout << std::endl << "**** Output options ****" << std::endl << std::endl;
            std::cout << "\toutput.tree.file={path}                         The phylogenetic tree file to write to. " << std::endl;
            std::cout << "\toutput.tree.format={Newick|Nexus|NHX}           The format of the output tree file. " << std::endl;
            std::cout << "\toutput.trees.file={path}                        The file that will contain multiple trees. " << std::endl;
            std::cout << "\toutput.trees.format={Newick|Nexus|NHX}          The format of the output tree file." << std::endl;
            std::cout << "\toutput.infos={{path}|none}                      Alignment information log file (site specific rates, etc)" << std::endl;
            std::cout << "\toutput.estimates={{path}|none}                  Write parameter estimated values. " << std::endl;
            std::cout << "\toutput.estimates.alias={boolean}                Write the alias names of the aliased parameters instead of their values (default: true). " << std::endl;


            std::cout << std::endl << "Verbosity:" << std::endl;
            std::cout << "\texport GLOG_v={0,1,2,3}" << std::endl << "\texport GLOG_minloglevel={INFO,FATAL,WARN}" << std::endl << std::endl;
            std::cout << "Examples:" << std::endl;
            std::cout << "Documentation can be found at https://bitbucket.org/acg-team/minijati/" << std::endl;
        }


        void banner() {

            LOG(INFO) << "*****************************************************************************************************************************************";
            LOG(INFO) << "* " << appName_ << "  *";
            LOG(INFO) << "* Authors: Lorenzo Gatti & Massimo Maiolo                                                                                               *";
            LOG(INFO) << "* ------------------------------------------------------------------------------------------------------------------------------------- *";
            LOG(INFO) << "* Based on Bio++ by J. Dutheil, B. Boussau, L. Guéguen, M. Groussin                                                                     *";
            LOG(INFO) << "* Inspired on codonPhyML (Gil M. et al.)                                                                                                *";
            LOG(INFO) << "* Inspired on PrographMSA (Szalkowski A. et al.)                                                                                        *";
            LOG(INFO) << "* Implements the Poisson Indel Model (Bouchard-Cote A. et al.)                                                                          *";
            LOG(INFO) << "*****************************************************************************************************************************************";


        }

    };


} //end of namespace bpp;



#endif //MINIJATI_JATIAPPLICATION_HPP
