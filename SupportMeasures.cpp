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
 * @file SupportMeasures.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 02 05 2018
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

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include "SupportMeasures.hpp"
#include "UnifiedTSHomogeneousTreeLikelihood_PIP.hpp"

using namespace bpp;

Boostrap::Boostrap(AbstractHomogeneousTreeLikelihood *tl,
                   const SiteContainer &data,
                   DiscreteDistribution *rDist,
                   tshlib::Utree *utree,
                   UtreeBppUtils::treemap *tm,
                   std::map<std::string, std::string>& params) {

    unsigned int nbBS = ApplicationTools::getParameter<unsigned int>("bootstrap.number", params, 0, "", true, 1);

    if (nbBS > 0) {
        ApplicationTools::displayResult("Number of bootstrap samples", TextTools::toString(nbBS));
        bool approx = ApplicationTools::getBooleanParameter("bootstrap.approximate", params, true, "", true, 2);
        ApplicationTools::displayBooleanResult("Use approximate bootstrap", approx);
        bool bootstrapVerbose = ApplicationTools::getBooleanParameter("bootstrap.verbose", params, false, "", true, 2);

        const Tree *initTree = &tl->getTree();

        std::string bsTreesPath = ApplicationTools::getAFilePath("bootstrap.output.file", params, false, false);
        std::ofstream *out = nullptr;

        if (bsTreesPath != "none") {
            ApplicationTools::displayResult("Bootstrap trees stored in file", bsTreesPath);
            out = new std::ofstream(bsTreesPath.c_str(), std::ios::out);
        }

        Newick newick;
        ParameterList paramsToIgnore = tl->getSubstitutionModelParameters();
        paramsToIgnore.addParameters(tl->getRateDistributionParameters());

        ApplicationTools::displayTask("Bootstrapping", true);
        std::vector<Tree *> bsTrees(nbBS);
        for (unsigned int i = 0; i < nbBS; i++) {
            ApplicationTools::displayGauge(i, nbBS - 1, '=');
            VectorSiteContainer *sample = SiteContainerTools::bootstrapSites(*tl->getData());
            if (!approx) {
                tl->getModel()->setFreqFromData(*sample);
            }

            AbstractHomogeneousTreeLikelihood *tlRep;
            //auto *tlRep = new UnifiedTSHomogeneousTreeLikelihood_PIP(initTree, sample, nullptr, nullptr, utree, &tm, true, params, "", false, false, false);


            tlRep->initialize();
            ParameterList parametersRep = tlRep->getParameters();
            if (approx) {
                parametersRep.deleteParameters(paramsToIgnore.getParameterNames());
            }

            //tl = dynamic_cast<AbstractHomogeneousTreeLikelihood *>(Optimizer::optimizeParameters(tlRep, nullptr, tlRep->getParameters(), params, "", true, true, 0));

            bsTrees[i] = new TreeTemplate<Node>(tlRep->getTree());
            if (out && i == 0) newick.write(*bsTrees[i], bsTreesPath, true);
            if (out && i > 0) newick.write(*bsTrees[i], bsTreesPath, false);
            delete tlRep;
            delete sample;
        }
        if (out) out->close();
        if (out) delete out;
        ApplicationTools::displayTaskDone();


        ApplicationTools::displayTask("Compute bootstrap values");
        //TreeTools::computeBootstrapValues(tl->getTree(), bsTrees);
        ApplicationTools::displayTaskDone();
        for (unsigned int i = 0; i < nbBS; i++) {
            delete bsTrees[i];
        }

    // Write resulting tree:
        PhylogeneticsApplicationTools::writeTree(*bestTree, params);
    }

}

