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
            std::cout << "Usage: miniJATI [options]" << std::endl;
            std::cout << std::endl << "### Alphabet options" << std::endl << std::endl;
            std::cout << "alphabet={DNA|RNA|Protein)|Codon(letter={DNA|RNA},type={Standard|EchinodermMitochondrial|InvertebrateMitochondrial|VertebrateMitochondrial})}" << std::endl;

            std::cout << "                                        The alphabet to use when reading sequences. DNA and RNA alphabet can in addition take" << std::endl;
            std::cout << "                                        an argument: " << std::endl;
            std::cout << "                                        bangAsgap={bool}" << std::endl;
            std::cout << "                                      Tell is exclamation mark should be considered as a gap character. The default " << std::endl;
            std::cout << "                                      is to consider it as an unknown character such as 'N' or '?'. " << std::endl;
            std::cout << "genetic_code={translation table}        Where ’translation table’ specifies the code to use, either as a text description," << std::endl;
            std::cout << "                                        or as the NCBI number. The following table give the currently implemented codes" << std::endl;
            std::cout << "                                        with their corresponding names: " << std::endl;
            std::cout << "                                        Standard                    1" << std::endl;
            std::cout << "                                        VertebrateMitochondrial     2" << std::endl;
            std::cout << "                                        YeastMitochondrial          3" << std::endl;
            std::cout << "                                        MoldMitochondrial           4" << std::endl;
            std::cout << "                                        InvertebrateMitochondrial   5" << std::endl;
            std::cout << "                                        EchinodermMitochondrial     9" << std::endl;
            std::cout << "                                        AscidianMitochondrial       13" << std::endl;
            std::cout << "                                        The states of the alphabets are in alphabetical order." << std::endl;

            std::cout << std::endl << "### [Input] Reading sequences" << std::endl << std::endl;

            std::cout << "input.sequence.file={path}                                The sequence file to use. (These sequences can also be not aligned). " << std::endl;
            std::cout << "input.sequence.format={format}                            The sequence file format. " << std::endl;
            std::cout << "input.sequence.sites_to_use={all|nogap|complete}          Tells which sites to use " << std::endl;
            std::cout << "input.sequence.remove_stop_codons={boolean}               Removes the sites where there is a stop codon (default: ’yes’)" << std::endl;
            std::cout << "input.sequence.max_gap_allowed=100%                       It specifies the maximum amount of gap allowed per site." << std::endl;
            std::cout << "input.sequence.max_unresolved_allowed=100%                It specifies the maximum amount of unresolved states per site." << std::endl;
            std::cout << "input.site.selection={list of integers}                   Will only consider sites in the given list of positions, in extended format :" << std::endl;
            std::cout << "                                                          positions separated with \",\", and \"i-j\" for all positions between i and j, " << std::endl;
            std::cout << "                                                          included." << std::endl;
            std::cout << "input.site.selection = {Sample(n={integer} [, replace={true}])}" << std::endl;
            std::cout << "                                                          Will consider {n} random sites, with optional replacement. " << std::endl;

            std::cout << std::endl << "The following formats are currently supported:" << std::endl;
            std::cout << "Fasta(extended={bool}, strictNames={bool})                The fasta format. The argument extended, default to 'no' " << std::endl;
            std::cout << "                                                          allows to enable the HUPO-PSI extension of the format." << std::endl;
            std::cout << "                                                          The argument strict_names, default to 'no', specifies" << std::endl;
            std::cout << "                                                          that only the first word in the fasta header is used as" << std::endl;
            std::cout << "                                                          a sequence names, the rest of the header being considered" << std::endl;
            std::cout << "                                                          as comments." << std::endl;
            std::cout << "Mase(siteSelection={chars})                               The Mase format (as read by Seaview and Phylo_win for instance)," << std::endl;
            std::cout << "                                                          with an optional site selection name. " << std::endl;
            std::cout << "Phylip(order={interleaved|sequential}, type={classic|extended}, split={spaces|tab})" << std::endl;
            std::cout << "                                                          The Phylip format, with several variations. " << std::endl;
            std::cout << "                                                          The argument order distinguishes between sequential and interleaved " << std::endl;
            std::cout << "                                                          format, while the option type distinguished between the plain old " << std::endl;
            std::cout << "                                                          Phylip format and the more recent extention allowing for sequence " << std::endl;
            std::cout << "                                                          names longer than 10 characters, as understood by PAML and PhyML." << std::endl;
            std::cout << "                                                          Finally, the split argument specifies the type of character that " << std::endl;
            std::cout << "                                                          separates the sequence name from the sequence content. " << std::endl;
            std::cout << "                                                          The conventional option is to use one (classic) or more (extended)" << std::endl;
            std::cout << "                                                          spaces, but tabs can also be used instead. " << std::endl;
            std::cout << "Clustal(extraSpaces={int})                                The Clustal format. " << std::endl;
            std::cout << "                                                          In its basic set up, sequence names do not have space characters," << std::endl;
            std::cout << "                                                          and one space splits the sequence content from its name. The parser" << std::endl;
            std::cout << "                                                          can however be configured to allow for spaces in the sequence names," << std::endl;
            std::cout << "                                                          providing a minimum number of space characters is used to split" << std::endl;
            std::cout << "                                                          the content from the name. Setting extraSpaces to 5 for instance, the" << std::endl;
            std::cout << "                                                          sequences are expected to be at least 6 spaces away for their names." << std::endl;
            std::cout << "Dcse()                                                    The DCSE alignment format. The secondary structure annotation will be ignored." << std::endl;
            std::cout << "Nexus()                                                   The Nexus alignment format. (Only very basic support is provided)" << std::endl;
            std::cout << "GenBank()                                                 The GenBank not aligned sequences format. " << std::endl;
            std::cout << "                                                          Very basic support: only retrieves the sequence content for now, " << std::endl;
            std::cout << "                                                          all features are ignored." << std::endl;


            std::cout << std::endl << "### [Input] Reading trees" << std::endl << std::endl;

            std::cout << "input.tree.file={path}                      The phylogenetic tree file to use." << std::endl;
            std::cout << "input.tree.format={Newick|Nexus|NHX}        The format of the input tree file." << std::endl;

            std::cout << std::endl << "**** Alignment options ****" << std::endl << std::endl;
            std::cout << "alignment=<bool>                                     [requested]" << std::endl;

            std::cout << std::endl << "**** Branch lengths initial values ****" << std::endl << std::endl;
            std::cout << "init.tree={user|random|distance}                  Set the method for the initial tree to use. " << std::endl;
            std::cout << "                                                  The user option allows you to use an existing file passed via input.tree.file" << std::endl;
            std::cout << "                                                  This file may have been built using another method like neighbor joining or " << std::endl;
            std::cout << "                                                  parsimony for instance. The random option picks a random tree, which is handy " << std::endl;
            std::cout << "                                                  to test convergence.  This may however slows down significantly the optimization  " << std::endl;
            std::cout << "                                                  process. " << std::endl;
            std::cout << "init.distance.matrix.file={path}                  A distance matrix can be supplied instead of being computed from the alignment." << std::endl;
            std::cout << "init.distance.method={wpgma|upgma|nj|bionj}       When distance method is required, the user can specify which algorithm to use." << std::endl;
            std::cout << "init.brlen.method={method description}            Set how to initialize the branch lengths. Available methods include:" << std::endl;
            std::cout << std::endl;
            std::cout << "Input(midpoint_root_branch={boolean})             Keep initial branch lengths as is. Additional argument specifies if the root " << std::endl;
            std::cout << "                                                  position should be moved to the midpoint position of the branch containing it. " << std::endl;
            std::cout << std::endl;
            std::cout << "Equal(value={float>0})                            Set all branch lengths to the same value, provided as argumemt. " << std::endl;
            std::cout << std::endl;
            std::cout << "Clock                                             Coerce to a clock tree." << std::endl;
            std::cout << std::endl;
            std::cout << "Grafen(height={{real>0}|input}, rho = {real>0})   Uses Grafen’s method to compute branch lengths." << std::endl;
            std::cout << "                                                  In Grafen’s method, each node is given a weight equal to the number of underlying  " << std::endl;
            std::cout << "                                                  leaves. The length of each branch is then computed as the difference of the weights" << std::endl;
            std::cout << "                                                  of the connected nodes, and further divided by the number of leaves in the tree. " << std::endl;
            std::cout << "                                                  The height of all nodes are then raised to the power of ’rho’, a user specified value. " << std::endl;
            std::cout << "                                                  The tree is finally scaled to match a given total height, which can be the original " << std::endl;
            std::cout << "                                                  one (height=input), or fixed to a certain value (usually height=1). A value of " << std::endl;
            std::cout << "                                                  rho=0 provides a star tree, and the greater the value of rho, the more recent the" << std::endl;
            std::cout << "                                                  inner nodes. " << std::endl;
            std::cout << std::endl;
            std::cout << "input.tree.check_root = {boolean}                 Tell if the input tree should be checked regarding to the presence of a root." << std::endl;

            std::cout << std::endl << "### Substitution model options" << std::endl << std::endl;
            std::cout << "model=<string>                                      A description of the substitution model to use, using the keyval syntax. " << std::endl;
            //std::cout << "model_setfreqsfromdata=<bool>                       [requested](inherited from bpp documentation)" << std::endl;
            //std::cout << "model_pip_lambda=<float>                            (if models_indels=true)" << std::endl;
            //std::cout << "model_pip_mu=<float>                                (if models_indels=true)" << std::endl;

            std::cout << std::endl << "The following nucleotide models are currently available:" << std::endl << std::endl;
            std::cout << "JC69 " << std::endl;
            std::cout << "K80([kappa={real>0}]) " << std::endl;
            std::cout << "F84([kappa={real>0}, theta={real]0,1[}, theta1={real]0,1[},theta2={real]0,1[} ,\"equilibrium frequencies\"])" << std::endl;
            std::cout << "HKY85([kappa={real>0}, theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[} ,\"equilibrium frequencies\"])" << std::endl;
            std::cout << "T92([kappa={real>0}, theta={real]0,1[} ,\"equilibrium frequencies\"])" << std::endl;
            std::cout << "TN93([kappa1={real>0}, kappa2={real>0}, theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[} ,\"equilibrium frequencies\"])" << std::endl;
            std::cout << "GTR([a={real>0}, b={real>0}, c={real>0}, d={real>0}, e={real>0}, theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[} ,\"equilibrium frequencies\"])" << std::endl;
            std::cout << "L95([beta={real>0}, gamma={real>0}, delta={real>0}, theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[} ,\"equilibrium frequencies\"])" << std::endl;
            std::cout << "SSR([beta={real>0}, gamma={real>0}, delta={real>0}, theta={real]0,1[}])" << std::endl;
            std::cout << "RN95([thetaR={real]0,1[}, thetaC={real]0,1[}, thetaG={real]0,1[}, kappaP={real[0,1[}, gammaP={real[0,1[}, sigmaP={real>1}, alphaP={real>1}])" << std::endl;
            std::cout << "RN95s([thetaA={real]0,0.5[}, gamma={real]0,0.5[}, alphaP={real>1}])" << std::endl;

            std::cout << std::endl << "The following protein models are currently available:" << std::endl << std::endl;
            std::cout << "JC69 " << std::endl;
            std::cout << "DSO78 " << std::endl;
            std::cout << "JTT92" << std::endl;
            std::cout << "WAG01" << std::endl;
            std::cout << "LG08" << std::endl;
            std::cout << "LLG08_EX2([relrate1={real]0,1[}, relproba1={real]0,1[}])" << std::endl;
            std::cout << "LLG08_EX3([relrate1={real]0,1[}, relrate2={real]0,1[}, relproba1={real]0,1[}, relproba2={real]0,1[}])" << std::endl;
            std::cout << "LLG08_EHO([relrate1={real]0,1[}, relrate2={real]0,1[}, relproba1={real]0,1[}, relproba2={real]0,1[}])" << std::endl;
            std::cout << "LLG08_UL2([relrate1={real]0,1[}, relproba1={real]0,1[}])" << std::endl;
            std::cout << "LLG08_UL3([relrate1={real]0,1[}, relrate2={real]0,1[}, relproba1={real]0,1[}, relproba2={real]0,1[}])" << std::endl;
            std::cout << "LGL08_CAT(nbCat={[10,20,30,40,50,60]}, [relrate1={real]0,1[}, relrate2={real]0,1[}, ..., relproba1={real]0,1[}, relproba2={real]0,1[}, ...] ))" << std::endl;
            std::cout << "LGL08_CAT_C{[1,...,nbCat]}(nbCat={[10,20,30,40,50,60]})" << std::endl;
            std::cout << "DSO78+F([theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[}, ... ,\"equilibrium frequencies\"])" << std::endl;
            std::cout << "JTT92+F([theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[}, ..., \"equilibrium frequencies\"])" << std::endl;
            std::cout << "WAG01+F([theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[}, ..., \"equilibrium frequencies\"])" << std::endl;
            std::cout << "LG08+F([theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[}, ..., \"equilibrium frequencies\"])" << std::endl;
            std::cout << "Empirical(name={chars}, file={path})" << std::endl;
            std::cout << "Empirical+F(name={chars}, file={path}, [theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[}, ...,  \"equilibrium frequencies\"])" << std::endl;

            std::cout << std::endl << "The following meta models are currently available:" << std::endl << std::endl;
            std::cout << "PIP(model={model description}, lambda={real>0}, mu={real>0} [, \"equilibrium frequencies\"]) " << std::endl;
            std::cout << "TS98(model={model description}, s1={real>0}, s2={real>0} [, \"equilibrium frequencies\"])" << std::endl;
            std::cout << "G01(model={model description}, rdist={rate distribution description}, mu={real>0} [, \"equilibrium frequencies\"])" << std::endl;
            std::cout << "RE08(model={model description}, lambda={real>0}, mu={real>0} [, \"equilibrium frequencies\"])" << std::endl;

            std::cout << std::endl << "### Frequencies distribution sets" << std::endl << std::endl;
            std::cout << "The following frequencies distributions are available:" << std::endl;
            std::cout << "Fixed()                                                                  All frequencies are fixed to their initial value and are not estimated. " << std::endl;
            std::cout << std::endl;
            std::cout << "GC(theta={real]0,1[})                                                    For nucleotides only, set the G content equal to the C content. " << std::endl;
            std::cout << std::endl;
            std::cout << "Full(theta1={real]0,1[}, theta2={real]0,1[}, ..., thetaN={real]0,1[})    Full parametrization. Contains N free parameters." << std::endl;
            std::cout << std::endl;

            std::cout << std::endl << "### Rate across site distribution" << std::endl << std::endl;
            std::cout << "rate_distribution={rate distribution description}                Specify the rate across sites distribution" << std::endl;

            std::cout << std::endl << "The following distributions are currently available:" << std::endl << std::endl;

            std::cout << "[Most used distributions]" << std::endl;
            std::cout << "Constant                                                         Uses a constant rate across sites" << std::endl;
            std::cout << std::endl;
            std::cout << "Gamma(n={int>=2}, alpha={float>0})                               A discretized gamma distribution of rates, with n classes, and a given shape," << std::endl;
            std::cout << "                                                                 with mean 1 (scale=shape)." << std::endl;
            std::cout << std::endl;
            std::cout << "Invariant(dist={rate distribution description}, p={real[0,1]})   A composite distribution allowing a special class of invariant site, with a " << std::endl;
            std::cout << "                                                                 probability p." << std::endl;
            std::cout << std::endl;
            std::cout << "[Standard distributions]" << std::endl;
            std::cout << "Beta(n={int>=2}, alpha={float>0}, beta={float>0})                A discretized beta distribution, with n classes, with standard parameters" << std::endl;
            std::cout << "                                                                 alpha and beta." << std::endl;
            std::cout << std::endl;
            std::cout << "Gamma(n={int>=2}, alpha={float>0}, beta={float>0})               A discretized gamma distribution, with n classes, a shape alpha and a rate " << std::endl;
            std::cout << "                                                                 beta, with parameters alpha and beta. " << std::endl;
            std::cout << std::endl;
            std::cout << "Gaussian(n={int>=1}, mu={float}, sigma={float>0})                a discretized gaussian distribution, with n classes, a mean mu and a " << std::endl;
            std::cout << "                                                                 standard deviation sigma, with parameters mu and sigma. " << std::endl;
            std::cout << std::endl;
            std::cout << "Exponential(n={int>=2}, lambda={float>0})                        a discretized exponential distribution, with n classes and parameter lambda." << std::endl;
            std::cout << std::endl;
            std::cout << "Simple(values={vector<double>}, probas={vector<double>} [, ranges={vector<parametername[min;val]>}])" << std::endl;
            std::cout << "                                                                 A discrete distribution with specific values (in values) and their respective " << std::endl;
            std::cout << "                                                                 non-negative probabibilities (in probas). The parameters are V1, V2, ..., Vn" << std::endl;
            std::cout << "                                                                 for all the values and the relative probabibility parameters are " << std::endl;
            std::cout << "                                                                 theta1, theta2, ..., thetan-1." << std::endl;
            std::cout << "                                                                 Optional argument {ranges} sets the allowed ranges of values taken by the " << std::endl;
            std::cout << "                                                                 parameters; usage is like 'ranges=(V1[0.2;0.9],V2[1.1;999])'." << std::endl;
            std::cout << std::endl;
            std::cout << "TruncExponential(n={int>=2}, lambda={float>0}, tp={float>0})     A discretized truncated exponential distribution, with n classes, parameter " << std::endl;
            std::cout << "                                                                 lambda and a truncation point tp. The parameters are lambda and tp. " << std::endl;
            std::cout << std::endl;
            std::cout << "Uniform(n={int>=1}, begin={float>0}, end={float>0})              A uniform distribution, with n classes in interval [begin,end]." << std::endl;
            std::cout << "                                                                 There are no parameters." << std::endl;
            std::cout << std::endl;
            std::cout << "[Mixture Distributions]" << std::endl;
            std::cout << "Invariant(dist={distribution description}, p={float>0})          A Mixture of a given discrete distributution and a 0 Dirac. p is the  " << std::endl;
            std::cout << "                                                                 probability of this 0 Dirac." << std::endl;
            std::cout << std::endl;
            std::cout << "Mixture(probas={vector<double>}, dist1={distribution description}, ..., distn={distribution description})" << std::endl;
            std::cout << "                                                                 A Mixture of discrete distributions with specific probabilities (in probas)" << std::endl;
            std::cout << "                                                                 and their respective desccriptions (in probas). The parameters are the relative " << std::endl;
            std::cout << "                                                                 probabibility parameters theta1, theta2, ..., thetan-1, and the parameters" << std::endl;
            std::cout << "                                                                 of the included distributions prefixed by Mixture.i_ where i is the order " << std::endl;
            std::cout << "                                                                 of the distribution. " << std::endl;
            std::cout << std::endl;

            std::cout << std::endl << "### Likelihood computation options" << std::endl << std::endl;

            std::cout << std::endl << "### Numerical parameter optimisation options" << std::endl << std::endl;
            std::cout << "This program allows to (re-)estimate numerical parameters, including\n"
                    "- Branch lengths\n"
                    "- Entries of the substitution matrices, included base frequencies values)\n"
                    "- Parameters of the rate distribution (currently shape parameter of the gamma law, proportion of invariant sites)." << std::endl << std::endl;

            std::cout << "optimization={method}" << std::endl;

            std::cout << std::endl << "The following methods are currently available:" << std::endl << std::endl;
            std::cout << "None                                                         No optimization is performed, initial values are kept 'as is'." << std::endl;
            std::cout << std::endl;
            std::cout << "FullD(derivatives={Newton|Gradient})                         Full-derivatives method." << std::endl;
            std::cout << "                                                             Branch length derivatives are computed analytically, others numerically." << std::endl;
            std::cout << "                                                             The derivatives arguments specifies if first or second order derivatives" << std::endl;
            std::cout << "                                                             should be used. In the first case, the optimization method used is the  " << std::endl;
            std::cout << "                                                             so-called conjugate gradient method, otherwise the Newton-Raphson method " << std::endl;
            std::cout << "                                                             will be used." << std::endl;
            std::cout << std::endl;
            std::cout << "D-Brent(derivatives={Newton|Gradient|BFGS}, nstep={int>0})   Branch lengths parameters are optimized using either the conjugate gradient, " << std::endl;
            std::cout << "                                                             the Newton-Raphson method or BFGS, other parameters are estimated using the " << std::endl;
            std::cout << "                                                             Brent method in one dimension. The algorithm then loops over all parameters" << std::endl;
            std::cout << "                                                             until convergence. The nstep arguments allow to specify a number of " << std::endl;
            std::cout << "                                                             progressive steps to perform during optimization. If nstep=3 and " << std::endl;
            std::cout << "                                                             precision=E-6, a first optimization with precision=E-2, will be performed, " << std::endl;
            std::cout << "                                                             then a round with precision set to E-4 and finally precision will be set to " << std::endl;
            std::cout << "                                                             E-6. This approach generally increases convergence time. " << std::endl;
            std::cout << std::endl;
            std::cout << "D-BFGS(derivatives={Newton|Gradient|BFGS}, nstep={int>0})    Branch lengths parameters are optimized using either the conjugate gradient,  " << std::endl;
            std::cout << "                                                             the Newton-Raphson method or BFGS, other parameters are estimated using the " << std::endl;
            std::cout << "                                                             BFGS method in one dimension. The algorithm then loops over all parameters" << std::endl;
            std::cout << "                                                             until convergence. The nstep arguments allow to specify a number of " << std::endl;
            std::cout << "                                                             progressive steps to perform during optimization. If nstep=3 and " << std::endl;
            std::cout << "                                                             precision=E-6, a first optimization with precision=E-2, will be performed, " << std::endl;
            std::cout << "                                                             then a round with precision set to E-4 and finally precision will be set to " << std::endl;
            std::cout << "                                                             E-6. This approach generally increases convergence time. " << std::endl;
            std::cout << std::endl;
            std::cout << "optimization.reparametrization=<bool>                        Tells if parameters should be transformed in order to remove constraints (for" << std::endl;
            std::cout << "                                                             instance positivie-only parameters will be log transformed in order to obtain" << std::endl;
            std::cout << "                                                             parameters defined from -inf to +inf). This may improve the optimization, " << std::endl;
            std::cout << "                                                             particularly for parameter-rich models, but the likelihood calculations " << std::endl;
            std::cout << "                                                             will take a bit more time." << std::endl;
            std::cout << "optimization.final={powell|simplex|bfgs}                     Optional final optimization step, useful if numerical derivatives are " << std::endl;
            std::cout << "                                                             to be used. Leave the field empty in order to skip this step." << std::endl;
            std::cout << "optimization.profiler={{path}|std|none}                      A file where to dump optimization steps (a file path or std for standard " << std::endl;
            std::cout << "                                                             output or none for no output)." << std::endl;
            std::cout << "optimization.message_handler={{path}|std|none}               A file where to dump warning messages." << std::endl;
            std::cout << "optimization.max_number_f_eval=<int>0>                       The maximum number of likelihood evaluations to perform." << std::endl;
            std::cout << "optimization.tolerance=<float>0>                             The precision on the log-likelihood to reach." << std::endl;
            std::cout << "optimization.constrain_parameter={list<chars=interval>}      A list of parameters on which the authorized values are limited to a given " << std::endl;
            std::cout << "                                                             interval. For example:" << std::endl;
            std::cout << "                                                             optimization.constrain_parameter = YN98.omega = [-inf;1.9[, *theta* = [0.1;0.7[, BrLen*=[0.01;inf]" << std::endl;
            std::cout << "optimization.ignore_parameter={list<chars>}                  A list of parameters to ignore during the estimation process. The parameter" << std::endl;
            std::cout << "                                                             name should include there 'namespace', that is their model name, for " << std::endl;
            std::cout << "                                                             instance K80.kappa,TN93.theta, GTR.a, Gamma.alpha, etc.'BrLen' will " << std::endl;
            std::cout << "                                                             ignore all branch lengths and 'Model' will ignore all model parameters." << std::endl;
            std::cout << "                                                             The ’*’ wildcard can be used, as in *theta* for all the parameters whose " << std::endl;
            std::cout << "                                                             name has theta in it." << std::endl;


            std::cout << std::endl << "### Topology optimisation options" << std::endl << std::endl;

            std::cout << "optimization.topology=<bool>                                                                Enable the tree topology estimation" << std::endl;
            std::cout << "optimization.topology.algorithm={greedy|hillclimbing|swarm}                                 Algorithm to use for topology estimation" << std::endl;
            std::cout << "optimization.topology.algorithm.operations={nni-search|spr-search|tbr-search|best-search}   Tree rearrangment operations to apply during tree-search" << std::endl;
            std::cout << "optimization.topology.algorithm.maxcycles=<int>                                             Max number of tree-search cycles" << std::endl;
            std::cout << "optimization.topology.likelihood={single,double}                                            Method for recomputing the likelihood after tree move" << std::endl;
            std::cout << "optimization.topology.algorithm.hillclimbing.startnodes=<int>                               Number of starting random nodes to use during hc" << std::endl;

            std::cout << std::endl << "### Output options" << std::endl << std::endl;

            std::cout << "output.tree.file={path}                         The phylogenetic tree file to write to. " << std::endl;
            std::cout << "output.tree.format={Newick|Nexus|NHX}           The format of the output tree file. " << std::endl;
            std::cout << "output.trees.file={path}                        The file that will contain multiple trees. " << std::endl;
            std::cout << "output.trees.format={Newick|Nexus|NHX}          The format of the output tree file." << std::endl;
            std::cout << "output.infos={{path}|none}                      Alignment information log file (site specific rates, etc)" << std::endl;
            std::cout << "output.estimates={{path}|none}                  Write parameter estimated values. " << std::endl;
            std::cout << "output.estimates.alias={boolean}                Write the alias names of the aliased parameters instead of their values (default: true). " << std::endl;

            std::cout << std::endl << "###  Verbosity" << std::endl << std::endl;
            std::cout << "export GLOG_v={0,1,2,3}" << std::endl << "export GLOG_minloglevel={INFO,FATAL,WARN}" << std::endl << std::endl;

            std::cout << std::endl << "### Examples" << std::endl << std::endl;

            std::cout << "Documentation can be found at https://bitbucket.org/acg-team/minijati/" << std::endl;
        }


        void banner() {

            auto host_name = boost::asio::ip::host_name();

            bpp::ApplicationTools::displayMessage("*****************************************************************************************************************************************");
            bpp::ApplicationTools::displayMessage("* " + appName_ + " by Lorenzo Gatti & Massimo Maiolo                                                                                      *");
            bpp::ApplicationTools::displayMessage("* Build on commit: " + appVersion_ + " on date: " + appBuild_ + "                           *");
            bpp::ApplicationTools::displayMessage("*****************************************************************************************************************************************");
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
