/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti & Massimo Maiolo
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
 * @file main.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 20 12 2017
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

#include <iostream>

<<<<<<< HEAD
int main() {
    std::cout << "Hello, World!" << std::endl;
    std::cout << "Hello, World!" << std::endl;
    std::cout << "Hello, World!" << std::endl;
=======
/*
* From Core:
*/
#include <Bpp/Io/OutputStream.h>

/*
* From SeqLib:
*/
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
/*
* From PhylLib:
*/
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/Nucleotide/K80.h>
#include <Bpp/Phyl/Model/Nucleotide/GTR.h>
#include <Bpp/Phyl/Distance/DistanceMethod.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/BipartitionList.h>



#include <Alignment.hpp>

#include "utils.hpp"



#include <glog/logging.h>
#include <gflags/gflags.h>



static bool ValidateString(const char* flagname, const std::string& path) {
    std::ifstream file(path);
    return file.good();
}

DEFINE_string(input_sequences, "path/to/fasta.fa", "File path containing the sequences [FASTA]");
DEFINE_validator(input_sequences, &ValidateString);

DEFINE_string(input_tree, "path/to/newick.nwk", "File path containing the tree file [NWK]");
DEFINE_validator(input_tree, &ValidateString);

//DECLARE_bool(big_menu);
DECLARE_string(input_sequences);
DECLARE_string(input_tree);




int main(int argc, char *argv[]) {

    FLAGS_alsologtostderr = true;
    google::InitGoogleLogging("JATI-minimal");
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    LOG(INFO) << "JATI-minimal v.0.0.1 build: cpp0001 " << std::endl;

    //------------------------------------------------------------------------------------------------------------------
    // LOAD MSA FROM FILE
    // Parse fasta file containing aligned DNA sequences


    auto alignment = new Alignment_DNA;

    LOG(INFO) << "Input file: " << FLAGS_input_sequences;

    bpp::Fasta seqReader;
    bpp::SequenceContainer *sequences = seqReader.readAlignment(FLAGS_input_sequences, &bpp::AlphabetTools::DNA_ALPHABET);
    std::vector<std::string> seqNames = sequences->getSequencesNames();

    for(int i=0; i<seqNames.size(); i++){
        std::string stringSeq = sequences->getSequence(seqNames.at(i)).toString();
        alignment->addSequence(seqNames.at(i),stringSeq);
    }
    alignment->getAlignmentSize();
    alignment->align_num_characters.resize((unsigned long) alignment->align_length);
    alignment->align_alphabetsize += 1; // DNA +1 per PIP



    bpp::SiteContainer *sites = new bpp::VectorSiteContainer(*sequences);

    size_t num_leaves = sequences->getNumberOfSequences();

    LOG(INFO) << "[Sequences in MSA] Leaves: " << num_leaves;
    std::string testSeq = sequences->getSequence(seqNames.at(0)).toString();
    bpp::Site testSite = nullptr;

    std::vector<int> countCharNumber(sites->getNumberOfSites(), 0);

    for (int i=0; i<sites->getNumberOfSites(); i++){
       testSite = sites->getSite((size_t) i);
        for (int j=0; j<testSite.size(); j++){
            if(testSite.getValue(j)>=0){
                countCharNumber.at(i) += 1;
            }

        }
    }

    alignment->align_num_characters = countCharNumber;

    delete sequences;

    //------------------------------------------------------------------------------------------------------------------
    // INIT ROOTED TREE

    auto * newickReader = new bpp::Newick(false); //No comment allowed!
    bpp::Tree *tree = nullptr;
    try {
        tree = newickReader->read(FLAGS_input_tree); // Tree in file MyTestTree.dnd
        cout << "Tree has " << tree->getNumberOfLeaves() << " leaves." << endl;
    } catch (bpp::Exception e) {
        cout << "Error when reading tree." << endl;
    }
    auto ttTree = bpp::TreeTemplate<bpp::Node>(*tree);

    auto utree = new Utree();
    UtreeBppUtils::convertUtree(&ttTree, utree);
    LOG(INFO) << "[Initial Utree Topology] " << utree->printTreeNewick(true);

    delete tree;
    delete newickReader;


    utree->prepareSetADesCountOnNodes((int) alignment->getAlignmentSize(), alignment->align_alphabetsize);
    UtreeUtils::associateNode2Alignment(alignment, utree);

    // Add the root
    utree->addVirtualRootNode();
    utree->rootnode->initialiseLikelihoodComponents((int) alignment->getAlignmentSize(), alignment->align_alphabetsize);


    double k80_kappa = 0.5;

    // Set the substitution model
    bpp::SubstitutionModel *model_jc69 = new bpp::JCnuc(&bpp::AlphabetTools::DNA_ALPHABET);
    bpp::SubstitutionModel *model_k80 = new bpp::K80(&bpp::AlphabetTools::DNA_ALPHABET, k80_kappa);
    bpp::SubstitutionModel *model_GTR = new bpp::GTR(&bpp::AlphabetTools::DNA_ALPHABET);

    // Fill Q matrix as for JC69
    Eigen::MatrixXd Q = MatrixBppUtils::Matrix2Eigen(model_jc69->getGenerator());
    std::cerr << Q << std::endl;

    // set Pi, steady state frequencies
    auto freqs = model_jc69->getFrequencies();
    Eigen::VectorXd pi = MatrixBppUtils::Vector2Eigen(freqs);

    double mu = 0.1;
    double lambda = 0.2;

    long rows = Q.rows();
    long cols = Q.cols();

    (Q).conservativeResize(rows+1, cols+1);

    rows = Q.rows();
    cols = Q.cols();

    // Fill Q matrix as for JC69+PIP
    for(int r = 0; r<rows-1; r++){
        for(int c=0; c<cols-1; c++ ){

            if(r==c){
                Q(r,c) =  Q(r,c)-mu;
            }else{
                Q(r,c) =   Q(r,c);
            }

        }
        Q(r,cols-1) = mu;
    }

    std::cerr << Q << std::endl;

>>>>>>> master
    return 0;
}