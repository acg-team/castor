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
 * @file cli_parser.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 04 01 2018
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
#ifndef MINIJATI_CLI_PARSER_HPP
#define MINIJATI_CLI_PARSER_HPP

#include <gflags/gflags.h>

static bool ValidateString(const char* flagname, const std::string& path) {
    std::ifstream file(path);
    return file.good();
}

static bool ValidateOptimTopology(const char* flagname, const std::string& option) {
    bool status = false;
    if(option == "full-search" || option == "smart-search" || option == "nni-search" || option == "spr-search" || option == "no-search"){
        status = true;
    }

    return status;
}

// 1. Input

// 1.1 MSA
DEFINE_string(input_sequences, "", "File path containing the sequences [FASTA]");
DEFINE_validator(input_sequences, &ValidateString);

// 1.2 Tree
DEFINE_string(input_tree, "", "File path containing the tree file [NWK]");
//DEFINE_validator(input_tree, &ValidateString);

// 2. Optimisations

// Topology
DEFINE_string(optim_topology, "full-search,smart-search,nni-search,spr-search,no-search", "Topology optimisation under a predefined scheme");
DEFINE_validator(optim_topology, &ValidateOptimTopology);


// 3. Models
DEFINE_string(model_substitution, "GTR", "Substitution model to apply on the input data");
DEFINE_bool(model_indels, false, "Extend the substitution model to include Insertion & Deletion events (PIP)");
DEFINE_bool(model_setfreqsfromdata, false, "Set substitution model frequencies from data");


// Compute likelihood
DEFINE_bool(lkmove_bothways, false, "Compute likelihood of the model for each topology rearrangment (apply and revert)");

// Compute alignment
DEFINE_bool(alignment, false, "Compute alignment of the given fasta file");



// 3.1 PIP params
DEFINE_double(lambda_PIP, 1.0, "insertion rate (lambda)");
DEFINE_double(mu_PIP, 1.0, "deletion rate (mu)");

// 4. Output

// 4.1 MSA
DEFINE_string(output_msa, "path/to/fasta.fa", "File path for output MSA [FASTA]");

// 4.2 LK
DEFINE_string(output_lk, "path/to/lk", "File path for output lk");


// Declarations
DECLARE_string(model_substitution);
DECLARE_bool(model_indels);
DECLARE_bool(model_setfreqsfromdata);
DECLARE_bool(lkmove_bothways);
DECLARE_bool(alignment);
DECLARE_string(input_sequences);
DECLARE_string(input_tree);

DECLARE_string(output_msa);
DECLARE_string(output_lk);
DECLARE_double(lambda_PIP);
DECLARE_double(mu_PIP);

#endif //MINIJATI_CLI_PARSER_HPP
