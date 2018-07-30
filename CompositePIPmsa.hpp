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
 * @file pPIP.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 19 02 2018
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

#include "PIPnode.hpp"

#ifndef MINIJATI_COMPOSITEPIPMSA_HPP
#define MINIJATI_COMPOSITEPIPMSA_HPP

namespace bpp {

    //*******************************************************************************************

    typedef std::string MSAcolumn_t; // MSA column type
    typedef std::vector<MSAcolumn_t> MSA_t; // MSA as vector of columns

    //*******************************************************************************************
    class PIPmsa {

    private:

    public:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

        double score_;

        std::vector<std::string>  seqNames_; // vector of strings (sequence names)

        MSA_t msa_;

        std::vector<int>  traceback_path_;

        std::vector<int> traceback_map_;
//        std::vector<int> traceback_mapL_;
//        std::vector<int> traceback_mapR_;

        std::vector<int> map_compressed_seqs_;

        std::vector<int> rev_map_compressed_seqs_;

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        PIPmsa() { score_ = -std::numeric_limits<double>::infinity(); }

        void _setSeqNameLeaf(std::string &seqName);

        void _setMSAleaf(const bpp::Sequence *sequence);

        void _setSeqNameNode(std::vector<std::string> &seqNamesL,
                             std::vector<std::string> &seqNamesR);

        int getNumSites();

        void _setTracebackPathleaves();

        std::vector<string> getSequenceNames(){ return seqNames_; };

    };
    //*******************************************************************************************
    class iPIPmsa {

    public:

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        void _setMSAleaves(bpp::Sequence &sequence, std::string &seqName);

        void add(PIPmsa *msa) {};

    };
    //*******************************************************************************************
    class PIPmsaSingle : public iPIPmsa {

    private:

        //***************************************************************************************
        // PRIVATE FIELDS
        //***************************************************************************************

    public:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

        PIPmsa *pipmsa;

//        int subMSAidxL_;
//        int subMSAidxR_;

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        PIPmsaSingle(){
//            subMSAidxL_ = 0;
//            subMSAidxR_ = 0;
        }

        void _compressMSA(const bpp::Alphabet *alphabet);

        void _build_MSA(MSA_t &msaL,MSA_t &msaR);

        int getMSAlength();

        int getCompressedMSAlength();

        PIPmsa *_getMSA(){ return pipmsa; }

        void add(PIPmsa *x);

        void _compress_lk_components( std::vector<double> &lk_down_not_compressed,
                                                    std::vector<std::vector<bpp::ColMatrix<double> > > &fv_data_not_compressed);

        ~PIPmsaSingle() {
            delete pipmsa;
        }
    };
    //*******************************************************************************************
    class PIPmsaComp : public iPIPmsa {

    private:

        //***************************************************************************************
        // PRIVATE FIELDS
        //***************************************************************************************

    public:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

        std::vector<PIPmsa *> pipmsa;

        std::vector<int> subMSAidxL_;
        std::vector<int> subMSAidxR_;

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        PIPmsaComp(){}

        void _compressMSA(const bpp::Alphabet *alphabet,int idx_sb);

        void _build_MSA(MSA_t &msaL,MSA_t &msaR,int idx_sb);

        int getMSAlength(int idx);

        int getCompressedMSAlength(int idx);

        PIPmsa *_getMSA(int idx){ return pipmsa.at(idx); }

        void add(PIPmsa *x);

        void _compress_lk_components(std::vector<double> &lk_down_not_compressed,
                                     std::vector<std::vector<bpp::ColMatrix<double> > > &fv_data_not_compressed,
                                     int idx_sb);

        ~PIPmsaComp() {
            for (std::vector<PIPmsa *>::const_iterator iter = pipmsa.begin(); iter != pipmsa.end(); ++iter) {
                delete *iter;
            }
        }
    };
    //*******************************************************************************************

}
//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************
namespace compositePIPmsaUtils {

    std::vector<std::string> siteContainer2sequenceVector(std::vector<bpp::MSAcolumn_t> &MSA);

    std::vector<int> reverse_map(std::vector<int> &m);

    bpp::SiteContainer *pPIPmsa2Sites(const bpp::Alphabet *alphabet,
                                      std::vector<std::string> &seqNames,
                                      std::vector<std::string> &MSA);

}
//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************

#endif //MINIJATI_COMPOSITEPIPMSA_HPP
