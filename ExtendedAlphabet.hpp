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
 * @file ExtendedAlphabet.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 11 01 2018
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
#ifndef MINIJATI_EXTENDEDALPHABET_HPP
#define MINIJATI_EXTENDEDALPHABET_HPP


#include <Bpp/Seq/Alphabet/LetterAlphabet.h>
#include <Bpp/Seq/Alphabet/NucleicAlphabetState.h>
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>


namespace bpp {


    /**
     * @brief This alphabet is used to deal with DNA_EXTENDED sequences.
     *
     * It supports all 4 nucleotides (A, T, G and C) with their standard denomination.
     * Gaps are coded by '-', unresolved characters are coded by 'X, N, O, 0 or ?'.
     * Extensive support for generic characters (e.g. 'P', 'Y', etc.) is provided.
     */
    class DNA_EXTENDED : public NucleicAlphabet {
    public:
        /**
         * @param exclamationMarkCountsAsGap If yes, '!' characters are replaced by gaps.
         * Otherwise, they are counted as unknown characters.
         */
        DNA_EXTENDED(bool exclamationMarkCountsAsGap = false);

        DNA_EXTENDED(const DNA_EXTENDED &bia) : NucleicAlphabet(bia) {}

        DNA_EXTENDED &operator=(const DNA_EXTENDED &bia) {
            NucleicAlphabet::operator=(bia);
            return *this;
        }

        DNA_EXTENDED *clone() const {
            return new DNA_EXTENDED(*this);
        }

        virtual ~DNA_EXTENDED() {}

    public:
        std::vector<int> getAlias(int state) const throw(BadIntException);

        std::vector<std::string> getAlias(const std::string &state) const throw(BadCharException);

        int getGeneric(const std::vector<int> &states) const throw(BadIntException);

        std::string getGeneric(const std::vector<std::string> &states) const throw(BadCharException);

        std::string getAlphabetType() const { return "DNA_EXTENDED"; }
        // return 4 : A, C, G, T (or U)
        unsigned int getSize() const { return 5; }

        // return 15 : gap isn't included, generic unresolved bases (N, X, ?, O, 0) count for one
        unsigned int getNumberOfTypes() const { return 16; }

        int getUnknownCharacterCode() const { return 15; }

        bool isUnresolved(int state) const { return state > 3; }

        bool isUnresolved(const std::string& state) const { return charToInt(state) > 3; }

        int getGapCharacterCode() const { return 4; }
    };


}

#endif //MINIJATI_EXTENDEDALPHABET_HPP
