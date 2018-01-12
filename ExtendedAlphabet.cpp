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
 * @file ExtendedAlphabet.cpp
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

#include <Bpp/Text/TextTools.h>
#include <Bpp/Utils/MapTools.h>

using namespace bpp;

// From STL:

using namespace std;

#include "ExtendedAlphabet.hpp"


DNA_EXTENDED::DNA_EXTENDED(bool exclamationMarkCountsAsGap) {
    // Alphabet content definition
    // all unresolved bases use n°14
    registerState(new NucleicAlphabetState(0, "A", 1, "Adenine"));
    registerState(new NucleicAlphabetState(1, "C", 2, "Cytosine"));
    registerState(new NucleicAlphabetState(2, "G", 4, "Guanine"));
    registerState(new NucleicAlphabetState(3, "T", 8, "Thymine"));
    registerState(new NucleicAlphabetState(4, "-", 16, "Gap"));
    registerState(new NucleicAlphabetState(5, "M", 3, "Adenine or Cytosine"));
    registerState(new NucleicAlphabetState(6, "R", 5, "Purine (Adenine or Guanine)"));
    registerState(new NucleicAlphabetState(7, "W", 9, "Adenine or Thymine"));
    registerState(new NucleicAlphabetState(8, "S", 6, "Cytosine or Guanine"));
    registerState(new NucleicAlphabetState(9, "Y", 10, "Pyrimidine (Cytosine or Thymine)"));
    registerState(new NucleicAlphabetState(10, "K", 12, "Guanine or Thymine"));
    registerState(new NucleicAlphabetState(11, "V", 7, "Adenine or Cytosine or Guanine"));
    registerState(new NucleicAlphabetState(12, "H", 11, "Adenine or Cytosine or Thymine"));
    registerState(new NucleicAlphabetState(13, "D", 13, "Adenine or Guanine or Thymine"));
    registerState(new NucleicAlphabetState(14, "B", 14, "Cytosine or Guanine or Thymine"));
    registerState(new NucleicAlphabetState(15, "N", 15, "Unresolved base"));
    registerState(new NucleicAlphabetState(15, "X", 15, "Unresolved base"));
    registerState(new NucleicAlphabetState(15, "O", 15, "Unresolved base"));
    registerState(new NucleicAlphabetState(15, "0", 15, "Unresolved base"));
    registerState(new NucleicAlphabetState(15, "?", 15, "Unresolved base"));

    if (exclamationMarkCountsAsGap)
        registerState(new NucleicAlphabetState(-1, "!", 0, "Frameshift"));
    else
        registerState(new NucleicAlphabetState(15, "!", 15, "Unresolved base"));

}

std::vector<int> DNA_EXTENDED::getAlias(int state) const throw(BadIntException) {
    if (!isIntInAlphabet(state))
        throw BadIntException(state, "DNA_EXTENDED::getAlias(int): Specified base unknown.");
    vector<int> v;
    const NucleicAlphabetState &st = getState(state);
    if (state == -1)
        v.push_back(-1);
    if (st.getBinaryCode() & 1)
        v.push_back(0);
    if (st.getBinaryCode() & 2)
        v.push_back(1);
    if (st.getBinaryCode() & 4)
        v.push_back(2);
    if (st.getBinaryCode() & 8)
        v.push_back(3);
    if (st.getBinaryCode() & 16)
        v.push_back(4);


    return v;
}

std::vector<std::string> DNA_EXTENDED::getAlias(const std::string &state) const throw(BadCharException) {
    string locstate = TextTools::toUpper(state);
    if (!isCharInAlphabet(locstate)) throw BadCharException(locstate, "DNA_EXTENDED::getAlias(int): Specified base unknown.");
    vector<int> vi = this->getAlias(this->charToInt(state));
    vector<string> v;
    for (unsigned int i = 0; i < vi.size(); i++)
        v.push_back(this->intToChar(vi[i]));
    return v;
}


int DNA_EXTENDED::getGeneric(const std::vector<int> &states) const throw(BadIntException) {
    int v = 0;
    for (size_t i = 0; i < states.size(); ++i) {
        if (!isIntInAlphabet(states[i])) throw BadIntException(states[i], "DNA_EXTENDED::getGeneric(const vector<int>& states): Specified base unknown.");
        v |= getState(states[i]).getBinaryCode();
    }
    return getStateByBinCode(v).getNum();
}



std::string DNA_EXTENDED::getGeneric(const std::vector<std::string> &states) const throw(BadCharException) {
    vector<int> vi;
    for (unsigned int i = 0; i < states.size(); ++i) {
        if (!isCharInAlphabet(states[i])) throw BadCharException(states[i], "DNA_EXTENDED::getGeneric(const vector<string>& states): Specified base unknown.");
        vi.push_back(this->charToInt(states[i]));
    }
    return intToChar(getGeneric(vi));
}

