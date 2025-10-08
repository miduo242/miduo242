/*
 *   This file is part of the OpenPhase (R) software library.
 *
 *   Copyright (c) 2009-2022 Ruhr-Universitaet Bochum,
 *                 Universitaetsstrasse 150, D-44801 Bochum, Germany
 *             AND 2018-2022 OpenPhase Solutions GmbH,
 *                 Universitaetsstrasse 136, D-44799 Bochum, Germany.
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   File created :   2014
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev;
 *                         Philipp Engels; Johannes Goerler
 *
 */

#include "Base/UserInterface.h"
#include "Info.h"
namespace openphase
{
using namespace std;

// ======================= Auxiliary functions ===============================//
string removeEntireWhiteSpaces(const string InString)
{
    string StringTemp = InString;
    StringTemp.erase(std::remove_if( StringTemp.begin(), StringTemp.end(),
     [](char c){ return (c =='\r' || c =='\t' || c == ' ' || c == '\n');}), StringTemp.end() );
     return StringTemp;
}

string removeLeadingTrailingWhiteSpaces(const string InString)
{
    std::string whitespaces (" \t");
    size_t first = InString.find_first_not_of(whitespaces);
    size_t last = InString.find_last_not_of(whitespaces);
    if (first == last) return "";
    return InString.substr(first, (last-first+1));
}

string replacePunctuationMarksWithComma(string& Instring)
{
	Instring.erase (std::remove (Instring.begin(), Instring.end(), ' '), Instring.end());
    std::replace(Instring.begin(), Instring.end(), '/', ',');
    std::replace(Instring.begin(), Instring.end(), '-', ',');
    std::replace(Instring.begin(), Instring.end(), '@', ',');
    std::replace(Instring.begin(), Instring.end(), ':', ',');
    return Instring;
}

int StringToInt(const std::string& str, const std::string name)
{
    int ivar = 0;

    try
    {
        ivar = std::stoi(str);
    }
    catch (const std::invalid_argument&)
    {
        Info::WriteExit("Argument for $" + name + " is invalid", "UserInterface", "ReadParameterI()");
        //throw;
        exit(3);
    }
    return ivar;
}

double StringToDouble(const std::string& str, const std::string name)
{
    double dvar = 0.0;

    try
    {
        dvar = std::stod(str);
    }
    catch (const std::invalid_argument&)
    {
        Info::WriteExit("Argument for $" + name + " is invalid", "UserInterface", "ReadParameterD()");
        //throw;
        exit(3);
    }
    return dvar;
}
// ====================== Auxiliary functions end ============================//

string UserInterface::MakeFileName(string Directory, string NameBase, int Index, string FileExtension)
{
    stringstream converter;
    converter << Index;
    string Count = converter.str();
    int n = Count.size();
    string Number;
    for(int i = 0; i < 8 - n; i++) Number.append("0");
    Number.append(Count);

    return Directory + NameBase + Number + FileExtension;
}

int UserInterface::FindModuleLocation(stringstream& sInp, const string module)
{
    std::string copy = sInp.str();
    std::stringstream Inp;
    Inp << copy;
    // Returns the location of the first line in the module specified in the parameters
    Inp.seekg(0, ios::beg);
    while (Inp.good())
    {
        string tmp;
        string ReadKey;

        getline(Inp, tmp, '@');
        if (Inp.eof()) break;
        Inp >> ReadKey;

        if (!ReadKey.compare(module))
        {
            streampos temp = Inp.tellg();
            //temp += 1;
            return temp;
        }
    }
    Inp.seekg(0, ios::beg);
    return 0;
}

int UserInterface::FindParameter(stringstream& sInp,
    const int location, const string Key)
{
    std::string copy = sInp.str();
    std::stringstream Inp;
    Inp << copy;
    streampos curPos(location+1);
    Inp.seekg(curPos);

    bool found = false;
    int keyLocation = 0;
    while (!Inp.eof() and !found)
    {
        string tmp;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            keyLocation = -1;
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                keyLocation = -1;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            keyLocation = Inp.tellg();
            keyLocation -= (ReadKey.size() + 1);
            found = true;
            break;
        }
    }
    Inp.clear();
    return keyLocation;
}

int UserInterface::FindParameterLocation(stringstream& sInp,
    const int location, const string Key)
{
    std::string copy = sInp.str();
    std::stringstream Inp;
    Inp << copy;
    streampos curPos(location+1);
    Inp.seekg(curPos);

    bool found = false;
    int keyLocation = 0;
    while (!Inp.eof() and !found)
    {
        string tmp;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            keyLocation = -1;
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                keyLocation = -1;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            //getline(Inp, tmp, ':');
            keyLocation = Inp.tellg();
            //keyLocation -= (ReadKey.size() + 1);
            found = true;
            break;
        }
    }
    Inp.clear();
    return keyLocation;
}

double UserInterface::ReadParameterD(stringstream& sInp, int currentLocation, string Key,
    const bool mandatory, const double defaultval)
{
    // Note the following optional arguments:
    // mandatory: if parameter not found, program will stop (default: true = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    std::string copy = sInp.str();
    std::stringstream Inp;
    Inp << copy;
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    double ReturnValue = defaultval;

    while (!Inp.eof())
    {
        string tmp;
        string tmp2;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            getline(Inp, tmp, ':');
            getline(Inp, tmp2);
            tmp = removeLeadingTrailingWhiteSpaces(tmp);
            tmp2 = removeEntireWhiteSpaces(tmp2);

            // Transfer string tmp2 to double (and check if this is possible at all)
            ReturnValue = StringToDouble(tmp2, tmp);

            // Checks if comment is given, otherwise use "Key" as output
            if (tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::Write(Key, ReturnValue, 5);
            }
            else
            {
                Info::Write(tmp, ReturnValue, 5);
            }
            found = true;
            break;
        }

    }

    if (not found)
    {
        if (mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            // Return defaultval instead
            ReturnValue = defaultval;
        }

    }
    Inp.clear();
    return ReturnValue;
}

dVector3 UserInterface::ReadParameterV3(stringstream& sInp, int currentLocation, string Key,
    const bool mandatory, const dVector3 defaultval)
{
    // Note the following optional arguments:
    // mandatory: if parameter not found, program will stop (default: true = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    std::string copy = sInp.str();
    std::stringstream Inp;
    Inp << copy;
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    dVector3 ReturnValue = defaultval;

    while (!Inp.eof())
    {
        string tmp;
        string tmp2;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            getline(Inp, tmp, ':');
            tmp = removeLeadingTrailingWhiteSpaces(tmp);

            ReturnValue.read(Inp);

            // Checks if comment is given, otherwise use "Key" as output
            if (tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::Write(Key, ReturnValue, 5);
            }
            else
            {
                Info::Write(tmp, ReturnValue, 5);
            }
            found = true;
            break;
        }
    }

    if (not found)
    {
        if (mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            // Return defaultval instead
            ReturnValue = defaultval;
        }
    }
    Inp.clear();
    return ReturnValue;
}

dMatrix3x3 UserInterface::ReadParameterM3x3(stringstream& sInp, int currentLocation, string Key,
    const bool mandatory, const dMatrix3x3 defaultval)
{
    // Note the following optional arguments:
    // mandatory: if parameter not found, program will stop (default: true = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    std::string copy = sInp.str();
    std::stringstream Inp;
    Inp << copy;
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    dMatrix3x3 ReturnValue = defaultval;

    while (!Inp.eof())
    {
        string tmp;
        string tmp2;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            getline(Inp, tmp, ':');
            tmp = removeLeadingTrailingWhiteSpaces(tmp);

            ReturnValue.read(Inp);

            // Checks if comment is given, otherwise use "Key" as output
            if (tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::Write(Key, ReturnValue, 5);
            }
            else
            {
                Info::Write(tmp, ReturnValue, 5);
            }
            found = true;
            break;
        }
    }

    if (not found)
    {
        if (mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            // Return defaultval instead
            ReturnValue = defaultval;
        }
    }
    Inp.clear();
    return ReturnValue;
}

dMatrix6x6 UserInterface::ReadParameterT6(stringstream& sInp, int currentLocation, string Key,
    const bool mandatory, const dMatrix6x6 defaultval)
{
    // Note the following optional arguments:
    // verbose: if true, found variable is printed to terminal, standard = true
    // mandatory: if parameter not found, program will stop (standard = true = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    std::string copy = sInp.str();
    std::stringstream Inp;
    Inp << copy;
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    dMatrix6x6 ReturnValue = defaultval;

    while (!Inp.eof())
    {
        string tmp;
        string tmp2;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            getline(Inp, tmp, ':');
            tmp = removeLeadingTrailingWhiteSpaces(tmp);

            ReturnValue.read(Inp);

            // Checks if comment is given, otherwise use "Key" as output
            if (tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::Write(Key, ReturnValue, 5);
            }
            else
            {
                Info::Write(tmp, ReturnValue, 5);
            }
            found = true;
            break;
        }
    }

    if (not found)
    {
        if (mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            // Return defaultval instead
            ReturnValue = defaultval;
        }
    }
    Inp.clear();
    return ReturnValue;
}

dMatrix6x6 UserInterface::ReadParameterM6x6(stringstream& sInp, int currentLocation, string Key,
    const bool mandatory, const dMatrix6x6 defaultval)
{
    // Note the following optional arguments:
    // mandatory: if parameter not found, program will stop (default: true = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    std::string copy = sInp.str();
    std::stringstream Inp;
    Inp << copy;
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    dMatrix6x6 ReturnValue = defaultval;

    while (!Inp.eof())
    {
        string tmp;
        string tmp2;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            getline(Inp, tmp, ':');
            tmp = removeLeadingTrailingWhiteSpaces(tmp);
            ReturnValue.read(Inp);

            // Checks if comment is given, otherwise use "Key" as output
            if (tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::Write(Key, ReturnValue, 5);
            }
            else
            {
                Info::Write(tmp, ReturnValue, 5);
            }
            found = true;
            break;
        }
    }

    if (not found)
    {
        if (mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            // Return defaultval instead
            ReturnValue = defaultval;
        }
    }
    Inp.clear();
    return ReturnValue;
}

string UserInterface::ReadParameterS(stringstream& sInp, int currentLocation, string Key,
    const bool mandatory, const string defaultval)
{
    // Note the following optional arguments:
    // mandatory: if parameter not found, program will stop (default: true = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    std::string copy = sInp.str();
    std::stringstream Inp;
    Inp << copy;
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    string ReturnValue;
    while (!Inp.eof())
    {
        string tmp;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            // Yet problem with blanks in name
            getline(Inp, tmp, ':');
            Inp >> ReturnValue;
            if (tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::WriteStandard(Key, ReturnValue);
            }
            else
            {
                tmp.erase(0, tmp.find_first_not_of(" \t"));
                tmp.erase(tmp.find_last_not_of(" \t") + 1, tmp.size());
                Info::WriteStandard(tmp, ReturnValue);
            }
            found = true;
            break;
        }
    }

    if (not found)
    {
        if (mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            // Return defaultval instead
            ReturnValue = defaultval;
        }
    }

    Inp.clear();
    return ReturnValue;
}

string UserInterface::ReadParameterK(stringstream& Inp, int currentLocation, string Key,
    const bool mandatory, const string defaultval)
{
    // Note the following optional arguments:
    // mandatory: if parameter not found, program will stop (default: true = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval

    string ReturnValue = ReadParameterF(Inp, currentLocation, Key,mandatory,defaultval);

    std::transform(ReturnValue.begin(), ReturnValue.end(), ReturnValue.begin(), ::toupper);
    ReturnValue.erase(std::remove (ReturnValue.begin(), ReturnValue.end(), ' '), ReturnValue.end());

    return ReturnValue;
}

vector<string> UserInterface::ReadParameterVS(stringstream& sInp, int currentLocation, string Key,
    const bool mandatory, const vector<string> defaultval)
{
    // Note the following optional arguments:
    // mandatory: if parameter not found, program will stop (default: true = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    std::string copy = sInp.str();
    std::stringstream Inp;
    Inp << copy;
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    vector<string> ReturnValue;
    string locString;

    while (!Inp.eof())
    {
        string tmp;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            // Yet problem with blanks in name
            getline(Inp, tmp, ':');
            getline(Inp, locString);
            std::transform(locString.begin(), locString.end(), locString.begin(), ::toupper);
            locString = removeLeadingTrailingWhiteSpaces(locString);
            locString = replacePunctuationMarksWithComma(locString);
            string tmp2;
            for(unsigned int i =0; i < locString.size(); i++)
            {

                if(locString[i] != ',')
                {
                    tmp2.push_back(locString[i]);
                }
                else
                {
                    ReturnValue.push_back(tmp2);
                    tmp2.clear();
                }
            }
            ReturnValue.push_back(tmp2);

            if (tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::WriteStandard(Key, locString);
            }
            else
            {
                tmp.erase(0, tmp.find_first_not_of(" \t"));
                tmp.erase(tmp.find_last_not_of(" \t") + 1, tmp.size());
                Info::WriteStandard(tmp, locString);
            }
            found = true;
            break;
        }
    }
    if (not found)
    {
        if (mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            // Return defaultval instead
            ReturnValue = defaultval;
        }
    }
    Inp.clear();
    return ReturnValue;
}

bool UserInterface::ReadParameterB(stringstream& sInp, int currentLocation, const string Key,
    const bool mandatory, const bool defaultval)
{
    // Note the following optional arguments:
    // mandatory: if parameter not found, program will stop (default: true = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    std::string copy = sInp.str();
    std::stringstream Inp;
    Inp << copy;
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    bool ReturnValue = false;
    string Answer;
    while (!Inp.eof())
    {
        string tmp;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            getline(Inp, tmp, ':');
            getline(Inp, Answer);

            std::transform(Answer.begin(), Answer.end(), Answer.begin(), ::toupper);
            //ReadKey.erase(std::remove(Answer.begin(), Answer.end(), ' '), Answer.end());

            Answer.erase(0, Answer.find_first_not_of("YESNO") + 1);
            if (Answer.find_last_of("YESNO") + 1 != Answer.size())
            {
                Answer.erase(Answer.find_last_of("YESNO") + 1, Answer.size());
            }
            //Inp >> ReturnValue;
            if (tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::WriteStandard(Key, Answer);
            }
            else
            {
                tmp.erase(0, tmp.find_first_not_of(" \t"));
                tmp.erase(tmp.find_last_not_of(" \t") + 1, tmp.size());
                Info::WriteStandard(tmp, Answer);
            }
            found = true;
            break;
        }
    }
    string tmp;
    for (unsigned int i = 0; i < Answer.size(); ++i)
        if (Answer[i] != ' ') tmp.push_back(Answer[i]);

    if (found and !tmp.compare("YES")) ReturnValue = true;
    if (found and !tmp.compare("NO")) ReturnValue = false;
    if (not found)
    {
        if (mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            ReturnValue = defaultval;
        }
    }
    Inp.clear();
    return ReturnValue;
}

string UserInterface::ReadParameterF(stringstream& sInp, int currentLocation, const string Key,
    const bool mandatory, const string defaultval)
{
    // Note the following optional arguments:
    // mandatory: if parameter not found, program will stop (default: true = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    std::string copy = sInp.str();
    std::stringstream Inp;
    Inp << copy;
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    string ReturnValue;
    string tFileName;
    while (!Inp.eof())
    {
        string tmp;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            getline(Inp, tmp, ':');
            getline(Inp, tFileName);
            if (tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::WriteStandard(Key, tFileName);
            }
            else
            {
                tmp.erase(0, tmp.find_first_not_of(" \t"));
                tFileName.erase(0, tFileName.find_first_not_of(" \t"));
                Info::WriteStandard(tmp, tFileName);
            }
            found = true;
            break;
        }
    }
    string tmp;
    for (unsigned int i = 0; i < tFileName.size(); i++)
        if (tFileName[i] != ' ') ReturnValue.push_back(tFileName[i]);

    if (not found)
    {
        if (mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            // Return defaultval instead
            ReturnValue = defaultval;
        }
    }

    Inp.clear();
    return ReturnValue;
}

string UserInterface::ReadParameterFW(stringstream& sInp, int currentLocation, const string Key,
    const bool mandatory, const string defaultval)
{
    // Note the following optional arguments:
    // mandatory: if parameter not found, program will stop (default: true = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    std::string copy = sInp.str();
    std::stringstream Inp;
    Inp << copy;
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    string ReturnValue;
    string tFileName;
    while (!Inp.eof())
    {
        string tmp;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;
        Inp >> noskipws;
        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            getline(Inp, tmp, ':');
            getline(Inp, tFileName);
            if (tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::WriteStandard(Key, tFileName);
            }
            else
            {
                tmp.erase(0, tmp.find_first_not_of(" \t"));
                tFileName.erase(0, tFileName.find_first_not_of(" \t"));
                Info::WriteStandard(tmp, tFileName);
            }
            found = true;
            break;
        }
    }
    string tmp;
    for (unsigned int i = 0; i < tFileName.size(); i++)
        ReturnValue.push_back(tFileName[i]);

    if (not found)
    {
        if (mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            // Return defaultval instead
            ReturnValue = defaultval;
        }
    }

    Inp.clear();
    return ReturnValue;
}

int UserInterface::ReadParameterI(stringstream& sInp, int currentLocation, string Key,
    const bool mandatory, int const defaultval)
{
    // Note the following optional arguments:
    // mandatory: if parameter not found, program will stop (default: true = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    std::string copy = sInp.str();
    std::stringstream Inp;
    Inp << copy;
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    int ReturnValue = defaultval;
    while (!Inp.eof())
    {
        string tmp;
        string tmp2;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            getline(Inp, tmp, ':');
            getline(Inp, tmp2);
            tmp = removeLeadingTrailingWhiteSpaces(tmp);
            tmp2 = removeEntireWhiteSpaces(tmp2);

            // Transfer string tmp2 to int (and check if this is possible at all)
            ReturnValue = StringToInt(tmp2, tmp);

            if (tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::Write(Key, ReturnValue);
            }
            else
            {
                Info::Write(tmp, ReturnValue);
            }
            found = true;
            break;
        }
    }

    if (not found)
    {
        if (mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            // Return defaultval instead
            ReturnValue = defaultval;
        }
    }
    Inp.clear();
    return ReturnValue;
}

char UserInterface::ReadParameterC(stringstream& sInp, int currentLocation, const string Key,
                              const bool mandatory, const char defaultval)
{
    // Note the following optional arguments:
    // mandatory: if parameter not found, program will stop (default: true = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    std::string copy = sInp.str();
    std::stringstream Inp;
    Inp << copy;
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    char ReturnValue;
    while (!Inp.eof())
    {
        string tmp;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if(!ReadKey.compare(Key))
        {
            getline(Inp, tmp, ':');
            Inp >> ReturnValue;
            if(tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::Write(Key, tmp);
            }
            else
            {
                tmp.erase(0,tmp.find_first_not_of(" \t"));
                Info::Write(tmp, ReturnValue);
            }
            found = true;
            break;
        }
    }

    if(not found)
    {
        if(mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            // Return defaultval instead ( = 0 if not defined explicitly)
            ReturnValue = defaultval;
        }
    }
    Inp.clear();
    return ReturnValue;
}

}//namespace openphase
