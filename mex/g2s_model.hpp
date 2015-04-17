#ifndef __G2S_MODEL_HPP__
#define __G2S_MODEL_HPP__

#include <string>
#include <vector>
#include <map>
#include "mex.h"
#include "modelsettings.hpp"

typedef bool (ModelSettings::*BoolFunction)(void);

typedef double (ModelSettings::*DoubleFunction)(void);

typedef std::vector<double> (ModelSettings::*DoubleVectorFunction)(void);

typedef int (ModelSettings::*IntegerFunction)(void);

typedef unsigned long (ModelSettings::*UnsignedLongFunction)(void);

typedef std::string (ModelSettings::*StringFunction)(void);

typedef std::vector<std::string> (ModelSettings::*StringVectorFunction)(void);


class G2SModel {

  public:
    G2SModel(std::string path);

    ~G2SModel() {
        delete this->_model_settings;
    }

    ModelSettings *getModelSettings();

    bool hasFunction(const std::string &command);

    bool isBoolFunction(const std::string &command);

    bool isDoubleFunction(const std::string &command);

    bool isDoubleVectorFunction(const std::string &command);

    bool isIntegerFunction(const std::string &command);

    bool isUnsignedLongFunction(const std::string &command);

    bool isStringFunction(const std::string &command);

    bool isStringVectorFunction(const std::string &command);

    void initiateCommand(const std::string &command, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

  private:
    ModelSettings *_model_settings;
    std::map<std::string, DoubleFunction> _double_functions;
    std::map<std::string, BoolFunction> _bool_functions;
    std::map<std::string, DoubleVectorFunction> _double_vector_functions;
    std::map<std::string, IntegerFunction> _integer_functions;
    std::map<std::string, UnsignedLongFunction> _unsigned_long_functions;
    std::map<std::string, StringFunction> _string_functions;
    std::map<std::string, StringVectorFunction> _string_vector_functions;

};


#endif
