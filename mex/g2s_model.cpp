#include "g2s_model.hpp"

#include "class_handle.hpp"
#include "xmlmodelfile.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "modelsettings.hpp"


G2SModel::G2SModel(std::string path) {
    //NRLib::LogKit::SetFileLog("g2s_mex.log", NRLib::LogKit::L_Low);
    //NRLib::LogKit::StartBuffering();

    XmlModelFile modelFile(path);
    if (modelFile.getParsingFailed()) {
        //NRLib::LogKit::EndLog();
        mexErrMsgIdAndTxt("G2SModel:Constructor", "Parsing of '%s' failed!", path.data());
    } else {
        //NRLib::LogKit::EndLog();
    }

    this->_model_settings = modelFile.getModelSettings();

    this->_bool_functions["useCornerPointInterpolationInDepth"] = &ModelSettings::GetUseCornerpointInterpol;
    this->_bool_functions["isRicker"] = &ModelSettings::GetRicker;

    this->_double_functions["getZeroThicknessLimit"] = &ModelSettings::GetZeroThicknessLimit;
    this->_double_functions["getTheta0"] = &ModelSettings::GetTheta0;
    this->_double_functions["getDTheta"] = &ModelSettings::GetDTheta;
    this->_double_functions["getThetaMax"] = &ModelSettings::GetThetaMax;
    this->_double_functions["getPeakFrequency"] = &ModelSettings::GetPeakFrequency;
    this->_double_functions["getWaveletScale"] = &ModelSettings::GetWaveletScale;
    this->_double_functions["getStandardDeviation"] = &ModelSettings::GetStandardDeviation;

    this->_double_vector_functions["getVsConstants"] = &ModelSettings::GetConstVs;
    this->_double_vector_functions["getVpConstants"] = &ModelSettings::GetConstVp;
    this->_double_vector_functions["getRhoConstants"] = &ModelSettings::GetConstRho;
    this->_double_vector_functions["getExtraParametersDefaultValue"] = &ModelSettings::GetExtraParameterDefaultValues;

    this->_string_functions["getEclipseFilename"] = &ModelSettings::GetEclipseFileName;
    this->_string_functions["getWaveletFileFormat"] = &ModelSettings::GetWaveletFileFormat;
    this->_string_functions["getWaveletFileName"] = &ModelSettings::GetWaveletFileName;
    this->_string_functions["getAreaFromSurface"] = &ModelSettings::GetAreaFromSurface;

    this->_string_vector_functions["getParameterName"] = &ModelSettings::GetParameterNames;
    this->_string_vector_functions["getExtraParameters"] = &ModelSettings::GetExtraParameterNames;

    this->_unsigned_long_functions["getSeed"] = &ModelSettings::GetSeed;

}

ModelSettings *G2SModel::getModelSettings() {
    return this->_model_settings;
}

bool G2SModel::isDoubleVectorFunction(const std::string &command) {
    return this->_double_vector_functions.find(command) != this->_double_vector_functions.end();
}

bool G2SModel::isDoubleFunction(const std::string &command) {
    return this->_double_functions.find(command) != this->_double_functions.end();
}

bool G2SModel::isBoolFunction(const std::string &command) {
    return this->_bool_functions.find(command) != this->_bool_functions.end();
}

bool G2SModel::isIntegerFunction(const std::string &command) {
    return this->_integer_functions.find(command) != this->_integer_functions.end();
}

bool G2SModel::isUnsignedLongFunction(const std::string &command) {
    return this->_unsigned_long_functions.find(command) != this->_unsigned_long_functions.end();
}

bool G2SModel::isStringFunction(const std::string &command) {
    return this->_string_functions.find(command) != this->_string_functions.end();
}

bool G2SModel::isStringVectorFunction(const std::string &command) {
    return this->_string_vector_functions.find(command) != this->_string_vector_functions.end();
}

bool G2SModel::hasFunction(const std::string &command) {
    return this->isDoubleFunction(command)
            || this->isBoolFunction(command)
            || this->isDoubleVectorFunction(command)
            || this->isIntegerFunction(command)
            || this->isUnsignedLongFunction(command)
            || this->isStringFunction(command)
            || this->isStringVectorFunction(command);
}

void G2SModel::initiateCommand(const std::string &command, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    std::string function_id = std::string("G2SModel:") + command;
    if (nlhs < 0 || nrhs < 2) {
        mexErrMsgIdAndTxt(function_id.data(), "Unexpected arguments. Expected at least 1 output and a minimum of 2 input arguments.");
    }


    if (this->isDoubleFunction(command)) {
        DoubleFunction function = this->_double_functions[command];
        double value = (this->getModelSettings()->*function)();
        plhs[0] = mxCreateDoubleScalar(value);
    } else if (this->isBoolFunction(command)) {
        BoolFunction function = this->_bool_functions[command];
        bool value = (this->getModelSettings()->*function)();
        plhs[0] = mxCreateLogicalScalar((mxLogical) value);
    } else if (this->isDoubleVectorFunction(command)) {
        DoubleVectorFunction function = this->_double_vector_functions[command];
        std::vector<double> values = (this->getModelSettings()->*function)();
        int length = values.size();

        if (nrhs == 3) {
            int index = (int) mxGetScalar(prhs[2]);

            if (index < 0 || index >= length) {
                mexErrMsgIdAndTxt(function_id.data(), "Expected an index in the range [0, %d] got: %d.", (length - 1), index);
            }

            plhs[0] = mxCreateDoubleScalar(values[index]);

        } else {
            plhs[0] = mxCreateDoubleMatrix((mwSize) length, (mwSize) 1, mxREAL);
            double *mex_values = mxGetPr(plhs[0]);

            for (int i = 0; i < length; i++) {
                mex_values[i] = values[i];
            }
        }
    } else if (this->isStringFunction(command)) {
        StringFunction function = this->_string_functions[command];
        std::string value = (this->getModelSettings()->*function)();
        plhs[0] = mxCreateString(value.data());
    } else if (this->isStringVectorFunction(command)) {
        StringVectorFunction function = this->_string_vector_functions[command];
        std::vector<std::string> values = (this->getModelSettings()->*function)();
        int length = values.size();

        if (nrhs == 3) {
            int index = (int) mxGetScalar(prhs[2]);

            if (index < 0 || index >= length) {
                mexErrMsgIdAndTxt(function_id.data(), "Expected an index in the range [0, %d] got: %d.", (length - 1), index);
            }

            plhs[0] = mxCreateString(values[index].data());
        } else {
            const mwSize dims[1] = {length};
            plhs[0] = mxCreateCellArray((mwSize) 1, dims);

            for (int i = 0; i < length; i++) {
                mxSetCell(plhs[0], i, mxCreateString(values[i].data()));
            }
        }
    } else if (this->isIntegerFunction(command)) {
        IntegerFunction function = this->_integer_functions[command];
        int value = (this->getModelSettings()->*function)();
        const mwSize dims[1] = {1};
        plhs[0] = mxCreateNumericArray(1, dims, mxINT32_CLASS, mxREAL);
        int *mex_values = (int *) mxGetData(plhs[0]);
        mex_values[0] = value;

    } else if (this->isUnsignedLongFunction(command)) {
        UnsignedLongFunction function = this->_unsigned_long_functions[command];
        unsigned long value = (this->getModelSettings()->*function)();
        const mwSize dims[1] = {1};
        plhs[0] = mxCreateNumericArray(1, dims, mxUINT64_CLASS, mxREAL);
        unsigned long *mex_values = (unsigned long *) mxGetData(plhs[0]);
        mex_values[0] = value;
    }

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Get the command string
    std::string command = mxArrayToString(prhs[0]);

    if (command == "new") {
        // Check parameters
        if (nlhs != 1) {
            mexErrMsgTxt("New: One output expected.");
        }
        if (nrhs < 2) {
            mexErrMsgTxt("Second input should be a filename string.");
        }
        // Return a handle to a new C++ instance
        char *path = mxArrayToString(prhs[1]);
        plhs[0] = convertPtr2Mat<G2SModel>(new G2SModel(path));
        return;
    }

    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2) {
        mexErrMsgTxt("Second input should be a class instance handle.");
    }

    // Delete
    if (command == "delete") {
        // Destroy the C++ object
        destroyObject<G2SModel>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2) {
            mexErrMsgIdAndTxt("G2SModel:delete", "Unexpected arguments ignored.");
        }
        return;
    }

    // Get the class instance pointer from the second input
    G2SModel *model_instance = convertMat2Ptr<G2SModel>(prhs[1]);

    if (model_instance->hasFunction(command)) {
        model_instance->initiateCommand(command, nlhs, plhs, nrhs, prhs);
    } else {
        // Got here, so command not recognized
        mexErrMsgIdAndTxt("G2SModel:mexFunction", "Command '%s' not recognized!", command.data());
    }
}
