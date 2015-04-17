#include "surface_container.hpp"

#include "class_handle.hpp"
#include "xmlmodelfile.hpp"
#include "modelsettings.hpp"
#include <seismic_geometry.hpp>
#include <seismic_parameters.hpp>
#include <seismic_regridding.hpp>
#include <seismic_forward.hpp>


SurfaceContainer::SurfaceContainer(std::string path) {
    _vpgrid = NULL;
    _vsgrid = NULL;
    _rhogrid = NULL;
    _zgrid = NULL;
    _twtgrid = NULL;
    _toptime = NULL;
    _bottime = NULL;
    _topeclipse = NULL;
    _boteclipse = NULL;

    XmlModelFile modelFile(path);
    if (modelFile.getParsingFailed()) {
        mexErrMsgIdAndTxt("SurfaceContainer:Constructor", "Parsing of '%s' failed!", path.data());
    } else {
        mexPrintf("Regridding started. This may take some time...!\n");
        _model_settings = modelFile.getModelSettings();
        _seismic_parameters = new SeismicParameters(_model_settings);
        SeismicRegridding::seismicRegridding(* _seismic_parameters);
        //SeismicForward::seismicForward(* _seismic_parameters);
        mexPrintf("Regridding successful!\n");
    }

}

mxArray * SurfaceContainer::createArrayFromGrid(NRLib::StormContGrid & grid) {
    size_t nx = grid.GetNI();
    size_t ny = grid.GetNJ();
    size_t nz = grid.GetNK();

    const mwSize dims[3] = {nx, ny, nz};
    mxArray * array = mxCreateNumericArray((mwSize) 3, dims, mxSINGLE_CLASS, mxREAL);
    mexMakeArrayPersistent(array);
    return array;
}

mxArray * SurfaceContainer::createArrayFromGrid(NRLib::RegularSurface<double> & grid) {
    size_t nx = grid.GetNI();
    size_t ny = grid.GetNJ();

    const mwSize dims[2] = {nx, ny};
    mxArray * array  = mxCreateNumericArray((mwSize) 2, dims, mxDOUBLE_CLASS, mxREAL);
    
    mexMakeArrayPersistent(array);
    return array;
}


void SurfaceContainer::copyValuesFromGridToArray(NRLib::StormContGrid & grid, mxArray * array) {
    size_t nx = grid.GetNI();
    size_t ny = grid.GetNJ();
    size_t nz = grid.GetNK();

    float * data = (float *) mxGetData(array);
    mwIndex * subs = (mwIndex *) mxCalloc(3, sizeof(mwIndex));

    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t k = 0; k < nz; k++) {
                subs[0] = i;
                subs[1] = j;
                subs[2] = k;
                mwIndex index = mxCalcSingleSubscript(array, (mwSize) 3, subs);
                data[index] = grid(i, j, k);
            }
        }
    }
    mxFree(subs);
}

void SurfaceContainer::copyValuesFromGridToArray(NRLib::RegularSurface<double> & grid, mxArray * array) {
    size_t nx = grid.GetNI();
    size_t ny = grid.GetNJ();

    double * data = (double *) mxGetData(array);
    mwIndex * subs = (mwIndex *) mxCalloc(2, sizeof(mwIndex));

    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            subs[0] = i;
            subs[1] = j;
            mwIndex index = mxCalcSingleSubscript(array, (mwSize) 2, subs);
            data[index] = grid(i, j);
        }
    }
    mxFree(subs);
}



mxArray * SurfaceContainer::zgrid() {
    if (_zgrid == NULL) {
        _zgrid = createArrayFromGrid(_seismic_parameters->zGrid());
        copyValuesFromGridToArray(_seismic_parameters->zGrid(), _zgrid);
    }
    return _zgrid;
}

mxArray * SurfaceContainer::vpgrid() {
    if (_vpgrid == NULL) {
        _vpgrid = createArrayFromGrid(_seismic_parameters->vpGrid());
        copyValuesFromGridToArray(_seismic_parameters->vpGrid(), _vpgrid);
    }
    return _vpgrid;
}

mxArray * SurfaceContainer::vsgrid() {
    if (_vsgrid == NULL) {
        _vsgrid = createArrayFromGrid(_seismic_parameters->vsGrid());
        copyValuesFromGridToArray(_seismic_parameters->vsGrid(), _vsgrid);
    }
    return _vsgrid;
}

mxArray * SurfaceContainer::rhogrid() {
    if (_rhogrid == NULL) {
        _rhogrid = createArrayFromGrid(_seismic_parameters->rhoGrid());
        copyValuesFromGridToArray(_seismic_parameters->rhoGrid(), _rhogrid);
    }
    return _rhogrid;
}

mxArray * SurfaceContainer::twtgrid() {
    if (_twtgrid == NULL) {
        _twtgrid = createArrayFromGrid(_seismic_parameters->twtGrid());
        copyValuesFromGridToArray(_seismic_parameters->twtGrid(), _twtgrid);
    }
    return _twtgrid;
}

mxArray * SurfaceContainer::toptime() {
    if(_toptime == NULL) {
        _toptime = createArrayFromGrid(_seismic_parameters->topTime());
        copyValuesFromGridToArray(_seismic_parameters->topTime(), _toptime);
    }
    return _toptime;
}

mxArray * SurfaceContainer::bottime() {
    if(_bottime == NULL) {
        _bottime = createArrayFromGrid(_seismic_parameters->bottomTime());
        copyValuesFromGridToArray(_seismic_parameters->bottomTime(), _bottime);
    }
    return _bottime;
}

mxArray * SurfaceContainer::topeclipse() {
    if(_topeclipse == NULL) {
        _topeclipse = createArrayFromGrid(_seismic_parameters->topEclipse());
        copyValuesFromGridToArray(_seismic_parameters->topEclipse(), _topeclipse);
    }
    return _topeclipse;
}

mxArray * SurfaceContainer::boteclipse() {
    if(_boteclipse == NULL) {
        _boteclipse = createArrayFromGrid(_seismic_parameters->bottomEclipse());
        copyValuesFromGridToArray(_seismic_parameters->bottomEclipse(), _boteclipse);
    }
    return _boteclipse;
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
    // Get the command string
    std::string command = mxArrayToString(prhs[0]);

    if (command == "new") {
        // Check parameters
        if (nlhs != 1) {
            mexErrMsgIdAndTxt("SurfaceContainer:new", "One output expected.");
        }
        if (nrhs < 2) {
            mexErrMsgIdAndTxt("SurfaceContainer:new", "Second input should be a filename string.");
        }
        // Return a handle to a new C++ instance
        char * path = mxArrayToString(prhs[1]);
        plhs[0] = convertPtr2Mat<SurfaceContainer>(new SurfaceContainer(path));
        return;
    }

    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2) {
        mexErrMsgIdAndTxt("SurfaceContainer:*", "Second input should be a class instance handle.");
    }

    // Delete
    if (command == "delete") {
        // Destroy the C++ object
        destroyObject<SurfaceContainer>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2) {
            mexErrMsgIdAndTxt("SurfaceContainer:delete", "Unexpected arguments ignored.");
        }
        return;
    }

    // Get the class instance pointer from the second input
    SurfaceContainer * model_instance = convertMat2Ptr<SurfaceContainer>(prhs[1]);


    if (command == "zgrid") {
        plhs[0] = model_instance->zgrid();
        return;

    } else if (command == "twtgrid") {
        plhs[0] = model_instance->twtgrid();
        return;

    } else if (command == "vpgrid") {
        plhs[0] = model_instance->vpgrid();
        return;

    } else if (command == "vsgrid") {
        plhs[0] = model_instance->vsgrid();
        return;

    } else if (command == "rhogrid") {
        plhs[0] = model_instance->rhogrid();
        return;

    } else if (command == "toptime") {
        plhs[0] = model_instance->toptime();
        return;

    } else if (command == "bottime") {
        plhs[0] = model_instance->bottime();
        return;

    } else if (command == "topeclipse") {
        plhs[0] = model_instance->topeclipse();
        return;

    } else if (command == "boteclipse") {
        plhs[0] = model_instance->boteclipse();
        return;

    } else {
        // Got here, so command not recognized
        mexErrMsgIdAndTxt("SurfaceContainer:mexFunction", "Command '%s' not recognized!", command.data());
    }

}

