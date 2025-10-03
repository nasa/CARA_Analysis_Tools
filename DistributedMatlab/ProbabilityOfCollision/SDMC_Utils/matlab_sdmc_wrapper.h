//====================================================================================
//
// Copyright (c) 2023-2025 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//
//====================================================================================
//
// Description:
//
//   Header definitions used by matlab_sdmc_wrapper.cpp defining the conversion
//   routines between C++ and Matlab.
//
//====================================================================================
//
// Initial version: Mar 2023; Latest update: Aug 2025
//
// ----------------- BEGIN CODE -----------------

#ifndef __MATLAB_SDMC_WRAPPER_HPP__
#define __MATLAB_SDMC_WRAPPER_HPP__

#include <string>
#include "mex.h"
#include "SdmcInterface.h"

using namespace std;

// The class that we are interfacing to
class matlab_sdmc_wrapper
{
    
private:
    string errMsg;
    
    struct SdmcInputs sdmcInputs;
    bool sdmcInputsInitialized;
    
    struct SdmcOutputs sdmcOutputs;
    bool sdmcOutputsInitialized;
    
    // Determines if the field exists in the mexStruct passed in
    //   Inputs:
    //     mexStruct - Pointer to a mex structure representing a Matlab
    //                 structure array
    //     index - Represents the indices into the mexStruct (0 indexed)
    //     fieldname - Name of the field to search in mexStruct
    //   Returns:
    //     bool - true if the field is found, false otherwise
    bool is_field(const mxArray *mexStruct, mwIndex index, string fieldname) {
        errMsg = "";
        if (NULL == mxGetField(mexStruct, index, fieldname.c_str())) {
            errMsg = "Field " + fieldname + " does not exist in structure";
            return false;
        }
        return true;
    }
    
    // Returns a pointer to an mxArray associated with the field name
    //   Inputs:
    //     mexStruct - Pointer to a mex structure representing a Matlab
    //                 structure array
    //     index - Represents the indices into the mexStruct (0 indexed)
    //     fieldname - Name of the field to search in mexStruct
    //   Returns:
    //     mxArray* - Pointer to an mxArray represented by the fieldname
    //   Errors:
    //     - when the fieldname does not exist in mexStruct
    //     - the field represented by fieldname is not a scalar value
    mxArray* get_scalar_field(const mxArray *mexStruct, mwIndex index, string fieldname) {
        if (!is_field(mexStruct, index, fieldname)) {
            mexErrMsgTxt(errMsg.c_str());
        }
        mxArray *field = mxGetField(mexStruct, index, fieldname.c_str());
        if (mxGetNumberOfElements(field) != 1) {
            errMsg = "Field " + fieldname + " must be a scalar value.";
            mexErrMsgTxt(errMsg.c_str());
        }
        return field;
    }
    
    // Returns a double precision floating point value from mexStruct that
    // is associated with fieldname
    //   Inputs:
    //     mexStruct - Pointer to a mex structure representing a Matlab
    //                 structure array
    //     index - Represents the indices into the mexStruct (0 indexed)
    //     fieldname - Name of the field to search in mexStruct
    //   Returns:
    //     double - Value associated with the fieldname
    //   Errors:
    //     - when the fieldname does not exist in mexStruct
    //     - the field represented by fieldname is not a scalar value
    double get_double_value(const mxArray *mexStruct, mwIndex index, string fieldname) {
        mxArray *field = get_scalar_field(mexStruct, index, fieldname);
#if MX_HAS_INTERLEAVED_COMPLEX
        mxDouble *valArray = mxGetDoubles(field);
#else
        mxDouble *valArray = mxGetPr(field);
#endif
        double val = valArray[0];
        return val;
    }
    
    // Returns an integer value from mexStruct that is associated with
    // fieldname
    //   Inputs:
    //     mexStruct - Pointer to a mex structure representing a Matlab
    //                 structure array
    //     index - Represents the indices into the mexStruct (0 indexed)
    //     fieldname - Name of the field to search in mexStruct
    //   Returns:
    //     int - Value associated with the fieldname
    //   Errors:
    //     - when the fieldname does not exist in mexStruct
    //     - the field represented by fieldname is not a scalar value
    int get_int_value(const mxArray *mexStruct, mwIndex index, string fieldname) {
        mxArray *field = get_scalar_field(mexStruct, index, fieldname);
#if MX_HAS_INTERLEAVED_COMPLEX
        mxDouble *valArray = mxGetDoubles(field);
#else
        mxDouble *valArray = mxGetPr(field);
#endif
        int val = (int)valArray[0];
        return val;
    }
    
    // Returns a string value from mexStruct that is associated with
    // fieldname
    //   Inputs:
    //     mexStruct - Pointer to a mex structure representing a Matlab
    //                 structure array
    //     index - Represents the indices into the mexStruct (0 indexed)
    //     fieldname - Name of the field to search in mexStruct
    //   Returns:
    //     string - Value associated with the fieldname
    //   Errors:
    //     - when the fieldname does not exist in mexStruct
    //     - the field represented by fieldname is not a scalar value
    //     - there are problems retrieving the value as a string
    string get_string_value(const mxArray *mexStruct, mwIndex index, string fieldname) {
        if (!is_field(mexStruct, index, fieldname)) {
            mexErrMsgTxt(errMsg.c_str());
        }
        mxArray *field = mxGetField(mexStruct, index, fieldname.c_str());
        int buflen = (mxGetN(field) * sizeof(mxChar)) + 1;
        char *buf = (char *)mxMalloc(buflen);
        int errStatus = mxGetString(field, buf, buflen);
        if (errStatus) {
            mxFree(buf);
            errMsg = "Errors encountered reading field '" + fieldname + "' from structure";
            mexErrMsgTxt(errMsg.c_str());
        }
        string val = string(buf);
        mxFree(buf);
        return val;
    }
    
    // Initializes the sdmcInputs structure with values from the stateInfo
    // structure
    //   Inputs:
    //     stateInfo - Pointer to a mex structure representing a Matlab
    //                 structure array which contains conjunction state
    //                 information. The stateInfo structure is populated by
    //                 the call_SDMC.m Matlab routine.
    void init_sdmc_inputs(const mxArray *stateInfo) {
        sdmcInputsInitialized = false;
        
        // Initialize the information for the primary satellite
        sdmcInputs.priSat.satno     = get_int_value(stateInfo, 0, "pri_satno");
        sdmcInputs.priSat.epoch     = get_double_value(stateInfo, 0, "pri_epoch");  // days since 1970 epoch
        sdmcInputs.priSat.posVel[0] = get_double_value(stateInfo, 0, "pri_pos_x");  // km
        sdmcInputs.priSat.posVel[1] = get_double_value(stateInfo, 0, "pri_pos_y");  // km
        sdmcInputs.priSat.posVel[2] = get_double_value(stateInfo, 0, "pri_pos_z");  // km
        sdmcInputs.priSat.posVel[3] = get_double_value(stateInfo, 0, "pri_vel_x");  // km/s
        sdmcInputs.priSat.posVel[4] = get_double_value(stateInfo, 0, "pri_vel_y");  // km/s
        sdmcInputs.priSat.posVel[5] = get_double_value(stateInfo, 0, "pri_vel_z");  // km/s
        sdmcInputs.priSat.cov[0]    = get_double_value(stateInfo, 0, "pri_cov_11"); // m^2
        sdmcInputs.priSat.cov[1]    = get_double_value(stateInfo, 0, "pri_cov_21"); // m^2
        sdmcInputs.priSat.cov[2]    = get_double_value(stateInfo, 0, "pri_cov_22"); // m^2
        sdmcInputs.priSat.cov[3]    = get_double_value(stateInfo, 0, "pri_cov_31"); // m^2
        sdmcInputs.priSat.cov[4]    = get_double_value(stateInfo, 0, "pri_cov_32"); // m^2
        sdmcInputs.priSat.cov[5]    = get_double_value(stateInfo, 0, "pri_cov_33"); // m^2
        sdmcInputs.priSat.cov[6]    = get_double_value(stateInfo, 0, "pri_cov_41"); // m^2/s
        sdmcInputs.priSat.cov[7]    = get_double_value(stateInfo, 0, "pri_cov_42"); // m^2/s
        sdmcInputs.priSat.cov[8]    = get_double_value(stateInfo, 0, "pri_cov_43"); // m^2/s
        sdmcInputs.priSat.cov[9]    = get_double_value(stateInfo, 0, "pri_cov_44"); // m^2/s^2
        sdmcInputs.priSat.cov[10]   = get_double_value(stateInfo, 0, "pri_cov_51"); // m^2/s
        sdmcInputs.priSat.cov[11]   = get_double_value(stateInfo, 0, "pri_cov_52"); // m^2/s
        sdmcInputs.priSat.cov[12]   = get_double_value(stateInfo, 0, "pri_cov_53"); // m^2/s
        sdmcInputs.priSat.cov[13]   = get_double_value(stateInfo, 0, "pri_cov_54"); // m^2/s^2
        sdmcInputs.priSat.cov[14]   = get_double_value(stateInfo, 0, "pri_cov_55"); // m^2/s^2
        sdmcInputs.priSat.cov[15]   = get_double_value(stateInfo, 0, "pri_cov_61"); // m^2/s
        sdmcInputs.priSat.cov[16]   = get_double_value(stateInfo, 0, "pri_cov_62"); // m^2/s
        sdmcInputs.priSat.cov[17]   = get_double_value(stateInfo, 0, "pri_cov_63"); // m^2/s
        sdmcInputs.priSat.cov[18]   = get_double_value(stateInfo, 0, "pri_cov_64"); // m^2/s^2
        sdmcInputs.priSat.cov[19]   = get_double_value(stateInfo, 0, "pri_cov_65"); // m^2/s^2
        sdmcInputs.priSat.cov[20]   = get_double_value(stateInfo, 0, "pri_cov_66"); // m^2/s^2
        sdmcInputs.priSat.dcpPosVelSensitivity[0] = get_double_value(stateInfo, 0, "pri_pos_sensitivity_x"); // m
        sdmcInputs.priSat.dcpPosVelSensitivity[1] = get_double_value(stateInfo, 0, "pri_pos_sensitivity_y"); // m
        sdmcInputs.priSat.dcpPosVelSensitivity[2] = get_double_value(stateInfo, 0, "pri_pos_sensitivity_z"); // m
        sdmcInputs.priSat.dcpPosVelSensitivity[3] = get_double_value(stateInfo, 0, "pri_vel_sensitivity_x"); // m/s
        sdmcInputs.priSat.dcpPosVelSensitivity[4] = get_double_value(stateInfo, 0, "pri_vel_sensitivity_y"); // m/s
        sdmcInputs.priSat.dcpPosVelSensitivity[5] = get_double_value(stateInfo, 0, "pri_vel_sensitivity_z"); // m/s
        sdmcInputs.priSat.dcpDensityForecastUncertainty = get_double_value(stateInfo, 0, "pri_forecast_uncertainty"); // unitless
        
        // Initialize the information for the secondary satellite
        sdmcInputs.secSat.satno     = get_int_value(stateInfo, 0, "sec_satno");
        sdmcInputs.secSat.epoch     = get_double_value(stateInfo, 0, "sec_epoch");  // days since 1970 epoch
        sdmcInputs.secSat.posVel[0] = get_double_value(stateInfo, 0, "sec_pos_x");  // km
        sdmcInputs.secSat.posVel[1] = get_double_value(stateInfo, 0, "sec_pos_y");  // km
        sdmcInputs.secSat.posVel[2] = get_double_value(stateInfo, 0, "sec_pos_z");  // km
        sdmcInputs.secSat.posVel[3] = get_double_value(stateInfo, 0, "sec_vel_x");  // km/s
        sdmcInputs.secSat.posVel[4] = get_double_value(stateInfo, 0, "sec_vel_y");  // km/s
        sdmcInputs.secSat.posVel[5] = get_double_value(stateInfo, 0, "sec_vel_z");  // km/s
        sdmcInputs.secSat.cov[0]    = get_double_value(stateInfo, 0, "sec_cov_11"); // m^2
        sdmcInputs.secSat.cov[1]    = get_double_value(stateInfo, 0, "sec_cov_21"); // m^2
        sdmcInputs.secSat.cov[2]    = get_double_value(stateInfo, 0, "sec_cov_22"); // m^2
        sdmcInputs.secSat.cov[3]    = get_double_value(stateInfo, 0, "sec_cov_31"); // m^2
        sdmcInputs.secSat.cov[4]    = get_double_value(stateInfo, 0, "sec_cov_32"); // m^2
        sdmcInputs.secSat.cov[5]    = get_double_value(stateInfo, 0, "sec_cov_33"); // m^2
        sdmcInputs.secSat.cov[6]    = get_double_value(stateInfo, 0, "sec_cov_41"); // m^2/s
        sdmcInputs.secSat.cov[7]    = get_double_value(stateInfo, 0, "sec_cov_42"); // m^2/s
        sdmcInputs.secSat.cov[8]    = get_double_value(stateInfo, 0, "sec_cov_43"); // m^2/s
        sdmcInputs.secSat.cov[9]    = get_double_value(stateInfo, 0, "sec_cov_44"); // m^2/s^2
        sdmcInputs.secSat.cov[10]   = get_double_value(stateInfo, 0, "sec_cov_51"); // m^2/s
        sdmcInputs.secSat.cov[11]   = get_double_value(stateInfo, 0, "sec_cov_52"); // m^2/s
        sdmcInputs.secSat.cov[12]   = get_double_value(stateInfo, 0, "sec_cov_53"); // m^2/s
        sdmcInputs.secSat.cov[13]   = get_double_value(stateInfo, 0, "sec_cov_54"); // m^2/s^2
        sdmcInputs.secSat.cov[14]   = get_double_value(stateInfo, 0, "sec_cov_55"); // m^2/s^2
        sdmcInputs.secSat.cov[15]   = get_double_value(stateInfo, 0, "sec_cov_61"); // m^2/s
        sdmcInputs.secSat.cov[16]   = get_double_value(stateInfo, 0, "sec_cov_62"); // m^2/s
        sdmcInputs.secSat.cov[17]   = get_double_value(stateInfo, 0, "sec_cov_63"); // m^2/s
        sdmcInputs.secSat.cov[18]   = get_double_value(stateInfo, 0, "sec_cov_64"); // m^2/s^2
        sdmcInputs.secSat.cov[19]   = get_double_value(stateInfo, 0, "sec_cov_65"); // m^2/s^2
        sdmcInputs.secSat.cov[20]   = get_double_value(stateInfo, 0, "sec_cov_66"); // m^2/s^2
        sdmcInputs.secSat.dcpPosVelSensitivity[0] = get_double_value(stateInfo, 0, "sec_pos_sensitivity_x"); // m
        sdmcInputs.secSat.dcpPosVelSensitivity[1] = get_double_value(stateInfo, 0, "sec_pos_sensitivity_y"); // m
        sdmcInputs.secSat.dcpPosVelSensitivity[2] = get_double_value(stateInfo, 0, "sec_pos_sensitivity_z"); // m
        sdmcInputs.secSat.dcpPosVelSensitivity[3] = get_double_value(stateInfo, 0, "sec_vel_sensitivity_x"); // m/s
        sdmcInputs.secSat.dcpPosVelSensitivity[4] = get_double_value(stateInfo, 0, "sec_vel_sensitivity_y"); // m/s
        sdmcInputs.secSat.dcpPosVelSensitivity[5] = get_double_value(stateInfo, 0, "sec_vel_sensitivity_z"); // m/s
        sdmcInputs.secSat.dcpDensityForecastUncertainty = get_double_value(stateInfo, 0, "sec_forecast_uncertainty"); // unitless
        
        // Initialize general conjunction inputs
        sdmcInputs.tca                = get_double_value(stateInfo, 0, "tca");              // days since 1970 epoch
        sdmcInputs.trajectoryMode     = get_int_value(stateInfo, 0,    "trajectory_mode");  // 0 = 2-body, 1 = rectilinear, 2 = rectiliniear (position deviations only)
        sdmcInputs.span               = get_double_value(stateInfo, 0, "span_days");        // days
        sdmcInputs.seed               = get_int_value(stateInfo, 0,    "seed");
        sdmcInputs.numTrials          = get_int_value(stateInfo, 0,    "num_trials");
        sdmcInputs.radius             = get_double_value(stateInfo, 0, "hbr_m");            // m
        sdmcInputs.maxRadius          = get_double_value(stateInfo, 0, "max_radius");       // m
        sdmcInputs.maxNumOutputTrials = get_int_value(stateInfo, 0,    "max_num_output_trials");
        string listFileName           = get_string_value(stateInfo, 0, "sdmc_list_file");
        if (listFileName.length() >= 255) {
            for (int i = 0; i < 256; i++) {
                sdmcInputs.listFileName[i] = ' ';
            }
        } else {
            strncpy(sdmcInputs.listFileName, listFileName.c_str(), listFileName.length());
            for (int i = listFileName.length(); i < 256; i++) {
                sdmcInputs.listFileName[i] = ' ';
            }
        }
        
        sdmcInputsInitialized = true;
    }
    
    // Frees the dynamically allocated memory used within the sdmcOutputs
    // data strcuture.
    void reset_sdmc_outputs(void) {
        if (sdmcOutputs.data != NULL) {
            delete[] sdmcOutputs.data;
            sdmcOutputs.data = NULL;
        }
    }
    
    // Initializes the sdmcOutputs structure by allocating memory according
    // to sizes specified in stateInfo. Initializes all values to 0.
    //   Inputs:
    //     stateInfo - Pointer to a mex structure representing a Matlab
    //                 structure array which contains conjunction state
    //                 information. The stateInfo structure is populated by
    //                 the call_SDMC.m Matlab routine.
    void init_sdmc_outputs(const mxArray *stateInfo) {
        sdmcOutputsInitialized = false;
        
        // Initialize the general outputs
        sdmcOutputs.numHits = 0;
        sdmcOutputs.numOutputTrials = 0;
        
        // Initialize the variable size outputs
        reset_sdmc_outputs();
        int nDataPoints = get_int_value(stateInfo, 0, "max_num_output_trials");
        sdmcOutputs.data = new SdmcOutputRecord[nDataPoints];
        for (int i = 0; i < nDataPoints; i++) {
            sdmcOutputs.data[i].hitIndicator = 0.0;
            sdmcOutputs.data[i].hitTime      = 0.0;
            sdmcOutputs.data[i].hitTimeMiss  = 0.0;
            sdmcOutputs.data[i].hitRadius    = 0.0;
            sdmcOutputs.data[i].pcaTime      = 0.0;
            sdmcOutputs.data[i].pcaRadius    = 0.0;
            for (int j = 0; j < POS_VEL_SIZE; j++) {
                sdmcOutputs.data[i].priPosVel[j] = 0.0;
                sdmcOutputs.data[i].secPosVel[j] = 0.0;
            }
        }
        
        sdmcOutputsInitialized = true;
    }
    
public:
    // Constructor
    matlab_sdmc_wrapper()
    {
        errMsg = "";
        
        sdmcInputsInitialized  = false;
        sdmcOutputsInitialized = false;
        
        sdmcOutputs.data = NULL;
        reset_sdmc_outputs();
    }
    
    // Destructor
    ~matlab_sdmc_wrapper()
    {
        reset_sdmc_outputs();
    }
    
    // Runs the SDMC library
    int run_sdmc(const mxArray *stateInfo)
    {
        init_sdmc_inputs(stateInfo);
        init_sdmc_outputs(stateInfo);
        
        if (!sdmcInputsInitialized || !sdmcOutputsInitialized) {
            errMsg = "Errors encountered initializing the state when running run_sdmc().";
            mexErrMsgTxt(errMsg.c_str());
        }
        
        return RunSdmc(&sdmcInputs, &sdmcOutputs);
    }
    
    // Returns SDMC outputs
    SdmcOutputs* get_sdmc_outputs(void) {
        return &sdmcOutputs;
    }
};
#endif // __MATLAB_SDMC_WRAPPER_HPP__

// ----------------- END OF CODE ------------------
//
// Please record any changes to the software in the change history 
// shown below:
//
// ----------------- CHANGE HISTORY ------------------
// Developer      |    Date     |     Description
// ---------------------------------------------------
// L. Baars       | 2023-Mar-13 | Initial release
// L. Baars       | 2025-Aug-06 | Minor documentation updates necessary for
//                                public release

//====================================================================================
//
// Copyright (c) 2023-2025 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//
//====================================================================================
