//====================================================================================
//
// Copyright (c) 2018, Oliver Woodford
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in
//       the documentation and/or other materials provided with the distribution
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
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
//   Mex C++ wrapper which creates an interface between the functions and structures
//   defined in SdmcInterface.h and Matlab.
//
//   Source code originated from Oliver Woodford's mex_class_wrapper GitHub project.
//
//====================================================================================
//
// References:
//
//   Oliver Woodford (2023). Example MATLAB class wrapper for a C++ class
//   (https://github.com/ojwoodford/mex_class_wrapper/releases/tag/v1.4.1),
//   GitHub. Retrieved March 13, 2023.
//
//====================================================================================
//
// Initial version: 2018; Latest update: Aug 2025
//
// ----------------- BEGIN CODE -----------------

#include "mex.h"
#include "matrix.h"
#include "class_handle.hpp"
#include "matlab_sdmc_wrapper.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get the command string
    char cmd[64];
	if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");
        
    // New
    if (!strcmp("new", cmd)) {
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("New: One output expected.");
        // Return a handle to a new C++ instance
        plhs[0] = convertPtr2Mat<matlab_sdmc_wrapper>(new matlab_sdmc_wrapper);
        return;
    }
    
    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2)
		mexErrMsgTxt("Second input should be a class instance handle.");
    
    // Delete
    if (!strcmp("delete", cmd)) {
        // Destroy the C++ object
        destroyObject<matlab_sdmc_wrapper>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }
    
    // Get the class instance pointer from the second input
    matlab_sdmc_wrapper* matlab_sdmc_wrapper_instance = convertMat2Ptr<matlab_sdmc_wrapper>(prhs[1]);
    
    // Call the various class methods
    // run_sdmc
    if (!strcmp("run_sdmc", cmd)) {
        // Check parameters
        if (nlhs > 1 || nrhs != 3)
            mexErrMsgTxt("run_sdmc: Unexpected arguments.");
        // Check the first argument, which should be the StateInfo struct
        // which should only have one row in it
        if (mxIsStruct(prhs[2]) == 0 || mxGetNumberOfElements(prhs[2]) != 1) {
            mexErrMsgTxt("run_sdmc: First argument must be a singular StateInfo struct.");
        }
        // Call the method
        int rtrnVal = matlab_sdmc_wrapper_instance->run_sdmc(prhs[2]);
        
        // No return values
        if (nlhs == 0)
            return;
        
        // Otherwise, return the return value
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
        mxDouble* pointer = mxGetDoubles(plhs[0]);
#else
        double* pointer = mxGetPr(plhs[0]);
#endif
        pointer[0] = rtrnVal;
        return;
    }
    // get_sdmc_outputs
    if (!strcmp("get_sdmc_outputs", cmd)) {
        // Check parameters
        if (nlhs > 2 || nrhs > 2)
            mexErrMsgTxt("get_sdmc_outputs: Unexpected arguments.");
        if (nlhs == 0)
            return;
        
        // Verify that there is data to return
        SdmcOutputs* sdmcOutputsPtr = matlab_sdmc_wrapper_instance->get_sdmc_outputs();
        int numRows = sdmcOutputsPtr->numOutputTrials;
        
        // Setup an array to return the data
        if (nlhs >= 1) {
            plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
            mxDouble* pointer = mxGetDoubles(plhs[0]);
#else
            double* pointer = mxGetPr(plhs[0]);
#endif
            pointer[0] = sdmcOutputsPtr->numHits;
        }
        if (nlhs == 2) {
            plhs[1] = mxCreateNumericMatrix(numRows, 18, mxDOUBLE_CLASS, mxREAL);
            
#if MX_HAS_INTERLEAVED_COMPLEX
            mxDouble* pointer = mxGetDoubles(plhs[1]);
#else
            double* pointer = mxGetPr(plhs[1]);
#endif
        
            // Fill the array with the returned data
            int numPosVel = 6;
            for (int i = 0; i < numRows; i++) {
                pointer[i]            = sdmcOutputsPtr->data[i].hitIndicator;
                pointer[i+numRows]    = sdmcOutputsPtr->data[i].hitTime;
                pointer[i+numRows*2]  = sdmcOutputsPtr->data[i].hitTimeMiss;
                pointer[i+numRows*3]  = sdmcOutputsPtr->data[i].hitRadius;
                pointer[i+numRows*4]  = sdmcOutputsPtr->data[i].pcaTime;
                pointer[i+numRows*5]  = sdmcOutputsPtr->data[i].pcaRadius;
                for (int j = 0; j < numPosVel; j++) {
                    pointer[i+numRows*(j+numPosVel)]   = sdmcOutputsPtr->data[i].priPosVel[j];
                    pointer[i+numRows*(j+2*numPosVel)] = sdmcOutputsPtr->data[i].secPosVel[j];
                }
            }
        }
        
        return;
    }
    
    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}

// ----------------- END OF CODE ------------------
//
// Please record any changes to the software in the change history 
// shown below:
//
// ----------------- CHANGE HISTORY ------------------
// Developer      |    Date     |     Description
// ---------------------------------------------------
// L. Baars       | 2023-Mar-13 | Copied code from GitHub site and modified
//                                the code to run with the SdmcInterface.h
//                                definitions.
// L. Baars       | 2025-Aug-06 | Minor documentation updates necessary for
//                                public release.

//====================================================================================
//
// Copyright (c) 2023-2025 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//
//====================================================================================
