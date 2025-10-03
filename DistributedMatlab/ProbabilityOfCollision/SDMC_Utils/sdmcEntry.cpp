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
//   Implements the RunSdmc method in C++ to call the sdmc_entry_() method supplied
//   by the Fortran SDMC library.
//
//====================================================================================
//
// Initial version: Mar 2023; Latest update: Aug 2025
//
// ----------------- BEGIN CODE -----------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sdmcEntry.h"
#include "SdmcInterface.h"

int RunSdmc (struct SdmcInputs  * sdmcInputs,
             struct SdmcOutputs * sdmcOutputs)
{
    // Create variables to be passed into sdmc_entry_(). The datatypes and
    // sizes must match their Fortran equivalents.
    
    // Double scalar and array inputs
    double pri_rep_eci[NSTATE+1];
    double sec_rep_eci[NSTATE+1];
    double pri_cov_uvw[NSTATE][NSTATE];
    double sec_cov_uvw[NSTATE][NSTATE];
    double pri_dsv_uvw[NPOSVEL+1];
    double sec_dsv_uvw[NPOSVEL+1];
    double ca_time = 0.0;
    double ca_span = 0.0;
    double radius  = 0.0;
    double maxrad  = 0.0;
    double cp_fact = 0.0;
    int    dataSize = sdmcInputs->maxNumOutputTrials*18;
    double datapntr[dataSize];

    // Integer scalar inputs
    int pri_satno  = 0;
    int sec_satno  = 0;
    int rn_seed    = 0;
    int num_trials = 0;
    int rectopt    = 0;
    int rtncode    = 0;
    int ds_size    = 0;
    int count_out  = 0;
    int count_hits = 0;
    int i          = 0;
    int j          = 0;
    int k          = 0;

    // Character array inputs
    char output_file[256] = " ";
    strncpy(output_file,sdmcInputs->listFileName,256);

    // Initialze the state arrays
    for (i = 0; i < NSTATE+1; i++)
    {
        pri_rep_eci[i] = 0.0;
        sec_rep_eci[i] = 0.0;
    }

    // Initialize the covariance matrices
    for (i = 0; i < NSTATE; i++)
    {
        for (j = 0; j < NSTATE; j++)
        {
            pri_cov_uvw[i][j] = 0.0;
            sec_cov_uvw[i][j] = 0.0;
        }
    }

    // Initialize the covariance cross-correlation information
    for (i = 0; i < NPOSVEL+1; i++)
    {
        pri_dsv_uvw[i] = 0.0;
        sec_dsv_uvw[i] = 0.0;
    }

    // Initialize the output data structure
    for (i = 0; i < dataSize; i++)
    {
        datapntr[i] = 0.0;
    }

    // Copy the object IDs
    pri_satno = sdmcInputs->priSat.satno;
    sec_satno = sdmcInputs->secSat.satno;

    // Copy epoch info as the first element of the state arrays
    pri_rep_eci[0] = sdmcInputs->priSat.epoch;
    sec_rep_eci[0] = sdmcInputs->secSat.epoch;

    // Copy the density forecast uncertainty as the first element of the
    // covariance cross-correlation information
    pri_dsv_uvw[0] = sdmcInputs->priSat.dcpDensityForecastUncertainty;
    sec_dsv_uvw[0] = sdmcInputs->secSat.dcpDensityForecastUncertainty;

    // Copy the rest of the state array and covariance cross-correlation
    // information
    for (j = 1; j <= NPOSVEL; j++)
    {
        pri_rep_eci[j] = sdmcInputs->priSat.posVel[j-1];
        sec_rep_eci[j] = sdmcInputs->secSat.posVel[j-1];

        // Convert into km units
        pri_dsv_uvw[j] = sdmcInputs->priSat.dcpPosVelSensitivity[j-1] * 1.0E-3;
        sec_dsv_uvw[j] = sdmcInputs->secSat.dcpPosVelSensitivity[j-1] * 1.0E-3;
    }

    // Copy the covariance matrices
    for (j = 0, i = 0; j < NPOSVEL; j++)
    {
        for (k = 0; k <= j; k++, i++)
        {
            // Populate both the upper and lower triangular portions,
            // convert into km units
            pri_cov_uvw[j][k] = sdmcInputs->priSat.cov[i] * 1.0E-6;
            pri_cov_uvw[k][j] = sdmcInputs->priSat.cov[i] * 1.0E-6;

            sec_cov_uvw[j][k] = sdmcInputs->secSat.cov[i] * 1.0E-6;
            sec_cov_uvw[k][j] = sdmcInputs->secSat.cov[i] * 1.0E-6;
        }
    }

    // Fill in the rest of the data inputs
    ca_time    = sdmcInputs->tca;
    ca_span    = sdmcInputs->span;
    rn_seed    = sdmcInputs->seed;
    num_trials = sdmcInputs->numTrials;
    radius     = sdmcInputs->radius    * 1.0E-3; // Convert to km
    maxrad     = sdmcInputs->maxRadius * 1.0E-3; // Convert to km
    cp_fact    = 1.0;
    rectopt    = sdmcInputs->trajectoryMode;
    ds_size    = sdmcInputs->maxNumOutputTrials;

    // Call the Fortran routine
#if defined(_WIN32) || defined(_WIN64)
    SDMC_ENTRY(
#else
    sdmc_entry_(
#endif
                pri_rep_eci,  sec_rep_eci,
                pri_cov_uvw,  sec_cov_uvw,
                pri_dsv_uvw,  sec_dsv_uvw,
               &pri_satno  , &sec_satno  ,
               &ca_time    , &ca_span    ,
               &rn_seed    , &num_trials ,
               &radius     , &maxrad,
               &cp_fact    , &rectopt,
               &ds_size    ,  datapntr,
               &count_out  , &count_hits,
                output_file, &rtncode);

    // Fill in summary SDMC outputs
    sdmcOutputs->numOutputTrials = count_out;
    sdmcOutputs->numHits         = count_hits;
    
    // Convert the output array into the output data structure
    for (i = 0; i < count_out; i++)
    {
        sdmcOutputs->data[i].hitIndicator = datapntr[i*18];
        sdmcOutputs->data[i].hitTime      = datapntr[i*18+1];
        sdmcOutputs->data[i].hitTimeMiss  = datapntr[i*18+2];
        sdmcOutputs->data[i].hitRadius    = datapntr[i*18+3];
        sdmcOutputs->data[i].pcaTime      = datapntr[i*18+4];
        sdmcOutputs->data[i].pcaRadius    = datapntr[i*18+5];
        if (maxrad < 0) {
            for (j = 0; j < POS_VEL_SIZE; j++)
            {
                sdmcOutputs->data[i].priPosVel[j] = datapntr[i*18+j+6];
                sdmcOutputs->data[i].secPosVel[j] = datapntr[i*18+j+12];
            }
        }
    }

  return rtncode;
}

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
