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
//   Header definitions used by sdmcEntry.cpp in order to interface directly with
//   the SDMC library.
//
//====================================================================================
//
// Initial version: Mar 2023; Latest update: Aug 2025
//
// ----------------- BEGIN CODE -----------------

#ifndef SDMC_ENTRY__H
#define SDMC_ENTRY__H

#include "SdmcInterface.h"

#define NSTATE  9
#define NPOSVEL 6

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32) || defined(_WIN64)
void SDMC_ENTRY(
#else
void sdmc_entry_(
#endif
                 double   pri_rep_eci[NSTATE+1],
                 double   sec_rep_eci[NSTATE+1],
                 double   pri_cov_uvw[NSTATE][NSTATE],
                 double   sec_cov_uvw[NSTATE][NSTATE],
                 double   pri_dsv_uvw[NPOSVEL+1],
                 double   sec_dsv_uvw[NPOSVEL+1],
                 int    * pri_satno  ,
                 int    * sec_satno  ,
                 double * ca_time_inp,
                 double * span_inp   ,
                 int    * seed_inp   ,
                 int    * num_trials_inp,
                 double * radius_inp ,
                 double * maxrad_inp ,
                 double * cp_fact_inp,
                 int    * rectopt_inp,
                 int    * ds_size_inp,
                 double   dataset_out[],
                 int    * count_out,
                 int    * count_hits,
                 char     output_file[256],
                 int    * return_code);

#ifdef __cplusplus
}
#endif

#endif

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
