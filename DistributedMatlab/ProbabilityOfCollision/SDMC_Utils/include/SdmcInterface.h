#ifndef SDMC_INTERFACE__H
#define SDMC_INTERFACE__H

//==============================================================================
// Copyright (c) 2023-2025 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//==============================================================================
// SdmcInterface.h
//
// Defines the structures and functions used to call SDMC as a library. This C
// header file definition has been structured such that it can be used with C
// and C++ linkages. See the embedded comments for a full definition of the
// inputs/outputs used by the library.
//==============================================================================

//==============================================================================
// Constants Definitions
//==============================================================================
// Hardcoded defaults for position, velocity, and covariance sizes

#define POS_VEL_SIZE  6
#define COV_SIZE      21

// Covariance definition, units = (m^2, m^2/s, m^2/s^2):
// [CR_R
//  CT_R     CT_T
//  CN_R     CN_T     CN_N
//  CRDOT_R  CRDOT_T  CRDOT_N  CRDOT_RDOT
//  CTDOT_R  CTDOT_T  CTDOT_N  CTDOT_RDOT  CTDOT_TDOT
//  CNDOT_R  CNDOT_T  CNDOT_N  CNDOT_RDOT  CNDOT_TDOT  CNDOT_NDOT]
//
// Component definitions:
// R = Radial
// T = Transverse Direction (a.k.a. In-Track)
// N = Normal Direction (a.k.a. Cross-Track)
// RDOT = Radial Inertial Velocity (vs relative to rotating RTN frame)
// TDOT = Transverse Inertial Velocity (vs relative to rotating RTN frame)
// NDOT = Normal Inertial Velocity (vs relative to rotating RTN frame)

//==============================================================================
// Data Structure Definitions
//==============================================================================

struct SatelliteInfo {
    int    satno;                              // Catalog satellite ID
    double epoch;                              // Pos/Vel/cov/dsv epoch time
                                               // (days since Jan 0.0 1970 UTC)
    double posVel[POS_VEL_SIZE];               // Inertial Pos/Vel at epoch (km, km/sec)
    double cov[COV_SIZE];                      // Covariance at epoch (see cov definition above)
    double dcpPosVelSensitivity[POS_VEL_SIZE]; // DCP RTN sensitivity vectors (m, m/s)*
    double dcpDensityForecastUncertainty;      // DCP density forecast uncertainty per a
                                               // relative density error sigma (unitless)*
};
// *Note: Zero values can be entered for these variables. When this occurs,
//        cross-correlation effects will not be included within the SDMC
//        calculation.

struct SdmcInputs {
    struct SatelliteInfo priSat;             // Primary satellite info at TCA
    struct SatelliteInfo secSat;             // Secondary satellite info at TCA
    double               tca;                // TCA (days since Jan 0.0 1970 epoch UTC)
    int                  trajectoryMode;     // 0 = 2-body, 1 = full rectilinear motion,
                                             // 2 = rectilinear motion (position deviations only)
    double               span;               // SDMC evaluation span (in days)
    int                  seed;               // random number generator seed (negative integer)
    int                  numTrials;          // number of SDMC trials to run (positive integer)
    double               radius;             // hard-body radius of a hit (m)
    double               maxRadius;          // maximum reporting radius (m), negative value
                                             // will be interpreted as abs(maxRadius) and will
                                             // cause priPosVel and secPosVel in SdmcOutputs
                                             // struct to be fully filled in
    int                  maxNumOutputTrials; // maximum number of trials which are allowed to
                                             // be reported in the SdmcOutputs.data array
    char                 listFileName[256];  // file for basic output (list) information
                                             // no output for first character blank or null
                                             // screen output if first 6 characters are STDOUT
                                             // file output otherwise to listFileName
};

struct SdmcOutputRecord {
    double hitIndicator;            // 1 = hit, 0 = miss (double to allow large Fortran array)
    double hitTime;                 // time of impingement (for hit) or global PCA** (for miss)
    double hitTimeMiss;             // time of hit miss distance or global PCA
    double hitRadius;               // hit miss distance or global PCA distance (m)
    double pcaTime;                 // time of global PCA
    double pcaRadius;               // global PCA distance (m)
    double priPosVel[POS_VEL_SIZE]; // Primary   sat inertial pos/vel (km,s) at hitTimeMiss above
    double secPosVel[POS_VEL_SIZE]; // Secondary sat inertial pos/vel (km,s) at hitTimeMiss above
};
// **Note: PCA = Point of Closest Approach and all times in structure are seconds since TCA

struct SdmcOutputs {
    int    numHits;                // Number of SDMC trials which resulted in a
                                   // PCA distance less than the hard-body radius
    int    numOutputTrials;        // Number of filled in records in the array
    struct SdmcOutputRecord *data; // Array of output records. The length of
                                   // allocated memory for this array should equal
                                   // SdmcInputs.maxNumOutputTrials.
};


//==============================================================================
// Function Definitions
//==============================================================================

// This extern definition allows the RunSdmc() function to be called from both
// C and C++ linkages
#ifdef __cplusplus
extern "C" {
#endif

// Runs Two-Body Monte Carlo (SDMC) for the number of trials indicated within
// the SdmcInputs structure. A data row in the SdmcOutputs structure will be
// generated for first N trials which have squared PCA distance less than or equal
// to the squared maxRadius. If the maxRadius is less than 0, then priPosVel and
// secPosVel will be populated with data. The maximum size (N) for numRecords will
// be SdmcInputs.maxNumOutputTrials; when allocating memory for SdmcOutputs, this
// should be the length used to allocate the data array. It is up to the caller
// of the RunSdmc() function to allocate/deallocate memory within SdmcOutputs.
//
// Inputs:
//    SdmcInputs - Pointer to structure containing input variables for running
//                 SDMC
//    SdmcOutputs - Pointer to structure containing allocated memory to store
//                  any outputs generated by SDMC
// Outputs:
//    runStatus - Indicates if any fatal errors where encountered while running
//                the RunSdmc() function.
//                0 = success, postive value = failure, negative value = warning
//                1 = Hyperbolic primary epoch elements
//                2 = Rectilinear primary epoch elements
//                3 = Retrograde equatorial primary epoch elements
//                4 = Hyperbolic secondary epoch elements
//                5 = Rectilinear secondary epoch elements
//                6 = Retrograde equatorial secondary epoch elements
//                7 = Max internal ephemeris points exceeded (should not occur)
//                8 = Too few internal ephemeris points (5)  (should not occur)
//                9 = Too many bad elements during MC trials
//               10 = Could not open output list file
//               -1 = NPD encountered for primary with remediation
//               -2 = NPD encountered for secondary with remediation
//               -3 = Small TCA relative velocity condition
//               -4 = Small epoch covariance condition
//               -5 = Bad elements encountered during MC trials

int RunSdmc(struct SdmcInputs  *sdmcInputs,
            struct SdmcOutputs *sdmcOutputs);

// Closes the extern definition
#ifdef __cplusplus
}
#endif

#endif

//==============================================================================
// Please record any changes to the software in the change history 
// shown below:
//
// ----------------- CHANGE HISTORY ------------------
// Developer      |    Date     |     Description
// ---------------------------------------------------
// L. Baars       | 2023-Feb-10 | Initial Development
// L. Baars       | 2025-Aug-06 | Minor updates to comments necessary for public
//                                release.
//==============================================================================
// Copyright (c) 2023-2025 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//==============================================================================
