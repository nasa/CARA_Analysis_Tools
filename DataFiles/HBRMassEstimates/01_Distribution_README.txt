This file provides some helpful information on the distribution and use of
NASA CARA hard body radius (HBR) and mass estimates. In addition, it provides
some instructions on merging NASA CARA HBR and mass estimates with ESA DISCOS
size and mass estimates.

================================================================================
Background:

The NASA CARA team uses RCS observation data from the Space Fence Kwajalein
(SFK) phased array sensor in order to estimate the size and mass of unknown
resident space objects. A rigorous filtering is performed which cleans the RCS
observation data from possible calibration issues and/or large departures from
monthly mean RCS data. Using NASA's Size Estimation Model (SEM), the cleaned
RCS data can then be used to estimate the sizes and masses of objects.

In addition to the RCS estimated outputs, NASA CARA has developed a lookup table
which includes CARA protected satellite values and European Space Agency (ESA)
Database and Information System Characterising Objects in Space (DISCOS)*
catalogs. If a resident space object does not appear in any of the
aforementioned lists, then default HBR and/or mass estimates are supplied.
Values within the lookup table are provided in descending order according to the
following priorities:
  1. CARA Protected Primary
  2. DISCOS DB
  3. SFK RCS Estimate
  4. Default

* Note: Due to data licensing constraints, NASA CARA cannot provide DISCOS data
        directly. However, a tool has been provided which will allow a user to
        download the data from the DISCOS website (requires setup of a user
        account) and merge the DISCOS data into the HBR and Mass lookup table.

For more information about the RCS estimation process, please review the
following references:

Baars, L., and Hall, D., "Processing Space Fence Radar Cross-Section Data to
  Produce Size and Mass Estimates," AAS Astrodynamics Specialist Conference,
  AAS Paper 22-586, Aug. 2022.

Hall, D., and Baars, L., "Satellite Collision and Fragmentation Probabilities
  Using Radar-Based Size and Mass Estimates," Journal of Spacecraft and Rockets,
  Volume 60, Issue 4, Jul. 2023, pg. 1319-1332.

================================================================================
Merging DISCOS Data:

As mentioned earlier, DISCOS data is not provided due to licensing constraints.
However, NASA CARA has developed a Matlab tool, MergDiscosData.m, which will
download sizing data from the DISCOS database and merge the data into the HBR
and mass estimates file. The files and directories required to run the tool
include:

  - MergeDiscosData.m - Driver function which downloads and merges DISCOS data
                        into the CSV file passed in.

  - connect_params_distrib.m - Contains connection parameters for the DISCOSweb
                               API, including personal access tokens to login
                               to the API.

  - GetDiscosData.m - Actual software which downloads DISCOS data and calculates
                      HBR sizes based on shape information provided by DISCOS.

  - DummyPasswordManager - Directory containing functions for a dummy password
                           manager.

    - GetConnectionIdx.m - Retrieves the connection index from
                           connect_params_distrib for the connection name passed
                           in (for cases where multiple connections are defined
                           within the file).

    - GetSecKey.m - Dummy function for getting a secret key.

    - GetPassword.m - Dummy function for retrieving a password defined by a
                      connection parameters file.

Note: GetDiscosData.m is a function used by NASA CARA to retrieve DISCOS data
      via automated processes. That software uses a password manager with
      encryption techonology in order to prevent the storage of plain-text
      passwords onto the file system. In order to prevent possible ITAR
      restrictions, NASA CARA has replaced the dependency on encrypted password
      storage with the supplied DummyPasswordManager. It is up to the end-user
      to implement their own encrypted password storage mechanism if it is
      required on their system.

--------------------------------------------------------------------------------
Setup Instructions:

1. Go to the DISCOSweb page and register for a new account:
   https://discosweb.esoc.esa.int/

2. Once your DISCOSweb account is approved, sign into your account, click on
   your name and then click on the "Access tokens" menu option.

   a. On the Personal access tokens page, click on the "Generate new token"
      button.

   b. Give the token a name and click on the "Generate token" button.

   c. *IMPORTANT* Copy the token displayed and save it off. The token will only
      be displayed this one time and cannot be retrieved once you navigate away
      from the page. If you lose this token you will have to revoke the token
      and generate a new one.

      Do not share your token with anyone else. Multiple people using the same
      token can cause "TOO MANY REQUESTS" errors discussed below.

3. Open the connect_params_distrib.m file and replace the <replace with token>
   text with the token you copied from the DISCOSweb page. Save the
   connect_params_distrib.m file.

Your installation of the tool is now ready to download and merge DISCOS data.
The steps in this section should only have to be run once.

--------------------------------------------------------------------------------
Running the Merge Tool:

1. Download the latest HBR_Mass_Estimates_yyyymmdd_ForDistribution.csv file from
   the NASA CARA GitHub site.

2. Open Matlab and run MergeDiscosData.m with the HBR/Mass estimates file name
   (full or relative path is OK) as the only parameter.

3. If no DISCOS data CSV file exists with the same "yyyymmdd" date stamp, the
   tool will download the DISCOS data into a file named
   discos_data_yyyymmdd.csv. The file will exist in the same directory as the
   HBR/Mass estimates file.

   Note: This process will take a long time (~30 minutes in late 2023) due to
         data throughput constraints imposed by the ESA website. If too much
         data is requested within a minute, the website will respond with a
         "TOO MANY REQUESTS" response and will not allow the token to be used
         until the request timeout has occurred (usually 1 minute). In order to
         prevent this from happening, the tool will automatically pause
         processing for 1 minute after 20 pages of data have been requested. It
         will then resume downloading another 20 pages of data and pause again.
         This will continue until all pages have been downloaded from the DISCOS
         website.

4. After the DISCOS data has been downloaded, the tool will merge the DISCOS
   data into the HBR/Mass estimates file and create a file named
   HBR_Mass_Estimates_yyyymmdd_Merged_NotForDistribution.csv. This file will be
   written to the same directory as the original HBR/Mass inputs file.

   Note: As part of signing up for the DISCOSweb account, the end-user falls
         under the same licensing agreement for using and distributing DISCOS
         data. The end-user is not allowed to redistribute the merged data,
         hence the "NotForDistribution" tag that has been added to the file
         name.

   a. The outputs from the Merge process should look something like the
      following:

>> MergeDiscosData('HBR_Mass_Estimates_20231101_ForDistribution.csv')
Downloading DISCOS data...success
Saved DISCOS data into file: discos_data_20231101.csv
Elapsed time is 2116.035076 seconds.
Merged 24383 HBR values from the DISCOS DB
Merged 25241 mass values from the DISCOS DB
Set 0 default HBR values due to missing or NaN DISCOS DB values
Set 0 default mass values due to missing or NaN DISCOS DB values
Wrote merged data into the file HBR_Mass_Estimates_20231101_Merged_NotForDistribution.csv

      The "Set x default ..." lines indicate instances where the original CSV
      file indicated that DISCOS data was supposed to be applied, but the query
      from the DISCOS website did not retrieve valid HBR or mass data for the
      object. This will occur more often when the date that the merge tool was
      run is significantly different from the date stamp of the HBR and mass
      estimates file.

This merge process should nominally be run once a month whenever a new HBR and
mass estimates file has been downloaded from the NASA CARA GitHub site.

================================================================================
Getting HBR and Mass Estimates:

The CARA NASA team has provided a convenient Matlab tool for retrieving HBR and
Mass estimates from the HBR/Mass estimates file, GetHBRMassEstimate.m. Simply
pass in the name of the merged HBR/Mass estimates file and a NORAD catalog ID
for the object in question and the tool will return HBR, HBRSigma, HBRSource,
Mass, MassSigma, and MassSource values for the object. See the documentation
within the tool for more information on the input and output values.

Possible data sources returned by the tool are the data sources described in the
background section above.

A particular use case would be to use the outputs from this tool in order to
supply the information needed to calculate the effective HBR as described in
equation 8 of [Baars and Hall, 2022]. The equation is provided here for
convenience.

effHBR = sqrt((Rbar_1 + Rbar_2)^2 + sigma_R1^2 + sigma_R2^2)

The HBR and HBRSigma values returned from the tool can be used as the Rbar and
sigma_R values in the equation, respectively. Specifically, the first call to
the tool would be for the primary satellite involved in the conjunction, thus
supplying Rbar_1 and sigma_R1 values. The second call to the tool would be for
the secondary satellite, providing the Rbar_2 and sigma_R2 values. Notably, if
the sigma value for the primary is 0, equation 8 reduces into equation 9 from
the same paper. If both the primary and secondary sigma values are 0, then the
equation reduces down to a simple combined HBR, which is the sum of the primary
HBR and secondary HBR.
