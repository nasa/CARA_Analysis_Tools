# Testing Information
The `CARA_PcMethod_Test_Conjunctions.xlsx` file contains the results of running various probability of collision (Pc) Methods (Pc2D, Nc2D, Nc3D, and SDMCPc) against the CDMs provided in this directory. The methods map to actual Matlab functions as follows:
| Method     | Matlab Function | Released                     |
| ---------- | --------------- | ---------------------------- |
| Pc2D_NoAdj | PcCircle.m      | Yes                          |
| Pc2D       | PcCircle.m      | Yes                          |
| Nc2D       | Pc2D_Hall.m     | Yes                          |
| Nc3D       | Pc3D_Hall.m     | Yes                          |
| SDMCPc     | Pc_SDMC.m       | Yes (Windows and Linux only) |

As of October 2025, all of the methods have been approved for public release and are available in the `DistributedMatlab/ProbabilityOfCollision` folder.

We have found that Matlab outputs can differ based on the Matlab version being used, the AVX Math Kernel Library (MKL) flags used, and the Operating System used. It is important to be aware of all of these parameters when comparing values.

Newer versions of Matlab can and will use updated versions of the BLAS and LAPACK low-level Fortran libraries. These updates generally represent algorithm enhancements for improved accuracy and/or performance. However, these updates can also prove to be problematic for numerical validation of results, and should be taken into consideration when comparing against the numbers provided in the spreadsheet. Different Matlab versions can account for Pc tolerance differences near 1e-8 or below.

The AVX Math Kernel Library can be controlled by setting the MKL_ENABLE_INSTRUCTIONS environment variable and restarting Matlab. Similarly to the Matlab versions, different MKL settings can account for Pc tolerance differences near 1e-12 or below. However, this setting can only be used on Intel-based processors (sorry Mac Silicon users).

We've also encountered very slight computation differences depending on the operating system used. When testing across Windows and Linux, we've seen Pc tolerance differences near 1e-15 or below. Testing on Macs with the new Silicon processors have produced Pc tolerance differences near 1e-12 or below.

The BLAS and LAPACK build info parameters are generally reflections of the Matlab version and MKL settings. In the Matlab prompt, the `version -blas` and `version -lapack` commands will show the settings on your machine.

The results reported within `CARA_PcMethod_Test_Conjunctions.xlsx` were generated using the following configuration:
| Item Description                | Value                 |
| ------------------------------- | --------------------- |
| Last Update                     | 07/03/2025            |
| Matlab Version                  | R2019b                |
| MKL_ENABLE_INSTRUCTIONS setting | AVX2                  |
| Operating System                | Windows 11 Enterprise |
| BLAS Build Info                 | Intel(R) Math Kernel Library Version 2018.0.3 Product Build 20180406 for Intel(R) 64 architecture applications, CNR branch AVX2 |
| LAPACK Build Info               | Intel(R) Math Kernel Library Version 2018.0.3 Product Build 20180406 for Intel(R) 64 architecture applications, CNR branch AVX2 Linear Algebra PACKage Version 3.7.0 |

Be advised that differences in the settings listed above will likely cause Pc calculation differences on your system.

# Calculating Pc Values from CDM Files and Sample Scripts

The accurate calculation of the Pc from CDM files requires a careful setup of the parameters provided within the CDM file. However, a direct read and calculation of Pc of the data within the CDM (including needed frame conversions) isn't the whole story since the time of closest approach (TCA) provided within the file is adjusted to the closest millisecond. The actual TCA of events does not conventiently line up with this millisecond boundary; therefore the states provided within a CDM need to be adjusted to match the actual TCA of the event.

Admittedly, the calculation of the Pc from a CDM file can seem a little convoluted. While this repository has provided the functions necessary to calculate a Pc from the CDM file, it did not provide a step-by-step set of instructions on how to properly perform this calculation...until now. Several scripts have been provided as examples for calculating and comparing Pc values:

| Script Name                          | Script Description |
| ------------------------------------ | ------------------ |
| `Pc2D_FromCDM.m`                     | Provides an example of how to read a CDM file, adjust states to TCA, and calculate the Pc value. |
| `PcMultiStep_FromCDM.m`              | Provides an example of how to read a CDM file, and calculate the Pc value. Note that adjusting TCA states are not needed since PcMultiStep already automatically does this. |
| `ComparePc2D_to_CARAValues.m`        | Provides an example of how to compare locally computed 2D-Pc values against CARA 2D-Pc values. In addition, includes Pc values obtained without the TCA adjustment. |
| `ComparePcMultiStep_to_CARAValues.m` | Provides an example of how to compare locally computed Pc values (Pc2D, Nc2D, Nc3D, or SDMCPc). |
| `GeneratePcMethodPcData.m`           | Provides an example of how the Pc values (and some associated data) were computed to generate `CARA_PcMethod_Test_Conjunctions.xlsx`. |

Please review the embedded documentation within the scripts for more details.

# Change History

| Developer | Update Date | Description |
| --------- | ----------- | ----------- |
| L. Baars  | 04/29/2025  | Initial version, added new README.md |
| L. Baars  | 07/03/2025  | Revised some wording in the 'Testing Information' section and added info for Mac Silicon users. Added new section on calculating Pc Values from CDM Files. |
| L. Baars  | 09/04/2025  | Added information for the PcMultiStep release. |