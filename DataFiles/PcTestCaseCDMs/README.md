# Testing Information
The `CARA_PcMethod_Test_Conjunctions.xlsx` file contains the results of running various Pc Methods (Pc2D, Nc2D, Nc3D, and PcSDMC) against the CDMs provided in this directory. The methods map to actual Matlab functions as follows:
| Method | Matlab Function | Released |
| ------ | --------------- | -------- |
| Pc2D   | PcCircle.m      | Yes      |
| Nc2D   | Pc2D_Hall.m     | No       |
| Nc3D   | Pc2D_Hall.m     | Yes      |
| PcSDMC | Pc_SDMC.m       | No       |

Some of the methods have been released for public release and are available in the `DistributedMatlab/ProbabilityOfCollision` folder. The other methods are going through the NASA public release process and will be released at a future date.

We have found that Matlab outputs can differ based on the Matlab version being used, the AVX Math Kernel Library (MKL) flags used, and the Operating System used. It is important to be aware of all of these parameters when comparing values.

Newer versions of Matlab can and will use updated versions of the BLAS and LAPACK low-level Fortran libraries. These updates generally represent algorithm enhancements for improved accuracy and/or performance. However, these updates can also prove to be problematic for numerical validation of results, and should be taken into consideration when comparing against the numbers provided in the spreadsheet. Different Matlab versions can account for Pc relative tolerance differences near 1e-8 or below.

The AVX Math Kernel Library can be controlled by setting the MKL_ENABLE_INSTRUCTIONS environment variable and restarting Matlab. Similarly to the Matlab versions, different MKL settings can account for Pc relative tolerance differences near 1e-8 or below.

We've also encountered very slight computation differences depending on the operating system used. When testing across Windows and Linux, we've seen Pc relative tolerance differences near 1e-15 or below.

The BLAS and LAPACK build info parameters are generally reflections of the Matlab version and MKL settings. In the Matlab prompt, the `version -blas` and `version -lapack` commands will show the settings on your machine.

The results reported within `CARA_PcMethod_Test_Conjunctions.xlsx` were generated using the following configuration:
| Item Description                | Value                 |
| ------------------------------- | --------------------- |
| Last Update                     | 4/25/2025             |
| Matlab Version                  | R2019b                |
| MKL_ENABLE_INSTRUCTIONS setting | AVX2                  |
| Operating System                | Windows 11 Enterprise |
| BLAS Build Info                 | Intel(R) Math Kernel Library Version 2018.0.3 Product Build 20180406 for Intel(R) 64 architecture applications, CNR branch AVX2 |
| LAPACK Build Info               | Intel(R) Math Kernel Library Version 2018.0.3 Product Build 20180406 for Intel(R) 64 architecture applications, CNR branch AVX2 Linear Algebra PACKage Version 3.7.0 |

Be advised that differences in the settings listed above will likely cause Pc calculation differences on your system.