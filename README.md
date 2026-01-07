# CARA Analysis Tools

This Software Development Kit (SDK) provides a set of tools and utilities that have been developed by NASA's Conjunction Assessment Risk Analysis (CARA) team. These tools and data products are designed to help users assess the risk of close approach events between Earth-orbiting satellites. This software is provided to the public under several NASA Open Source Software Agreements (see the `NOSA_GSC-18593-1.pdf`, `NOSA_GSC-18848-1.pdf`, and `NOSA_GSC-19374-1.pdf` files to view the agreements).

## Git LFS
With the January 2026 release all `*.mat` files within the baseline are controlled with the Git Large File Storage (`lfs`) extension. This extension enables the storage of large binary files within a git repository in an efficient manner. *Git for Windows* includes `lfs` automatically, but other git installations will require a separate install. Please visit https://git-lfs.com/ for download and install instructions.

> Note: Without `lfs` installed any `*.mat` files in the repository will download with text to links to the needed files. `lfs` is required in order to get the actual files within the repository.

After downloading and installing `lfs`, run the following commands in a terminal (*Git Bash* on Windows, or a regular terminal on other operating systems) to activate `lfs` and make sure all `*.mat` files are actually downloaded:

```bash
git lfs install
git lfs pull
```

> Note: These commands are only needed one time to install the `lfs` utility and then to fix any files that have been downloaded as links. After this, any updates to `lfs` controlled files within the repository will automatically download in the correct format.

## Directory/File Structure

The following directories and files will be found at the top level of the SDK:

- `DataFiles` - Contains common static data that are used by one or more tools in the SDK baseline. This data may be used by unit tests or by the tools themselves.
- `DistributedMatlab` - Contains Matlab code and tools that are available for public release. See the `README.md` file within the `DistributedMatlab` directory for a more comprehensive description of the tools provided.
- `References` - Contains reference documents which provide the theoretical basis for many of the tools provided in the SDK.
- `NOSA_GSC-18593-1.pdf`, `NOSA_GSC-18848-1.pdf`, and `NOSA_GSC-19374-1.pdf` - NASA Open Source License Agreements.
- `README.md` - Help document describing the SDK.

## Help Documentation

Most top-level functions should provide help text. From the Matlab prompt, simply use the `help` or `doc` commands to view the help text. Help documentation can generally be found within a directory's `README.md` file or within a `Documentation` subdirectory. If a piece of code or documentation lists a set of references, then the referenced document should appear within the `References` directory.

## Unit Tests

A number of the functions provided in the SDK also come with a set of automated unit tests. After downloading the SDK, it might be a good idea to run all of the unit tests using the following Matlab command within the SDK folder:

```matlab
results = runtests('IncludingSubfolders',true);
```

Matlab will then find and run all unit tests included within the SDK. This will ensure that you have all of the necessary toolboxes needed to run the CARA code. 

If any unit tests fail, a table of failed tests will be displayed to the screen at the end of the run. Running the `disp(results)` command will show an overview of the total tests passed, failed, or had incomplete runs as well as an overall runtime. Running the `table(results)` command will show similar information for each test case that was run.