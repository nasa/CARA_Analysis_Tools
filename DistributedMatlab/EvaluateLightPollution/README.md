# Evaluate Light Pollution Tool
The Evaluate Light Pollution tool assesses the risk that a satellite constellation will increase light pollution to levels detrimental to ground-based astronomical observations. The following sections provide an overview of the key scripts and directories in a logical order to help you get started.

## Documentation
The 'Documentation' folder is a good starting point. It contains a slide package with an overview of the analysis theory, description of the tool's functionality, and instructions on how to run it.

## References
To gain a deeper understanding of the research underlying the light pollution evaluation, read the 'Hall_2023_ConstellationLightPollutionEvaluation_JAS' paper in the 'References' folder one directory level up.

## EvaluateLightPollution.m
After reviewing the research paper and software documentation, open 'EvaluateLightPollution.m', the core function that runs the analysis. The help documentation at the top of the script provides an overview of the tool, including a detailed description of inputs and outputs.

## RunExamples.m
To see the code execution, check out 'RunExamples.m'. The script includes multiple example scenarios for predicting constellation light pollution. Each example has a significant runtime, so consider running only a few that interest you to start. These examples showcase different configuration options and typical outputs. 

## params
The 'params' directory contains input parameter files for 'EvaluateLightPollution.m'. Reviewing these input examples will give you a better understanding of the different capabilities of the tool.

## output_with_distribution
The 'output_with_distribution' directory contains sample outputs generated from example scenarios, showing the types of results the tool produces.

## src
The 'src' directory houses most of the functions used by 'EvaluateLightPollution.m'. The remaining functions are shared with other tools and reside in the 'Utils' folder one directory level up from 'EvaluateLightPollution'.

## unit_test
The unit test is used to verify that outputs remain accurate after modifications to the code are made.