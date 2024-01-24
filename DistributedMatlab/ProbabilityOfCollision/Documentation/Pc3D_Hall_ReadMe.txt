
The Matlab function Pc3D_Hall.m calculates satellite conjunction 
probability of collision (Pc) estimates using the 3D-Nc method,
as formulated by

	D.Hall, "Expected Collision Rates for Tracked Satellites"
	Journal of Spacecraft and Rockets, Vol. 58, No. 3,
	pp. 715–728, May–June 2021, https://doi.org/10.2514/1.A34919

For additional information on the 3D-Nc method, see

 D.Hall, L.Baars, S.Casali, "A Multistep Probability of Collision
 Computational Algorithm" AAS 23-398, 2023.
	
Also, see the file PDF file included with this distribution.

Two usage/verification examples for the Pc3D_Hall.m function
are contained in:

	Pc3D_Hall_Example1.m
	Pc3D_Hall_Example2.m

Also, the function

	AlfanoTestCases_Pc3D.m 

shows how the 3D-Nc method of Pc estimation can be applied to the
twelve test cases provided by Sal Alfano in his 2009 paper

	S.Alfano, "Satellite Conjunction Monte Carlo Analysis"
	AAS 09-233, 2009.

The software in this distribution is covered by the

	NASA OPEN SOURCE SOFTWARE AGREEMENT
	
as explained in detail in the file NOSA_GSC-18593-1.pdf file included with
this distribution, and subject to the following copyright

	Copyright (c) 2020 United States Government as represented by the
	Administrator of the National Aeronautics and Space Administration.
	All Rights Reserved.

Doyle T. Hall
Omitron, Inc.
Senior Conjunction Assessment Research Scientist
Supporting the NASA CARA Analysis Team
Last updated: 2023 September 27