# Matlab Code Information

This directory contains all of the Software Development Kit (SDK) Matlab code. It has been divided into directories based on major pieces of functionality used within the Conjunction Assessment (CA) enterprise. The following sections provide a brief synopsis of the functionality contained herein.

## Collision Consequence

Risk is properly considered as the combination of likelihood and consequence; but conjunction assessment has usually limited itself to the consideration of only collision likelihood. When considered from an orbital regime protection perspective, the focus shifts to the question of the amount of debris that a collision might produce (the “consequence”). An operational algorithm for determining the expected amount of debris production should a conjunction result in a collision, has been proposed and previously validated.

## Covariance Realism

This tool examines residual sets to determine whether their associated covariances reasonably represent the residuals’ expected distribution. The typical application is to examine, for a particular spacecraft, sets of predicted ephemerides that provide both predicted states and predicted position covariances, in the presence of a definitive (as-flown) ephemeris that will allow residuals between the predicted and definitive ephemerides to be calculated.

## Monte Carlo Pc

A method that may be employed as a method of determining the probability of collision is to perform a Monte Carlo simulation of both the primary and secondary object states at the time of close approach and statistically determine the probability of collision based on the number of trials which violate a predetermined proximity threshold. A Monte Carlo simulation is a computational technique that allows a probabilistic process to be modelled using random sampling from a known multivariate distribution. As the process of orbit determination yields both a best estimate of a satellite state and its associated uncertainty, the problem of collision probability lends itself well to using Monte Carlo sampling methodology.

## OD Quality Assessment

The Orbit Determination Quality Assessment (ODQA) algorithm is intended to analyze an input Conjunction Data Message (CDM) and assess whether the controls and inputs used in the orbit determination process are sufficient to have general confidence in output data products, or if these inputs require further review.

## Probability of Collision

The Probability of Collision (Pc) provides a metric to assess the risk that a satellite close approach event leads to an actual collision. There are a number of formulations for calculating Pc provided within this directory.

## Single Covariance Max Pc

Frisbee (2015) proposed a method by which the maximum possible probability of collision could be determined for a close approach event for which only one object has position uncertainty information. This is of particular use in determining whether an encounter may be of risk to an asset, as the maximum probability of collision may be below an actionable threshold. To determine the maximum probability of collision, the covariance ellipsoid of the object possessing a covariance matrix is mapped to the conjunction plane and distended so that the Mahalanobis distance between the two objects is equal to one. To do this, the covariance of the secondary object is oriented along a one dimensional position uncertainty along the miss vector between the two objects.

## Utils

General purpose utilities that are used across the CARA SDK baseline.