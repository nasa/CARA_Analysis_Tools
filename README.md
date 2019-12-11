# CARA_Analysis_Tools
The tools provided here are free and open-source.

**Conjunction Consequence Assessment:**
> Risk is properly considered as the combination of likelihood and consequence; but conjunction assessment has usually limited itself to the consideration of only collision likelihood. When considered from an orbital regime protection perspective, the focus shifts to the question of the amount of debris that a collision might produce (the “consequence”). An operational algorithm for determining the expected amount of debris production should a conjunction result in a collision, has been proposed and previously validated.

**Monte Carlo from TCA:**
> A method that may be employed as a method of determining the probability of collision is to perform a Monte Carlo simulation of both the primary and secondary object states at the time of close approach and statistically determine the probability of collision based on the number of trials which violate a predetermined proximity threshold.
A Monte Carlo simulation is a computational technique that allows a probabilistic process to be modelled using random sampling from a known multivariate distribution. As the process of orbit determination yields both a best estimate of a satellite state and its associated uncertainty, the problem of collision probability lends itself well to using Monte Carlo sampling methodology.

**Single Covariance Max Pc:**
> Frisbee (2015) proposed a method by which the maximum possible probability of collision could be determined for a close approach event for which only one object has position uncertainty information. This is of particular use in determining whether an encounter may be of risk to an asset, as the maximum probability of collision may be below an actionable threshold. To determine the maximum probability of collision, the covariance ellipsoid of the object possessing a covariance matrix is mapped to the conjunction plane and distended so that the Mahalanobis distance between the two objects is equal to one. To do this, the covariance of the secondary object is oriented along a one dimensional position uncertainty along the miss vector between the two objects.

**Two-Dimension Pc:**
> One method that may be employed as a method of determining the probability of collision is by transforming a close approach event from a three dimensional problem to a two dimensional problem which greatly simplifies the calculation of the probability of collision. This calculation is widely used to characterize and analyze close approach events and determine resultant probabilities of collision as a result of mitigation actions.
