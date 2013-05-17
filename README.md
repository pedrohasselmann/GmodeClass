GmodeClass
==========

G-mode Classification method (Gavrishin et al. 1992) for Python 2.7.2 using Numpy 1.5.1, Scipy 0.9 and Matplotlib 1.0.1 .

Introduction
------------

Asteroid taxonomy is a important tool for understanding the spectral diversity and to map
the compositional distribution in the Main Belt. Currently, there is a great number of databases
with many asteroids allowing us to collectively assess their physical properties. Thus, our goal is
to provide a new taxonomy for large asteroid samples like the Sloan Digital Sky Survey Moving
Object Catalog. We expect to verify the strength of the current asteroid taxonomies and ﬁnd
new spectral patterns that may give origin to new asteroid types.

To do so, we decided to use the G-mode statistical classiﬁcation method, which, without any
prior information, detects groups in a space of relevant sample variables. The G-mode chooses
a initial seed within the sample to calculate the G-parameter for each element and through a null
hypothesis test selects the class members. This procedure is repeated over multiple runs until
the class estimators are stabilized. The method then moves over to the identiﬁcation of the next
class, and so on. At the end of the procedure, the G-mode calculates the euclidean distances
between the groups and veriﬁes the relevant variables for the class recognition. If there is any
insigniﬁcant variable, the method is reloaded.

The identification starts with finding the closest points in the sample, 
which is a consumable problem when the data is not planar. 
Therefore, to achieve satisfying results in a human bearable time, we implemented the method, 
which was previously written in FORTRAN 77, in PYTHON 2.7.2 and using Numpy, Scipy and Matplotlib packages. 
The Numpy was used for array and matrix manipulation and Matplotlib for plot control. 
The Scipy had a import role in speeding up G-mode,  Scipy.cKD-Tree and Numpy.histogramdd were applied to 
find the initial seeds from which the clusters are going to evolve. 
Scipy was also used to quickly produce dendograms showing the distances between the clusters.
