:author: Pedro Henrique Hasselmann
:email: hasselmann at on dot br
:institution: Observatorio Nacional, Rio de Janeiro, Brazil

GmodeClass
==========

   The original G-mode was a clustering method developed by A. I. Gavrishin in the 70's for geochemical classification of rocks, 
   but was also applied to asteroid photometry, cosmic rays, lunar sample and planetary science spectroscopy data. 
   In this work, we used an adapted version to classify the asteroid photometry from SDSS Moving Objects Catalog. 
   The method works by identifying normal distributions in a multidimensional space of variables. 
   The identification starts by locating a set of points with smallest mutual distance in the sample, 
   which is a problem when data is not planar. Here we present a modified version of the G-mode algorithm,
   which was previously written in FORTRAN 77, in Python 2.7 and using NumPy, SciPy and Matplotlib packages. 
   The NumPy was used for array and matrix manipulation and Matplotlib for plot control. 
   The Scipy had a import role in speeding up G-mode, ``Scipy.spatial.distance.mahalanobis`` was chosen as distance estimator and 
   ``Numpy.histogramdd`` was applied to find the initial seeds from which clusters are going to evolve. 
   Scipy was also used to quickly produce dendrograms showing the distances among clusters.

   Finally, results for Asteroids Taxonomy and tests for different sample sizes and implementations are presented.

Introduction
------------

The clusters can be identified using the G-mode multivariate clustering method, designed by A. I. Gavrishin and Coradini [Cor76]_. 
The algorithm was originally written in FORTRAN V by Cor77_ to classify geochemical samples [Cor76_, Bia80_], but is also applicable to a wide range of astrophysical fields, 
as Small Solar System Bodies [Bar87_, Bir96_, Ful08_, Per10_], disk-resolved remote sensing [Pos80_, Tos05_, Cor08_, Ley10_, Tos10_], cosmic rays [Gio81]_ and quasars [Cor83]_. 
In 1987, Bar87_ used original G-mode implementation to classify measurements of asteroids made by the Eight-Color Asteroid Survey [Zel85]_ and 
IRAS geometric albedos [Mat86]_ to produce a taxonomic scheme. Using a sample of 442 asteroids with 8 variables, they recognized 18 classes using a confidence level
of 97.7 %. Those classes were grouped to represent the asteroid taxonomic types. G-mode also identified that just 3 variables
were enough to characterize the asteroid taxonomy.  

The G-mode classifies *N* elements into *Nc* unimodal clusters containing *Na* elements each. Elements are described by *M* variables. 
This method is unsupervised, which allows an automatic identification of clusters without any *a priori* knowledge of sample distribution. 
For that, user must control only one critical parameter for the classification, the confidence levels :math:`q_{1}` or
its corresponding critical value :math:`G_{q1}`. Smaller this parameter get, more clusters are resolved and smaller their spreads are.

So, we choose this method to classify the asteroid observations from Sloan Digital Sky Moving Object Catalog, 
the largest data set on photometry containing around 400,000 moving object entries, 
due to its previous success on asteroid taxonomy, unsupervision and lower number of input parameters. 
However, we were aware the computational limitation we were going to face, since the method never was applied to samples larger than 10,000 elements [Ley10]_
and its last implementation was outdated. Therefore, the G-mode used here follows an adapted version of the original method published by Gav92_, 
briefly described by Ful00_ and reviewed by Tos05_ . 
Median central tendency and absolute deviation estimators, a faster initial seed finder and statistical whitening were introduced to produce a more 
robust set of clusters and optimize the processing time. The coding was performed using Python 2.7 with support of Matplotlib, NumPy and SciPy packages [*]_. 
The algorithm can be briefly summarized by two parts: the first one is the cluster recognition and 
the second evaluates each variable in the classification process.


References
----------

.. [Abr72] Abramowitz, M. & Stegun, I. A. 
           *Handbook of Mathematical Functions Handbook of Mathematical Functions*. New York: Dover, 1972.

.. [Ham74] Hampel, F. R. 
           *The Influence Curve and its Role in Robust Estimation*. Journal ofthe American Statistical Association, 1974, 69, 383-393.

.. [Cor76] Coradini, A.; Fulchignoni, M. & Gavrishin, A. I. 
           *Classification of lunar rocks and glasses by a new statistical technique*. The Moon, 1976, 16, 175-190.

.. [Cor77] Coradini, A.; Fulchignoni, M.; Fanucci, O. & Gavrishin, A. I. 
           *A FORTRAN V program for a new classification technique: the G-mode central method*. Computers and Geosciences, 1977, 3, 85-105.

.. [Bia80] Bianchi, R.; Coradini, A.; Butler, J. C. & Gavrishin, A. I. 
           *A classification of lunar rock and glass samples using the G-mode central method*. Moon and Planets, 1980, 22, 305-322.

.. [Pos80] Poscolieri, M. 
           *Statistical reconstruction of a Martian scene - G-mode cluster analysis results from multispectral data population*. 
           Societa Astronomica Italiana, 1980, 51, 309-328.

.. [Gio81] Giovannelli, F.; Coradini, A.; Polimene, M. L. & Lasota, J. P. 
           *Classification of cosmic sources - A statistical approach*. Astronomy and Astrophysics, 1981, 95, 138-142.

.. [Cor83] Coradini, A.; Giovannelli, F. & Polimene, M. L. 
           *A statistical X-ray QSOs classification International*. Cosmic Ray Conference, 1983, 1, 35-38.

.. [Tho84] Tholen, D. J. 
           *Asteroid taxonomy from cluster analysis of Photometry*. Arizona Univ., Tucson., 1984.

.. [Zel85] Zellner, B.; Tholen, D. J. & Tedesco, E. F. 
           *The eight-color asteroid survey - Results for 589 minor planets*. Icarus, 1985, 61, 355-416.

.. [Mat86] Matson, D. L.; Veeder, G. J.; Tedesco, E. F.; Lebofsky, L. A. & Walker, R. G. 
           *IRAS survey of asteroids*. Advances in Space Research, 1986, 6, 47-56.

.. [Bar87] Barucci, M. A.; Capria, M. T.; Coradini, A. & Fulchignoni, M. 
           *Classification of asteroids using G-mode analysis*. Icarus, 1987, 72, 304-324.

.. [Gav92] Gavrishin, A. I.; Coradini, A. & Cerroni, P. 
           *Multivariate classification methods in planetary sciences*. Earth Moon and Planets, 1992, 59, 141-152.

.. [Bir96] Birlan, M.; Barucci, M. A. & Fulchignoni, M. 
           *G-mode analysis of the reflection spectra of 84 asteroids*. Astronomy and Astrophysics, 1996, 305, 984-+.

.. [Fuk96] Fukugita, M.; Ichikawa, T.; Gunn, J. E.; Doi, M.; Shimasaku, K. & Schneider, D. P. 
           *The Sloan Digital Sky Survey Photometric System*. Astrophisical Journal, 1996, 111, 1748-+.

.. [She97] Shevlyakov, G. L. 
           *On robust estimation of a correlation coefficient*. Journal of Mathematical Sciences, Vol. 83, No. 3, 1997.
 
.. [Lup99] Lupton, R. H.; Gunn, J. E. & Szalay, A. S. 
           *A Modified Magnitude System that Produces Well-Behaved Magnitudes, Colors, and Errors Even for Low Signal-to-Noise Ratio Measurements*. 
           Astrphysical Journal, 1999, 118, 1406-1410.
           
.. [Ful00] Fulchignoni, M.; Birlan, M. & Antonietta Barucci, M. 
           *The Extension of the G-Mode Asteroid Taxonomy*. Icarus, 2000, 146, 204-212.

.. [Ive01] Ivezić, v. Z.; Tabachnik, S.; Rafikov, R.; Lupton, R. H.; Quinn, T.; Hammergren, M.; Eyer, L.; Chu, J.; Armstrong, J. C.; Fan, X.; Finlator, K.; 
           Geballe, T. R.; Gunn, J. E.; Hennessy, G. S.; Knapp, G. R.; Leggett, S. K.; Munn, J. A.; Pier, J. R.; Rockosi, C. M.; Schneider, D. P.; 
           Strauss, M. A.; Yanny, B.; Brinkmann, J.; Csabai, I.; Hindsley, R. B.; Kent, S.; Lamb, D. Q.; Margon, B.; McKay, T. A.; Smith, J. A.; Waddel, P.; York, D. G. & the SDSS Collaboration.
           *Solar System Objects Observed in the Sloan Digital Sky Survey Commissioning Data*. Astrophysical Journal, 2001, 122, 2749-278.

.. [Tos05] Tosi, F.; Coradini, A.; Gavrishin, A. I.; Adriani, A.; Capaccioni, F.; Cerroni, P.; Filacchione, G. & Brown, R. H. 
           *G-Mode Classification of Spectroscopic Data*. Earth Moon and Planets, 2005, 96, 165-197.

.. [Cor08] Coradini, A.; Tosi, F.; Gavrishin, A. I.; Capaccioni, F.; Cerroni, P.; Filacchione, G.; Adriani, A.; Brown, R. H.; Bellucci, G.; 
           Formisano, V.; D'Aversa, E.; Lunine, J. I.; Baines, K. H.; Bibring, J.-P.; Buratti, B. J.; Clark, R. N.; Cruikshank, D. P.; Combes, M.; 
           Drossart, P.; Jaumann, R.; Langevin, Y.; Matson, D. L.; McCord, T. B.; Mennella, V.; Nelson, R. M.; Nicholson, P. D.; Sicardy, B.; Sotin, C.; 
           Hedman, M. M.; Hansen, G. B.; Hibbitts, C. A.; Showalter, M.; Griffith, C. & Strazzulla, G. 
           *Identification of spectral units on Phoebe*. Icarus, 2008, 193, 233-251.

.. [Ful08] Fulchignoni, M.; Belskaya, I.; Barucci, M. A.; de Sanctis, M. C. & Doressoundiram, A. Barucci, M. A.,
           *Transneptunian Object Taxonomy*. The Solar System Beyond Neptune, 2008, 181-192.

.. [Per10] Perna, D.; Barucci, M. A.; Fornasier, S.; DeMeo, F. E.; Alvarez-Candal, A.; Merlin, F.; Dotto, E.; Doressoundiram, A. & de Bergh, C. 
           *Colors and taxonomy of Centaurs and trans-Neptunian objects*. Astronomy and Astrophysics, 2010, 510, A53+.

.. [Ive10] Ivezic, Z.; Juric, M.; Lupton, R. H.; Tabachnik, S.; Quinn, T. & Collaboration, T. S. 
           *SDSS Moving Object Catalog V3.0*. 
           NASA Planetary Data System, 2010, 124.

.. [Ley10] Leyrat, C.; Fornasier, S.; Barucci, A.; Magrin, S.; Lazzarin, M.; Fulchignoni, M.; Jorda, L.; Belskaya, I.; Marchi, S.; Barbieri, C.; Keller, U.; Sierks, H. & Hviid, S. 
           *Search for Steins surface inhomogeneities from OSIRIS Rosetta images*. 
           Planetary and Space Science, 2010, 58, 1097-1106.

.. [Tos10] Tosi, F.; Turrini, D.; Coradini, A. & Filacchione, G. 
           *Probing the origin of the dark material on Iapetus*. Monthly Notices of the Royal Astronomical Society, 2010, 403, 1113-1130.
           
.. [Car10] Carvano, J. M.; Hasselmann, P. H.; Lazzaro, D. & Mothé-Diniz, T. 
           *SDSS-based taxonomic classification and orbital distribution of main belt asteroids*. 
           Astronomy and Astrophysics, 2010, 510, A43+.

