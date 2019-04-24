SIRFS release 1.0

This is the code for "Color Constancy, Intrinsic Images, and Shape Estimation", Jonathan T. Barron and Jitendra Malik, ECCV 2012.

http://www.cs.berkeley.edu/~barron/
http://www.cs.berkeley.edu/~barron/BarronMalikECCV2012.pdf
http://www.cs.berkeley.edu/~barron/BarronMalikECCV2012_supplementary.pdf
http://www.cs.berkeley.edu/~barron/BarronMalikECCV2012.bib

This release is all undocumented and uncommented research code. It might crash, or destroy your computer, or become sentient. Email me with questions and comments: barron@eecs.berkeley.edu, I'll help you out.

This code includes:
 - A heavily modified version of Eero Simoncelli's matlab code for pyramid image processing: http://www.cns.nyu.edu/~lcv/software.html
 - A slightly modified version of Mark Schmidt's optimization toolbox: http://www.di.ens.fr/~mschmidt/Software/minFunc.html
 - A heavily modified version of Carl Edward Rasmussen's "checkgrad" function for checking the analytical gradient of a function.
 - Tudor Dima's morphological image processing function.
 - An extended version of the MIT Intrinsic Image's dataset: http://people.csail.mit.edu/rgrosse/intrinsic/, to which we've added ground-truth shape and illumination data, as described in: http://www.cs.berkeley.edu/~barron/BarronMalikCVPR2012.pdf.
 - Our own "natural illumination" version of the MIT dataset, as described in: http://www.cs.berkeley.edu/~barron/BarronMalikECCV2012.pdf
 - Some pre-processed spherical harmonic illuminations of natural environment maps taken from the sIBL Archive: http://www.hdrlabs.com/sibl/archive.html


to use the code:

(Optionally) run go_TrainModel, to rebuild our priors on shape, reflectance, and illumination. Part of this takes a while, as it does lots of cross-validation, and I've attached the .mat files of the models so it's not necessary.

run go_MIT(), which takes as an argument a string that determines which sort of experiment you want to run. Some examples:


Run on the "natural" illumination version of the "raccoon" object, with unknown illumination:
go_MIT('params.EVAL_NAMES = {''Craccoon''}; params.SOLVE_LIGHT = 1;');

Run on the standard MIT dataset version of the "raccoon" object, using the known illumination:
go_MIT('params.EVAL_NAMES = {''Lraccoon''}; params.SOLVE_LIGHT = 0;');

Run on the entire training set of the standard MIT dataset using a shape-from-contour algorithm:
go_MIT('params.EVAL_NAMES = MIT_TEST_EXPAND({''L''}); params.SOLVE_LIGHT = 0; params.EXPERIMENT = 151;');

Many, many parameter settings can be seen in CONSTANTS.m. Most of these parameters settings are deprecated (again, this is research code) and most of the rest shouldn't be touched. The parameters in "params.MULT_OPTS.saifs" control the weights of each prior, and can be fiddled with (but are set to be optimal on the training set). These are pretty much the only "parameters" of the model.

This code is written around the MIT dataset. To get the code to run on arbitrary images, you'll need to gut go_MIT and pull out the important parts, and possibly gut parts of other functions that reference the ground-truth data.