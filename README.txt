This MATLAB software accompanies our CVPR'13 paper:

K. Karsch, Z. Liao, J. Rock, J.T. Barron and D. Hoiem. "Boundary Cues for 3D
Object Shape Recovery." CVPR 2013.

Descriptions of parameters and implementation details can be found in the 
paper and supplementary material (available at 
http://kevinkarsch.com). Please cite our paper 
if you use our code for your research.

This folder contains a modified version of Jon Barron's SIRFS code. For the 
original version of Jon's code, please visit his website 
(http://www.cs.berkeley.edu/~barron/).

-----

Sample usage:

img = im2double(imread('dog.png'));
labeldata = imlabel(img); %See imlabel.m for usage info
save('dog.mat', 'labeldata');
MakeStudyData( 'dog.mat' );

This will produce dog-result.{mat,pdf}, which have been included with this 
package for comparison.