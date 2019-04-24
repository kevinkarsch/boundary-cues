This dataset contains boundary cue annotations for several PASCAL object 
categories and the MIT Intrinsic Image dataset. The annotations and their
fieldnames are listed below. For more details, please see our CVPR'13 paper:

Boundary Cues for 3D Object Shape Recovery 
Karsch, K.; Liao, Z.; Rock, J.; Barron, J.T.; Hoiem, D. 
Computer Vision and Pattern Recognition, 2013.

Each folder contains  pairs of images and MATLAB .mat files. The .mat files
contain annotations for the corresponding images, and have the following
fields:

im           - Image data (same as the corresponding jpg image)
extbd        - Exterior boundary (silhoutte) contours
intbd        - Interior boundary (self occlusion) contours
contact_pts* - Locations where the object contacts the ground
horizonline* - Line segment that lies along the horizon [x1 y1 x2 y2]
paraline*    - Sets of parallel line segments in 3D (vanishing lines)
convexcnr    - Surface fold/crease contour (convex wrt camera)
concavecnr   - Surface fold/crease contour (concave wrt camera)

All contours are stored as cell arrays (one contour per array element). Each
contour is a list of vertices are stored as an Nx2 set of [x,y] locations in 
the image. We follow MATLAB's image coordinate system (origin is top left).

* These annotations were only partially collected and not used in the paper.


Also included is our labelling tool. To label your own images, see imlabel.m.
This function can also be used for visualizing existing annoations (examples 
below).

labeldata = imlabel(img_file); %Opens annotator (see imlabel for instructions)
imlabel(img_file, labeldata); %Displays annotations in labeldata on img_file

