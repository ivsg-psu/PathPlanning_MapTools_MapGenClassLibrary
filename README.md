# PathPlanning_MapTools_MapGenClassLibrary

The class library for functions to generate new maps, started from Seth Tau's original work.
<!-- PROJECT LOGO -->
<br>
<p align="center">
  <h2 align="center">
    PathPlanning_MapTools_MapGenClassLibrary
  </h2>
  <pre align="center">
        <img src=".\Images\MapGenClass.jpg" alt="main mapgen picture" width="960" height="540">
        <!-- figcaption>Fig.1 - The typical progression of map generation.</figcaption -->
        <font size="-2">Photo by <a href="https://unsplash.com/@marjan_blan?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Marjan Blan | @marjanblan</a> on <a href="https://unsplash.com/photos/6bXvYyAYVrE?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a>
    </font>
  </pre>
</p>

<p align="left">
    The MapGen class library hosts tools to build maps useful for studying path planning, autonomy, etc. The purpose of this code is to generate, modify, and analyze maps of obstacle fields.  The obstacles are stored as polytopes in vertex representation.  For information on this procedure and its utilities please see <a href="https://etda.libraries.psu.edu/catalog/18657sat5340">Seth Tau's thesis</a> or <a href="http://gvsets.ndia-mich.org/documents/AAIR/2022/Path-Free%20Estimation%20of%20Navigation%20Distance%20using%20Obstacle%20Shape%20Statistics%20and%20Density.pdf"> this paper </a>.
    <br>
    <a href="https://github.com/ivsg-psu/ivsg_master/wiki/Path-Planning#map-tools"><strong>Explore the docs</strong></a>
    <a href="https://github.com/ivsg-psu/PathPlanning_MapTools_MapGenClassLibrary">View Demo</a>
    <a href="https://github.com/ivsg-psu/PathPlanning_MapTools_MapGenClassLibrary/issues">Report Bug</a>
    <a href="https://github.com/ivsg-psu/PathPlanning_MapTools_MapGenClassLibrary/issues">Request Feature</a>
</p>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
    <ol>
        <li>
            <a href="#about-the-project">About the Project</a>
        </li>
        <li>
        <a href="#getting-started">Getting Started</a>
        <ul>
                <li><a href="#installation">Installation</a></li>
        </ul>
        </li>
        <li><a href="structure">Repo Structure</a>
            <ul>
                <li><a href="#directories">Top-Level Directories</li>
                <li><a href="#dependencies">Dependencies</li>
            </ul>
        </li>
        <li><a href="#functions">Functions</li>
            <ul>
                <li><a href="#single-polytope-functions">Single Polytope Functions</li>
                <ul>
                    <li><a href="#fcn_mapgen_iswithinaabb">fcn_MapGen_isWithinABBB - a method to test if points are within an AABB</li>
                    <li><a href="#fcn_mapgen_plotpolytopes">fcn_MapGen_plotPolytopes - a method to plot polytope arrays</li>
                    <li><a href="#fcn_mapgen_polytopecentroidandarea">fcn_MapGen_polytopeCentroidAndArea - calculates polytope centroid and area from vertices</li>
                    <li><a href="#fcn_mapgen_fillpolytopefieldsfromvertices">fcn_MapGen_fillPolytopeFieldsFromVertices - calculates polytope properties from vertices</li>
                </ul>
                <li><a href="#polytope-field-generation-functions">Polytope Field Generation Functions</li>
                <ul>
                    <li><a href="#fcn_mapgen_generatepolysfromvoronoiaabb">fcn_MapGen_generatePolysFromVoronoiAABB - a method to generate polytope fields by "manaully" bounding the Voronoi diagram</li>
                    <li><a href="#fcn_mapgen_tilepoints">fcn_MapGen_tilePoints - a method to tile points</li>
                </ul>
            </ul>
        <li><a href="#usage">Usage</a></li>
            <ul>
                <li><a href="#examples">Examples</li>
                <li><a href="#definition-of-endpoints">Definition of Endpoints</li>
            </ul>
        <li><a href="#license">License</a></li>
        <li><a href="#contact">Contact</a></li>
    </ol>
</details>

***
<!-- ABOUT THE PROJECT -->
## About The Project

The MapGen class library hosts tools to build maps useful for studying path planning, autonomy, etc. The result of these operations is usually a map representation of the expected obstacle environment, and for nearly all the mapping work in the team, maps are stored as polytopes.

Thus, a challenge when testing navigration algorithms is to generate polytope maps either from random inputs or from real data. Most of this library is coded to support random procedurally-generated map creation. As shown in Fig. 1 below, the process typically begins with user-selected seed points. These are then used to generate the Voronoi boundaries - the lines that are equidistant from other points. Next, an axis-aligned bounding box (AABB) is applied that crops the Voronoi boundaries to a specific working area. After this, the statistics for polytopes can be found (such as the mean of all verticies). These polytopes can then be used a the starter polytopes for shrinkage operations including shrinking polytopes to a particular radius, or shrinking equally from each polytope edge. The use of a procedural creation method to create obstacle maps allows the controlled re-creation of these maps across algorithms, and thereby permits analysis of map properties to discover relationships between vehicle mobility and obstacle characteristics.

The most common usage of the maps from this repository is to generate maps for use in testing path planners and path following controllers.  This repository also includes tools for modifying generated maps such as changing obstalce traversal costs, obstacle size, and obstacle size distribution.  Finally, some obstacle field analysis tools are included as well.

* Inputs:
  * desired number of obstalces
  * obstacle size distribution
* Outputs
  * Polytope struct array representing an obstacle field that can be plotted or passed to a path planner

<pre align="center">
    <img src=".\Images\MapGenProgression.png " alt="typical progression of map generation" width="960" height="450">
    <figcaption>Fig.1 - The typical progression of map generation.</figcaption>
    <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>

***
<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Installation

1. Make sure to run MATLAB 2020b or higher. Why? The "digitspattern" command used in the DebugTools utilities, one that is used for argument checking in many functions, was released late 2020. And this command is used heavily in the Debug routines. If debugging is shut off, then earlier MATLAB versions will likely work, and this has been tested back to 2018 releases.

2. Clone the repo

   ```sh
   git clone https://github.com/ivsg-psu/PathPlanning_MapTools_MapGenClassLibrary
   ```

3. Run the main code in the root of the folder (`script_demo_MapGenLibrary.m`). This will download the required utilities for this code, unzip the zip files into a Utilities folder (.\Utilities), and update the MATLAB path to include the Utility locations. This install process will only occur the first time. Note: to force the install to occur again, delete the Utilities directory and clear all global variables in MATLAB (type: "clear global *").
4. Confirm it works! Run `script_demo_MapGenLibrary`. If the code works, the script should run without errors. This script produces numerous example images such as those in this README file.

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>

<!-- STRUCTURE OF THE REPO -->
***

## Structure

### Directories

The following are the top level directories within the repository:
<ul>
    <li>/Documents folder: Descriptions of the functionality and usage of the various MATLAB functions and scripts in the repository.</li>
    <li>/Functions folder: The majority of the code functionalities are implemented in this directory. All functions as well as test scripts are provided.</li>
    <!--TODO uncomment this when dependency management is added to this repo: <li>/Utilities folder: Dependencies that are utilized but not implemented in this repository are placed in the Utilities directory. These can be single files but are most often folders containing other cloned repositories.</li>-->
  <li>/testFixtures folder: This includes .mat workspace files that some test scripts use as inputs so that they can opperate on a known workspace without having to run numerous preprocessing steps and other functions.  This saves time and prevents test flaking.</li>
</ul>

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>

### Dependencies

* [Errata_Tutorials_DebugTools](https://github.com/ivsg-psu/Errata_Tutorials_DebugTools) - The DebugTools repo is used for the initial automated folder setup, and for input checking and general debugging calls within subfunctions. The repo can be found at: [https://github.com/ivsg-psu/Errata_Tutorials_DebugTools](https://github.com/ivsg-psu/Errata_Tutorials_DebugTools)

    <!--TODO uncomment this when dependency management is added to this repo: Each should be installed in a folder called "Utilities" under the root folder, namely `./Utilities/DebugTools/`. If you wish to put this code in different directories, the main call stack in `script_demo_MapGenLibrary` can be easily modified with strings specifying the different location, but the user will have to make these edits directly. -->

    <!--TODO uncomment this when dependency management is added to this repo: For ease of getting started, the zip files of the directories used - without the .git repo information, to keep them small - are included in this repo.-->

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>

<!-- FUNCTION DEFINITIONS -->
## Functions

### Single Polytope Functions

The following are basic functions that operate on single polytopes

#### fcn_MapGen_isWithinAABB

The function  **fcn_MapGen_isWithinAABB** checks if the points are within the given axis-aligned bounding box, AABB, returning a vector of 1' or 0's the same length as the nubmer of rows of points. Each point must be strictly within the AABB - e.g. this function returns "false" if a point is on the "wall" of the AABB.

<pre align="center">
<img src=".\Images\fcn_MapGen_isWithinAABB.png " alt="AABB example" width="420" height="315">
<figcaption>fcn_MapGen_isWithinAABB - applied over the unit square with random points.</figcaption>
<!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>
***

#### fcn_MapGen_plotPolytopes

The function  **fcn_MapGen_plotPolytopes** plots a polytope array. Note: it only requires the verticies subfield to exist in order to work. For example:

```MATLAB
polytopes(1).vertices = [0 0; 4 2; 2 4; 0 0];
fcn_MapGen_plotPolytopes(polytopes,fig_num,'r-',line_width);
```

will give the following figure

<pre align="center">
<img src=".\Images\fcn_MapGen_plotPolytopes.png " alt="plotting example" width="420" height="315">
<figcaption>fcn_MapGen_plotPolytopes - an example of plotting a polytope of three points.</figcaption>
<!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>
***

#### fcn_MapGen_polytopeCentroidAndArea

The function  **fcn_MapGen_polytopeCentroidAndArea** calculates the centroid and area of a polytope, given the closed-form X and y points for the polytope vertices.

For example:

```MATLAB
fig_num = 2; % Define the figure

% Define some points
x = [3; 4; 2; -1; -2; -3; -4; -2; 1; 2; 3];
y = [1; 2; 2; 3; 2; -1; -2; -3; -3; -2; 1];

% Call the function
[Centroid,Area] = fcn_MapGen_polytopeCentroidAndArea([x,y],fig_num);

```

will give the following figure

<pre align="center">
<img src=".\Images\fcn_MapGen_polytopeCentroidAndArea.png " alt="centroid and area calculation example" width="420" height="315">
<figcaption>fcn_MapGen_polytopeCentroidAndArea - an example of calculating centroid and area of a polytope.</figcaption>
<!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>
***

#### fcn_MapGen_fillPolytopeFieldsFromVertices

The function  **fcn_MapGen_fillPolytopeFieldsFromVertices** fills in all fields of a polytope given the vertices.

For example:

```MATLAB
fig_num = 3; % Define the figure number
clear polytopes % Clear the variable
polytopes(1).vertices = [0 0; 4 2; 2 4; 0 0]; % Fill in the verticies
% Call the function to fill in all other fields
polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes,fig_num);
% Show that a new field now exists
assert(isequal(round(polytopes(1).max_radius,4),2.8284));
```

will give a structure with all the fields filled in

```sh
>> polytopes

polytopes =

  struct with fields:

       vertices: [4x2 double]
             xv: [0 4 2]
             yv: [0 2 4]
      distances: [3x1 double]
           mean: [2 2]
           area: 6
     max_radius: 2.8284
     min_radius: 2
    mean_radius: 2.2761
          radii: [3x1 double]
           cost: 0.0688
```

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>
***

### Polytope Field Generation Functions

There are several functions build to generate polytope fields.

#### fcn_MapGen_generatePolysFromVoronoiAABB

The function  **fcn_MapGen_generatePolysFromVoronoiAABB** calculates a polytope field given a set of seed points by calculating the Voronoi diagram and manually bounding it. This is an older method, and one prone to error as there are known edge cases where this fails. Also, it tends to be a slow process. The method's steps are illustrated by the following figure:

<pre align="center">
    <img src=".\Images\MapGenProgression.png " alt="typical progression of map generation" width="960" height="450">
    <figcaption>Fig.1 - The typical progression of map generation.</figcaption>
    <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

The steps are:

1. Given a set of seed points, calculate the Voronoi boundaries.
2. Given the Vornoi boundaries, contain them by imposing an AABB, including calculating new polytopes as needed.
3. With the new polytopes, calculate the statistics.
4. Shrink the polytopes, if necessary, to produce obstacles either by shrinking to a radius, or
5. Shrink the polytopes from their edges

The following code illustrates a basic call:

```Matlab
% Demos fcn_MapGen_generatePolysFromVoronoiAABB
% and script_test_fcn_MapGen_generatePolysFromVoronoiAABB
fig_num = 10;

% pull halton set
halton_points = haltonset(2); % Generate the 2-D Halton set
points_scrambled = scramble(halton_points,'RR2'); % scramble values


% pick values from halton set
Halton_range = [1 100]; % Define number of points we want (e.g 1 to 100)
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
seed_points = points_scrambled(low_pt:high_pt,:);
[V,C] = voronoin(seed_points); % Calculate the Voronoi diagram

AABB = [0 0 1 1]; % Define the axis-aligned bounding box
stretch = [1 1]; % Define the stretch factor on each axis

% fill polytopes from tiling, and give a figure number to show results
fcn_MapGen_generatePolysFromVoronoiAABB(seed_points,V,C,AABB, stretch,fig_num);
```

And it produces the following plot:

<pre align="center">
<img src=".\Images\fcn_MapGen_generatePolysFromVoronoiAABB.png " alt="AABB calculation of polytope fields" width="420" height="315">
<figcaption>fcn_MapGen_generatePolysFromVoronoiAABB - an example of calculating a polytope field constrained by an Axis-Aligned Bounding Box.</figcaption>
<!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>
***

#### fcn_MapGen_tilePoints

The function  **fcn_MapGen_tilePoints** takes as an input set of Nx2 vector of points that specify a points in an area "X". The function returns a tiling of the points by repeating them but with coordinate displacements along the Axis-aligned Bounding Box (AABB), a certain tile "depth". For example, if a region "X" is specified and a tiling depth of 1 is used, this returns tiling points that make a 1-unit boundary around X, as:

```sh
    Y Y Y
    Y X Y
    Y Y Y
```

For a tiling depth of 2, then a boundary of 2 repetitions are formed around the region "X" as:

```sh
    Y Y Y Y Y
    Y Y Y Y Y
    Y Y X Y Y
    Y Y Y Y Y
```

Note: a tile depth of 0 returns simply X, the center region.

The points are ordered such that the resulting matrix is Kx2, with K = (N*(2d+1)^2), and where d is the depth. Thus a depth of 2, which produces a 5x5 tiling, will produce N*25x2 matrix. The original points are in the middle-most portion of the matrix. Specifically, the resulting tile sets, S1 to SK, are organized as follows, using the depth=1 case as an example:

```sh
    S3 S6 S9
    S2 S5 S8
    S1 S4 S7
```

Thus, S1 contains the points on rows (1..N), S2 is ((N+1)...(2N)), S3 is ((2N+1)...(3N)), etc.

Below is an example output illustrating how the center original points are tiled:

<pre align="center">
<img src=".\Images\fcn_MapGen_tilePoints_BasicCall.png " alt="Tiling points example" width="420" height="315">
<figcaption>fcn_MapGen_tilePoints - the basic call.</figcaption>
<!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Note that the tiling does not require the points to be inside the AABB. If, for example, the center points are offset from the AABB, then the points are still tiled, resulting in points that are offset from the tiling pattern:

<pre align="center">
<img src=".\Images\fcn_MapGen_tilePoints_OffsetPoints.png " alt="Tiling points example with offsets" width="420" height="315">
<figcaption>fcn_MapGen_tilePoints - showing point offsets.</figcaption>
<!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Here's an example with a tiling depth of 2, e.g. the number of cells that will "fence" the original points.

<pre align="center">
<img src=".\Images\fcn_MapGen_tilePoints_DepthOf2.png " alt="Tiling points example, depth of 2" width="420" height="315">
<figcaption>fcn_MapGen_tilePoints - a depth of 2.</figcaption>
<!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>

#### fcn_MapGen_generatePolysFromVoronoiAABBWithTiling

The function  **fcn_MapGen_generatePolysFromVoronoiAABBWithTiling** is similar to **fcn_MapGen_generatePolysFromVoronoiAABB** in that it calculates a polytope field given a set of seed points by calculating the Voronoi diagram and manually bounding it. However, it "bounds" the polytope by creating a tiling of seed points around the original set, using the AABB for the tiling, before generating the Voronoi diagram. This prevents edge polytopes from extending to infinity or having close-off issues. Further, the resulting polytope field has a perfect tiling. The method's steps are illustrated by the following figure:

<pre align="center">
    <img src=".\Images\fcn_MapGen_generatePolysFromVoronoiAABBWithTiling.png " alt="steps in creating a Voronio polytope field via tiling" width="960" height="450">
    <figcaption>fcn_MapGen_generatePolysFromVoronoiAABBWithTiling - creates a Voronio polytope field via tiling.</figcaption>
    <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

The steps are:

1. Given a set of seed points, tile these to create a "border" of one copy thick around the original points
2. Calculate the Voronoi boundaries for the tiled points.
3. Extract the middle set and use these to create a polytope field.
4. Shrink the polytopes, if necessary, to produce obstacles either by shrinking to a radius, and
5. Note that the polytopes "tile", e.g. they align along their borders.

The following code illustrates a basic call:

```Matlab
%% pull halton set
halton_points = haltonset(2);
points_scrambled = scramble(halton_points,'RR2'); % scramble values

%% pick values from halton set
Halton_range = [1801 1901];
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
seed_points = points_scrambled(low_pt:high_pt,:);

%% fill polytopes from tiling
fig_num = 1;
AABB = [0 0 1 1];
stretch = [1 1];
polytopes = fcn_MapGen_generatePolysFromVoronoiAABBWithTiling(seed_points,AABB, stretch,fig_num);
```

And it produces the following plot:

<pre align="center">
<img src=".\Images\fcn_MapGen_generatePolysFromVoronoiAABBWithTiling_Example.png " alt="creating a polytope field via tiling" width="420" height="315">
<figcaption>fcn_MapGen_generatePolysFromVoronoiAABBWithTiling - an example of calculating a polytope field constrained by an Axis-Aligned Bounding Box and tiling.</figcaption>
<!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

This mapping method is particularly good for moving obstacles, such as a single polytope:

<pre align="center">
<img src=".\Images\testAnimated_onePoly_25Hz.gif " alt="animated single polytope" width="420" height="315">
<figcaption>A test animation of a single moving polytope, 25 Hz.</figcaption>
<!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

Or many polytopes:

<pre align="center">
<img src=".\Images\testAnimated_allPoly_25Hz.gif " alt="animated many polytopes" width="420" height="315">
<figcaption>A test animation of a many moving polytopes, 25 Hz.</figcaption>
<!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>
***

#### Basic Support Functions


`fcn_MapGen_increasePolytopeVertexCount.m` : Some path planners may be restricted to navigating to points that are polytope vertices.  To provide a field that gives these planners higher resolution paths, it may be desirable to have colinear vertices on polytope sides.  This function accomplishes that.

<img src=".\Images\fcn_MapGen_increasePolytopeVertexCount.png " alt="vertex interpolation example">

`fcn_MapGen_polytopesPredictLengthCostRatio.m` : This is an example of a map analysis tool.  Given a field of polytopes, without planning a path, this will approximate the typical distance cost-to-go for that field.  As expected, this increases with denser fields and larger obstacles.

<img src=".\Images\fcn_MapGen_polytopesPredictLengthCostRatio.png " alt="predicted and measured length cost ratio">

`fcn_MapGen_polytopesPredictLengthCostRatioStraightPath.m` : This performs a similar analysis to what was described above, but rather than approxiting the behavior of a path planner that routes around obstacles, this approximates the cost resulting path planner that routes through all obstacles in a straight line, traversing all encountered obstacles.

<img src=".\Images\fcn_MapGen_polytopesPredictLengthCostRatioStraightPath.png " alt="predicted and measured length cost ratio for straight path">

`fcn_MapGen_polytopesPredictUnoccupancyRatio.m` : This function, given a polytope field, both before and after shrinking, will use different methods to predict linear and area unoccupancy, i.e. the amount of free length or area relative to the amount of total length or area.

<img src=".\Images\fcn_MapGen_polytopesPredictUnoccupancyRatio.png " alt="estimated occupancy ratio">

`fcn_MapGen_polytopesRadiusDistributions.m` : This function rotates polytopes about their centroids to produce a radial probability of occupation and therefore an effective polytope width, independent of position and orientation.

<img src=".\Images\fcn_MapGen_polytopesRadiusDistributions_diagram.png " alt="diagram of statistical polytope representation">

<img src=".\Images\fcn_MapGen_polytopesRadiusDistributions_plot.png " alt="plot of polytope effective depths">

`script_ui_manuallyDefineMapLayers.m` : This script starts a command line and GUI utility that allows the user to load an image of their choosing, often a satellite map, and draw potentially overlapping polytopes on this image in multiple layers.  The layers will then be flattened into a non-overlapping, convex polytope struct array.

Example drawn polytopes:

<img src=".\Images\script_ui_manuallyDefineMapLayers.png " alt="polytopes manually drawn">

Example path planned through this:

<img src=".\Images\script_ui_manuallyDefineMapLayers_path.png " alt="polytopes manually drawn and used for path planning">

`fcn_MapGen_flattenPolytopeMap.m` : This takes a potentially overlapping set of polytopes and breaks them into non-overlapping triangles to enforce convexity.  Note, polytope traversal costs will be linearlly combined.
Example overlapping polytopes:

<img src=".\Images\fcn_MapGen_flattenPolytopeMap.png " alt="overlapping polytopes on a satellite map">

Example non-overlapping but potentially convex polytopes:

<img src=".\Images\fcn_MapGen_flattenPolytopeMap_2.png " alt="non-overlapping flattened polytopes">

Example non-convex and non-overlapping triangular polytopes:

<img src=".\Images\fcn_MapGen_flattenPolytopeMap_3.png " alt="triangulated flattened polytopes">

`fcn_MapGen_polytopesExpandEvenlyForConcave.m`: can be used as an alternative to `fcn_MapGen_polytopesExpandEvenly`.  The existing function works well for convex polytopes but unexpected behavior can occur with concave polytopes.  The new function uses Matlab's polyshape object and polybuffer method to create a buffer around the polytope that is robust to concave input polytopes.

![image](https://github.com/ivsg-psu/PathPlanning_MapTools_MapGenClassLibrary/assets/82562987/501b21a1-d26e-45c0-b6ce-5d08b0839f35)

`fcn_MapGen_calculateConvexHullOverlapRatio.m`: can be used to calculate the area of obstacle convex hulls that overlap relative to the total area covered by obstacles.  This is useful for assessing how much obstacles protrude into concavities of other obstacles.

![image](https://github.com/user-attachments/assets/f5cefffa-9149-472b-8990-cafb9c6ee310)

`fcn_MapGen_generateVoronoiDiagramBetweenPolytopes.m`: implemented as a wrapper for MATLAB's `voronoi` function to call `voronoi` on polytope vertices.  If you're looking to create the medial axis graph, this is not the function to use.  This is useful for comparing the voronoi diagram formed from vertices to the medial axis around the obstacles.  To generate the medial axis between polytopes without these errant edges, you may wish to use the `fcn_MedialAxis_*` functions in [PathPlanning_GridFreePathPlanners_BoundedAStar](https://github.com/ivsg-psu/PathPlanning_GridFreePathPlanners_BoundedAStar)

`fcn_MapGen_isCrossingAABB.m`: checks if line segments intersect the axis-aligned bounding boxes of polytopes.  This is useful for partially modifying the visibility graph by adding or removing edges while only checking for collisions with polytopes whose AABB contacts the graph edge.

`fcn_MapGen_snapInteriorPointToVertex.m`: checks if points are interior to polytopes.  If they are, it overwrites these points with the nearest polytope vertex of the containing polytope.  This is useful when dilating polytopes such that if dilation overlaps a point of interest, that point can be moved to the outside of the polytope.  See [here ](https://github.com/ivsg-psu/PathPlanning_GridFreePathPlanners_BoundedAStar/blob/main/Functions/script_test_polytope_canyon_replan_with_dilation.m) for an example of this in use.

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>

#### Core Functions

`fcn_MapGen_fillPolytopeFieldsFromVertices.m` : This is the core function for this repo that creates a polytope struct from an input set of vertices.  An example of the fields in this struct is shown below:

```matlab
       vertices: [4×2 double]
         xv: [0 4 2]
         yv: [0 2 4]
  distances: [3×1 double]
       mean: [2 2]
       area: 6
 max_radius: 2.8284
 min_radius: 2
mean_radius: 2.2761
      radii: [3×1 double]
       cost: 0.9058
```

</ul>
Each of the functions has an associated test script, using the convention

```sh
script_test_fcn_fcnname
```

where fcnname is the function name as listed above.

As well, each of the functions includes a well-documented header that explains inputs and outputs. These are supported by MATLAB's help style so that one can type:

```sh
help fcn_fcnname
```

for any function to view function details.

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>

<!-- USAGE EXAMPLES -->
## Usage
<!-- Use this space to show useful examples of how a project can be used.
Additional screenshots, code examples and demos work well in this space. You may
also link to more resources. -->

### Examples

1. Run the main script to set up the workspace and demonstrate main outputs, including the figures included here:

   ```sh
   script_demo_MapGenLibrary
   ```

2. After running the main script to define the included directories for utility functions, one can then navigate to the Functions directory and run any of the functions or scripts there as well. All functions for this library are found in the Functions sub-folder, and each has an associated test script. Run any of the various test scripts, such as:

   ```sh
   script_test_fcn_MapGen_plotPolytopes
   ```

For more examples, please refer to the [Documentation](https://github.com/ivsg-psu/PathPlanning_MapTools_MapGenClassLibrary/tree/main/Documents)

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>

## Major release versions

This code is still in development (alpha testing)

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>

<!-- CONTACT -->
## Contact

Sean Brennan - sbrennan@psu.edu

Project Link: [https://github.com/ivsg-psu/PathPlanning_MapTools_MapGenClassLibrary](https://github.com/ivsg-psu/PathPlanning_MapTools_MapGenClassLibrary)

<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>

<!-- MARKDOWN LINKS & IMAGES -->

<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links
[contributors-shield]: https://img.shields.io/github/contributors/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.svg?style=for-the-badge
[contributors-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.svg?style=for-the-badge
[forks-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/network/members
[stars-shield]: https://img.shields.io/github/stars/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.svg?style=for-the-badge
[stars-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/stargazers
[issues-shield]: https://img.shields.io/github/issues/ivsg-psu/reFeatureExtraction_Association_PointToPointAssociationpo.svg?style=for-the-badge
[issues-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/issues
[license-shield]: https://img.shields.io/github/license/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.svg?style=for-the-badge
[license-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/blob/master/LICENSE.txt -->
