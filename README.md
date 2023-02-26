
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
                <li><a href="#polytope-field-feneration-functions">Polytope Field Generation Functions</li>
                <ul>
                    <li><a href="#fcn_MapGen_tilePoints">fcn_MapGen_tilePoints - a method to tile points</li>
                    <li><a href="#dependencies">Dependencies</li>
                </ul>            </ul>
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

1.  Make sure to run MATLAB 2020b or higher. Why? The "digitspattern" command used in the DebugTools utilities, one that is used for argument checking in many functions, was released late 2020. And this command is used heavily in the Debug routines. If debugging is shut off, then earlier MATLAB versions will likely work, and this has been tested back to 2018 releases.

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

* [Errata_Tutorials_DebugTools](https://github.com/ivsg-psu/Errata_Tutorials_DebugTools) - The DebugTools repo is used for the initial automated folder setup, and for input checking and general debugging calls within subfunctions. The repo can be found at: https://github.com/ivsg-psu/Errata_Tutorials_DebugTools



    <!--TODO uncomment this when dependency management is added to this repo: Each should be installed in a folder called "Utilities" under the root folder, namely `./Utilities/DebugTools/`. If you wish to put this code in different directories, the main call stack in `script_demo_MapGenLibrary` can be easily modified with strings specifying the different location, but the user will have to make these edits directly. -->
    
    <!--TODO uncomment this when dependency management is added to this repo: For ease of getting started, the zip files of the directories used - without the .git repo information, to keep them small - are included in this repo.-->


<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>

<!-- FUNCTION DEFINITIONS -->
## Functions

### Polytope Field Generation Functions

#### fcn_MapGen_tilePoints 
The function  **fcn_MapGen_tilePoints** takes as an input set of Nx2 vector of points that specify a points in an area "X". The function returns a tiling of the points by repeating them but with coordinate displacements along the Axis-aligned Bounding Box (AABB), a certain tile "depth". For example, if a region "X" is specified and a tiling depth of 1 is used, this returns tiling points that make a 1-unit boundary around X, as:


    Y Y Y
    Y X Y
    Y Y Y

For a tiling depth of 2, then a boundary of 2 repetitions are formed around the region "X" as:

    Y Y Y Y Y
    Y Y Y Y Y
    Y Y X Y Y
    Y Y Y Y Y
    Y Y Y Y Y

Note: a tile depth of 0 returns simply X, the center region.

The points are ordered such that the resulting matrix is Kx2, with K = (N*(2d+1)^2), and where d is the depth. Thus a depth of 2, which produces a 5x5 tiling, will produce N*25x2 matrix. The original points are in the middle-most portion of the matrix. Specifically, the resulting tile sets, S1 to SK, are organized as follows, using the depth=1 case as an example:

    S3 S6 S9
    S2 S5 S8
    S1 S4 S7

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


**Basic Support Functions**

 `fcn_MapGen_increasePolytopeVertexCount.m` : Some path planners may be restricted to navigating to points that are polytope vertices.  To provide a field that gives these planners higher resolution paths, it may be desirable to have colinear vertices on polytope sides.  This function accomplishes that.
 
 ![image](https://user-images.githubusercontent.com/82562987/220429815-071f8a28-5fe6-4926-951d-3dd930c65a31.png)

`fcn_MapGen_polytopesPredictLengthCostRatio.m` : This is an example of a map analysis tool.  Given a field of polytopes, without planning a path, this will approximate the typical distance cost-to-go for that field.  As expected, this increases with denser fields and larger obstacles.

![image](https://user-images.githubusercontent.com/82562987/220430629-a30a9bc6-6a53-4e06-ba24-abe5e245a896.png)

`fcn_MapGen_polytopesPredictLengthCostRatioStraightPath.m` : This performs a similar analysis to what was described above, but rather than approxiting the behavior of a path planner that routes around obstacles, this approximates the cost resulting path planner that routes through all obstacles in a straight line, traversing all encountered obstacles.

![image](https://user-images.githubusercontent.com/82562987/220430650-5c13292e-096d-4939-af95-be387e37791e.png)

`fcn_MapGen_polytopesPredictUnoccupancyRatio.m` : This function, given a polytope field, both before and after shrinking, will use different methods to predict linear and area unoccupancy, i.e. the amount of free length or area relative to the amount of total length or area.

![image](https://user-images.githubusercontent.com/82562987/220431115-425e5ebe-efeb-4a6f-88a7-08246f4b51a5.png)

`fcn_MapGen_polytopesRadiusDistributions.m` : This function rotates polytopes about their centroids to produce a radial probability of occupation and therefore an effective polytope width, independent of position and orientation.

![image](https://user-images.githubusercontent.com/82562987/220431404-42fd9f6d-db9e-44d7-b3aa-7295b846ad3f.png)

![image](https://user-images.githubusercontent.com/82562987/220431555-5a21073c-3911-4546-a1cf-5be8888a681d.png)

`script_ui_manuallyDefineMapLayers.m` : This script starts a command line and GUI utility that allows the user to load an image of their choosing, often a satellite map, and draw potentially overlapping polytopes on this image in multiple layers.  The layers will then be flattened into a non-overlapping, convex polytope struct array.

Example drawn polytopes:

![image](https://user-images.githubusercontent.com/82562987/220432752-9c57c185-efa2-4538-85be-dd7ab84aca57.png)


Example path planned through this:

![image](https://user-images.githubusercontent.com/82562987/220432801-c217fd61-cc9e-46ba-9b3c-3a8a51197eee.png)


`fcn_MapGen_flattenPolytopeMap.m` : This takes a potentially overlapping set of polytopes and breaks them into non-overlapping triangles to enforce convexity.  Note, polytope traversal costs will be linearlly combined.
Example overlapping polytopes:

![image](https://user-images.githubusercontent.com/82562987/220432234-47b632c4-4c71-4686-86b6-225d30961457.png)

Example non-overlapping but potentially convex polytopes:

![image](https://user-images.githubusercontent.com/82562987/220432274-1a65b2b9-2646-4c4e-9fe0-818f44a4b744.png)


Example non-convex and non-overlapping triangular polytopes:

![image](https://user-images.githubusercontent.com/82562987/220432320-4d819071-6962-421e-807d-e0005cd5b356.png)


<a href="#pathplanning_maptools_mapgenclasslibrary">Back to top</a>

**Core Functions**

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
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.svg?style=for-the-badge
[contributors-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.svg?style=for-the-badge
[forks-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/network/members
[stars-shield]: https://img.shields.io/github/stars/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.svg?style=for-the-badge
[stars-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/stargazers
[issues-shield]: https://img.shields.io/github/issues/ivsg-psu/reFeatureExtraction_Association_PointToPointAssociationpo.svg?style=for-the-badge
[issues-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/issues
[license-shield]: https://img.shields.io/github/license/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.svg?style=for-the-badge
[license-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/blob/master/LICENSE.txt








