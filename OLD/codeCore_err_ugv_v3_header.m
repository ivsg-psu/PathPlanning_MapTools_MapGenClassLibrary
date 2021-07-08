    Detailed description is as follows:
    - scanning_results is a Nx3 matrix, N being the number of received
    lidar points or sensor scans
        -- entry 1: R (range or distance of laser point)
        -- entry 2: beta (scanning angle of respective laser point)
        -- entry 3: kappa (heading of vehicle wrt Euclidian line)
    
    - position_uncertainty is a vector of constants
        -- entry 1: dx0 (uncertainty of sensor location in x direction)
        -- entry 2: dy0 (uncertainty of sensor location in y direction)
        -- entry 3: dz0 (uncertainty of sensor location in z direction)
    
    -angular_uncertainty is a vector of constants
        -- entry 1: domega (uncertainty of sensor pointing angle around the
        x axis)
        -- entry 1: dpsi (uncertainty of sensor pointing angle around the
        y axis)
        -- entry 1: dkappa (uncertainty of sensor pointing angle around the
        z axis)
    
    -laser_uncertainty is a vector of constants
        -- entry 1: dbeta (uncertainty of scanning angle around x axis)
        -- entry 2: dR (uncertainty of scanning range value)
    
Enter in units of meters and degrees.
    
This function is designed to output a 'bubble' or 'shadow' of uncertainty
that pertains to how well a sensor on a UGV perceives the objects within
its scanning range. The inputs to this function are primarily error values
specified by the sensor itself, as well as the distance and scanning angles
to the object or point (e.g. for lidar, consider one reflection from an
emitted pulse at a given scanning angle). Error values will be defined on a
global (map) coordinate system, not the local (sensor/vehicle) system.
    
- x,y,z : local coordinate system of sensor, beam origin [meters]
    - x : direction of forward travel
    - y : direction of sideways travel (i.e. left, right)
    - z : altitude or distance above ground
- X,Y,Z : global coordinate system of ground or targeted object, nadir [meters]
    - X : intended path of travel (parallel to Euclidian distance line)
    - Y : deviations from intended path occur along Y (for flat ground)
    - Z : global reference to ground for altitude/height
- omega : roll, rotation around x axis [degrees]
- psi : pitch, rotation around y axis [degrees]
- kappa : yaw, rotation around z axis, angle from Euclidian line to actual
    line of travel, should include the heading uncertainty [degrees]
- beta : scan angle from x0,y0,z0 to X0,Y0,Z0 [degrees]
- R : range or distance from perceived object [meters]
    
Assumptions:
    - flat terrain
    - rake scanning across y axis (or y-z plane but not considering
    elevation or angle), path of laser is parallel to y-x plane (flat)
    - assumes that scanning frequency (scanning paths completed per second)
    is much greater than vehicle velocity (essentially, the vehicle moves a
    negligible distance between each return of the laser within a scan)
    - ignoring other error sources: time offset (scanner to clock),
    calibration offset or misalignment between sensors, possible errors in
    the transformation (post-processing)in the local coordinate system, 
    number, distribution, and distance of GPS reference stations, quality
    of the GPS/INS postprocessing, correction of the relative errors
    through block adjustment of the scan strips, etc.

*Total Error is actually magnitude error, calculated by taking the square
root of the sum of the square of individual error values.
