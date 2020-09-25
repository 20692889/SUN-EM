function [x,y,z] = transformSpericalCoordinateToCartesian(r,theta_radians,phi_radians)
    %transformSpericalCoordinateToCartesian
    %   Usage:
    %           [x,y,z] = transformSpericalCoordinateToCartesian(r,theta,phi)
    %
    %   Input Arguments:
    %       r, theta_radians, phi_radians
    %           The spherical co-ordinate. Note, angle is in radians
    %
    %   Output Arguments:
    %       x,y,z
    %           The cartesian point
    %
    %   Description:
    %       Transforms Spherical Coordinates to Cartesian
    %
    %   =======================
    %   Written by Danie Ludick on 2018.05.04
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za
    %   =======================

    x = r*sin(phi_radians)*cos(theta_radians);
    y = r*sin(phi_radians)*sin(theta_radians);
    z = r*cos(phi_radians);
    
end