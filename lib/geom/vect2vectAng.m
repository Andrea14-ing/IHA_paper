function [ angle ] = vect2vectAng( V1,V2)
% Angle between two vectors
% Output in deg. range [0, 180]

[m1, v1] = mvarray(V1);
[m2, v2] = mvarray(V2);
dt = dot(v1,v2,2);
angle = acos(dt);
angle = angle.*180./pi;
end
