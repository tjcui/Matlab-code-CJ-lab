function ODE2D_Quiver_Neuron(low1,step1,high1,low2,step2,high2,w11,w12,w21,w22)
% function ODE2D_Quiver_Neuron(low1,step1,high1,low2,step2,high2,w11,w12,w21,w22)
%
% produce a direction field for the system
%
% dv1/dt = -5*v1(t) + w11*atan(v1(t)) + w12*atan(v2(t))
% dv2/dt = -1*v2(t) + w21*atan(v1(t)) + w22*atan(v2(t))
%
% on the rectangle [low1,high1] x [low2,high2] with grid spacing "step"
% Define the grid points:
[v1,v2] = meshgrid([low1:step1:high1],[low2:step2:high2]);
% Compute the derivatives at each grid point:
dv1 = -5*v1 + w11*atan(v1) + w12*atan(v2);
dv2 = -1*v2 + w21*atan(v1) + w22*atan(v2);
% Put the (dv1,dv2) arrow at each grid point (v1,v2)
% and scale (time units) to make the arrows pretty:
quiver(v1,v2,dv1,dv2);
return
1