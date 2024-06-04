function [max_s] = bladeMatlab(a1,a2,t)

smodel = createpde('structural','static-solid');
importGeometry(smodel,'Blade.stl');
msh = generateMesh(smodel,'Hmax',t);

E = 227E9; % in Pa
CTE = 12.7E-6; % in 1/K
nu = 0.27; 

structuralProperties(smodel,'YoungsModulus',E, ...
                            'PoissonsRatio',nu, ...
                            'CTE',CTE);
structuralBC(smodel,'Face',3,'Constraint','fixed');

p1 = a1*1e6; %in Pa
p2 = a2*1e6; %in Pa

structuralBoundaryLoad(smodel,'Face',11,'Pressure',p1); % Pressure side
structuralBoundaryLoad(smodel,'Face',10,'Pressure',p2);  % Suction side 

Rs = solve(smodel);

max_s = max(Rs.VonMisesStress/1e6); % mean, max, min, median, 

end
