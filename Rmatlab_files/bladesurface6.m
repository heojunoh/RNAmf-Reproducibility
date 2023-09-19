Rs1 = solve(smodel1);
pdeplot3D(smodel1, 'ColorMapData', Rs1.VonMisesStress, 'Deformation', Rs1.Displacement, 'DeformationScaleFactor', 100);
Rs2 = solve(smodel2);
pdeplot3D(smodel2, 'ColorMapData', Rs2.VonMisesStress, 'Deformation', Rs2.Displacement, 'DeformationScaleFactor', 100);
Rs3 = solve(smodel3);
pdeplot3D(smodel3, 'ColorMapData', Rs3.VonMisesStress, 'Deformation', Rs3.Displacement, 'DeformationScaleFactor', 100);
