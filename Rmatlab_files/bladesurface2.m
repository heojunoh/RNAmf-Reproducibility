msh = generateMesh(smodel1, 'Hmax', 0.05);
pdeplot3D(smodel1);

msh = generateMesh(smodel2, 'Hmax', 0.025);
pdeplot3D(smodel2);

msh = generateMesh(smodel3, 'Hmax', 0.0125);
pdeplot3D(smodel3);
