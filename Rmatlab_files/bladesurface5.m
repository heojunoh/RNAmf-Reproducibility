p1 = 5e5;
p2 = 4.5e5;
structuralBoundaryLoad(smodel1, 'Face', 11, 'Pressure', p1);
structuralBoundaryLoad(smodel1, 'Face', 10, 'Pressure', p2);

structuralBoundaryLoad(smodel2, 'Face', 11, 'Pressure', p1);
structuralBoundaryLoad(smodel2, 'Face', 10, 'Pressure', p2);

structuralBoundaryLoad(smodel3, 'Face', 11, 'Pressure', p1);
structuralBoundaryLoad(smodel3, 'Face', 10, 'Pressure', p2);
