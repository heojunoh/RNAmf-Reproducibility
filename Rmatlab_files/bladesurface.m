smodel1 = createpde('structural', 'static-solid');
importGeometry(smodel1, 'Blade.stl');
pdegplot(smodel1, 'FaceLabels', 'on', 'FaceAlpha', 0.5);

smodel2 = createpde('structural', 'static-solid');
importGeometry(smodel2, 'Blade.stl');
pdegplot(smodel2, 'FaceLabels', 'on', 'FaceAlpha', 0.5);

smodel3 = createpde('structural', 'static-solid');
importGeometry(smodel3, 'Blade.stl');
pdegplot(smodel3, 'FaceLabels', 'on', 'FaceAlpha', 0.5);
