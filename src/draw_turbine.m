n6049_path = "data\NACA 6409\n6409.dat";

file = fopen(n6049_path, 'r');

format = '%f %f';
sizeA = [inf inf];
A = fscanf(file, format, sizeA);

fclose(file);