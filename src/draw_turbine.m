n6049_path = "C:\Users\user\Desktop\Wind Turbine Analysis\data\NACA 6409\n6409.dat";

file = fopen("n6409.txt", 'r');

format = '%f %f';
sizeA = [inf inf];
A = fscanf(file, format, sizeA);

fclose(file);