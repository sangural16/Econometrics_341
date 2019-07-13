[x4, y4] = meshgrid(-2:.1:2, -2:.1:2);
    z = x4.*exp(-x4.^2 - y4.^2);
    surf(x4,y4,z)
    title('Plotting a function of two variables')