x1 = 0:pi/100:2*pi;
y1 = sin(x1);
plot(x1, y1,'--')
xlabel('x')
ylabel('sinx')

y2 = sin(x1-0.25);
y3 = sin(x1-0.5);


plot(x1,y1,x1,y2,x1,y3)