
mn = 50;
N = 1;  
 
[x,y,z]=meshgrid(linspace(-10,10,mn));
val =  x + y + z + 8;
 
isosurface(x,y,z,val,0)
 
xlabel( 'x' );
ylabel( 'y' );
zlabel( 'z' );
 
axis equal
grid on