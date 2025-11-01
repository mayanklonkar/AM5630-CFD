% Domain length
lx=1; 
ly=1;
% no. of elements
n_x=input('nx'); 
n_y=input('ny');

% This specifies the the stretch direction 
% 1 for positive, -1 for negative and 0 for no stretch
e=input("x direction");
w=input('y direction');

x_size= lx/n_x;
y_size= ly/n_y;

% For equidistant mesh
dx= zeros(1,n_x);
dy= zeros(1,n_y);

z_x=zeros(1,n_x);
z_y=zeros(1,n_y);
z_x2=ones(1,n_x);
z_y2=ones(1,n_y);

% Defining grid array
grid_x=zeros(1,n_x+1);
grid_y=zeros(1,n_y+1);

%Defining nodes array
node_x=zeros(n_x,n_y);
node_y=zeros(n_x,n_y);

%Stretch Factor
st_x=0.9;
st_y=0.9;

%Length of first element using GP
dx(1)=(lx*(1-st_x)/(1-(st_x^n_x)));
dy(1)=(ly*(1-st_y)/(1-(st_y^n_y)));

for i=1:n_x

dx(i+1)=dx(i)*st_x;

end

for i=1:n_y

dy(i+1)=dy(i)*st_y;

end

%loop for grid
if e==-1

    for i=1:n_x
 
    grid_x(i+1)=grid_x(i)+dx((n_x+1)-i);
 
    end
end

if e==1
    for i=1:n_x
    % dx(i+1)=dx(i)*st_x;
    grid_x(i+1)=grid_x(i)+dx(i);
 
    end
end
if e==0


    for i=1:n_x
    % dx(i+1)=dx(i)*st_x;
    grid_x(i+1)=grid_x(i)+x_size;
 
    end

end



if w==1
    for i=1:n_y
       
        grid_y(i+1)=grid_y(i)+dy(i);
    
    end
end

if w==-1
     for i=1:n_y
       
        grid_y(i+1)=grid_y(i)+dy((n_y+1)-i);
     end
end

if w==0
     for i=1:n_y
  
    grid_y(i+1)=grid_y(i)+y_size;
 
     end
end

   

%Loop for node location
for j=1:n_x
    for k=1:n_y
       node_x(j,k)=(grid_x(j)+grid_x(j+1))/2;
       node_y(j,k)=(grid_y(k)+grid_y(k+1))/2;
    end
end

% Boundary point 

[X,Y]=meshgrid(grid_x,grid_y);
figure
plot(X,grid_y); hold on
plot(transpose(grid_x),transpose(Y));

% Plotting node points
plot(node_x, node_y, 'k.');

% Plotting corner points
plot(node_x,z_y,'k.');
plot(z_x,node_y,'k.');
plot(node_x,z_y2,'k.');
plot(z_x2,node_y,'k.');

xlabel('L_x')
ylabel('L_y')
title('Stretch Mesh')