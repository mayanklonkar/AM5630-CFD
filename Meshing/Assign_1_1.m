% Domain length
lx=1; 
ly=1;
% no. of elements
n_x=input('nx'); 
n_y=input('ny');

% For equidistant mesh
x_size= lx/n_x;
y_size= ly/n_y;

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


%loop for grid
for i=1:n_x
 
    grid_x(i+1)=grid_x(i)+x_size;

end

for i=1:n_y

    grid_y(i+1)=grid_y(i)+y_size;

end

%Loop for node location
for j=1:n_x
    for k=1:n_y
       node_x(j,k)=(grid_x(j)+grid_x(j+1))/2;
       node_y(j,k)=(grid_y(k)+grid_y(k+1))/2;
    end
end

[X,Y]=meshgrid(grid_x,grid_y);
figure
plot(X,grid_y); hold on
plot(transpose(grid_x),transpose(Y));
plot(node_x, node_y, 'k.');

plot(node_x,z_y,'k.');
plot(z_x,node_y,'k.');
plot(node_x,z_y2,'k.');
plot(z_x2,node_y,'k.');

xlabel('L_x')
ylabel('L_y')
title('Symmetric Mesh')



