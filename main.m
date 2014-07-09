%Tests for Simultaenous Orthogonal Rotation Angles

close all
clear all
clc

%Build visualization first for SEQUENTIAL global rotations%

%% 1) NO ROTATION DATUM FOR VISUALIZATON %

x = [0 1; 0 0; 0 0]; %X unit vector
y = [0 0; 0 1; 0 0]; %Y unit vector
z = [0 0; 0 0; 0 1]; %Z unit vector

%Plot standard unit vectors
figure(1);
hold on;
plot3(x(1,:),x(2,:),x(3,:),'r--');
plot3(y(1,:),y(2,:),y(3,:),'b--');
plot3(z(1,:),z(2,:),z(3,:),'g--');
axis equal
axis([-1,1,-1,1,-1,1]);
view(-33, -24);

xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');

%% 2) BUILD ROTATION MATRICES %

thetax = pi/10;
cx = cos(thetax);
sx = sin(thetax);

thetay = pi/10;
cy = cos(thetay);
sy = sin(thetay);

Rx = [1   0   0   ;
      0   cx  -sx ;
      0   sx  cx ];

Ry = [cy  0   sy ;
      0   1   0  ;
      -sy 0   cy];
        
%% 3) ROTATE UNIT VECTORS FROM 1). ABOUT GLOBAL X THEN GLOBAL Y %

xb = x;
yb = y;
zb = z;

xb(:,2) = Ry*Rx*xb(:,2); %Ry then Rx suggests global rotations in reverse order
yb(:,2) = Ry*Rx*yb(:,2);
zb(:,2) = Ry*Rx*zb(:,2);

figure(2);
hold on;
plot3(xb(1,:),xb(2,:),xb(3,:),'r');
plot3(yb(1,:),yb(2,:),yb(3,:),'b');
plot3(zb(1,:),zb(2,:),zb(3,:),'g');
axis equal
axis([-1,1,-1,1,-1,1]);
view(-33, -24);

%% 4) ROTATE UNIT VECTORS FROM 1). ABOUT GLOBAL Y THEN GLOBAL Z (REVERSE OF 3))%

xc = x;
yc = y;
zc = z;
 
xc(:,2) = Rx*Ry*xc(:,2); %Rx then Ry suggests global rotations in reverse order
yc(:,2) = Rx*Ry*yc(:,2);
zc(:,2) = Rx*Ry*zc(:,2);

figure(3);
hold on;
plot3(xc(1,:),xc(2,:),xc(3,:),'r');
plot3(yc(1,:),yc(2,:),yc(3,:),'b');
plot3(zc(1,:),zc(2,:),zc(3,:),'g');
axis equal
axis([-1,1,-1,1,-1,1]);
view(-33, -24);

%% 5) BUILD INCREMENTAL ROTATION MATRICES %

N = 10^7;
nthetax = thetax/N;
ncx = cos(nthetax);
nsx = sin(nthetax);

nthetay = thetay/N;
ncy = cos(nthetay);
nsy = sin(nthetay);

nRx = [1   0    0   ;
       0   ncx  -nsx;
       0   nsx  ncx ];

nRy = [ncy  0   nsy ;
       0    1   0   ;
       -nsy 0   ncy ];

xd = x;
yd = y;
zd = z;

xe = xd;
ye = yd;
ze = zd;

Rd = eye(3);
Re = eye(3);

for i = 1:N
    
    Rd = nRx*nRy*Rd;
    Re = nRy*nRx*Re;
    %xd(:,2) = nRx*nRy*xd(:,2); %Order shouldn't really matter here
    %yd(:,2) = nRx*nRy*yd(:,2);
    %zd(:,2) = nRx*nRy*zd(:,2);
end

xd(:,2) = Rd*xd(:,2); %Order shouldn't really matter here
yd(:,2) = Rd*yd(:,2);
zd(:,2) = Rd*zd(:,2);
    
xe(:,2) = Re*xe(:,2); %Order shouldn't really matter here
ye(:,2) = Re*ye(:,2);
ze(:,2) = Re*ze(:,2);
    
figure(4);
hold on;
plot3(xd(1,:),xd(2,:),xd(3,:),'r');
plot3(yd(1,:),yd(2,:),yd(3,:),'b');
plot3(zd(1,:),zd(2,:),zd(3,:),'g');
axis equal
axis([-1,1,-1,1,-1,1]);
view(-33, -24);

figure(5);
hold on;
plot3(xe(1,:),xe(2,:),xe(3,:),'r');
plot3(ye(1,:),ye(2,:),ye(3,:),'b');
plot3(ze(1,:),ze(2,:),ze(3,:),'g');
axis equal
axis([-1,1,-1,1,-1,1]);
view(-33, -24);

%% 6) SORA MATRICES %

phi = [thetax; thetay; 0];
k = phi/norm(phi);
zeta = norm(phi);

c = cos(zeta);
s = sin(zeta);
v = 1-c;

R = [ k(1)^2*v + c        , k(1)*k(2)*v - k(3)*s, k(1)*k(3)*v + k(2)*s;
      k(1)*k(2)*v + k(3)*s, k(2)^2*v + c        , k(2)*k(3)*v - k(1)*s;
      k(1)*k(3)*v - k(2)*s, k(2)*k(3)*v + k(1)*s, k(3)^2*v + c        ];