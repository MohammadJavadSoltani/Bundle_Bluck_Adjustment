% Mohammad Javad Soltani 9822663 
#Calc:EOP + Coods(Tie+ Z for horizantal GCPs + X,Y for Vertical GCPs)   
tic
close all; clear; clc; format long g;
%% First We Load Our Parameters Coumputed In Last Project And Our Data
load('parameters1.mat');%EOP_Images as "result_images"
load('parameters2.mat');%Tie Coordinates as "result_coordinates"
load('data.mat'); load('data_control_points.mat');%Main DATA
len_con = length(data_control_points); len_data = length(data);
plot(result_coordinates(:,2), result_coordinates(:,3), '*m', 'MarkerSize', 8);%Computed Control Points.
hold on;%we plot Tie coordinates this soon becuse it will change douring the App.
%% These Variables Are Primary Values
data_control_points(4, 4) = 641.272969739521;%surface1(x)
data_control_points(10, 4) = 657.615318655257;%surface1(y)
data_control_points(6, 2) = 1.817848254013564e+03;%vertical
data_control_points(6, 3) = 3.830066680064776e+03;%vertical  -->increase accuracy
Full_control_point = [data_control_points(1 : 3, :);data_control_points(5, :);...
                      data_control_points(7 : 9, :); data_control_points(11, :)];
Vertical_control_point = [data_control_points(6, :)];
Surface_control_point = [data_control_points(4, :); data_control_points(10, :)];

%% Now We Will Create Our Main Matrix
result_images(:, 8) = result_images(:, 6);%prevuse result was 1*6. 
result_images(:, 6 : 7) = zeros(14, 2);%creat omega and phi columns.
len_data = length(data);
syms X0 Y0 Z0 Omega Phi Kappa X Y Z x y
R = [cos(Kappa), -sin(Kappa), 0; sin(Kappa), cos(Kappa), 0; 0, 0, 1] * ...
    [cos(Phi), 0, sin(Phi); 0, 1, 0; -sin(Phi), 0, cos(Phi)] * ...
    [1, 0, 0; 0, cos(Omega), -sin(Omega); 0, sin(Omega), cos(Omega)];
%% build m, n, q, r, s
temp = R * [(X - X0); (Y - Y0); (Z - Z0)];%temp = A' = R(k)*R(phi)*R(omega)*(X-X0 ; Y-Y0 ; Z - Z0)
m = temp(1);%r11(X-X0) + r12(Y-Y0) + r13(Z-Z0) ; 
n = temp(2);%r21(X-X0) + r22(Y-Y0) + r23(Z-Z0) ;
q = temp(3);%r31(X-X0) + r32(Y-Y0) + r33(Z-Z0) ;
clear temp;
r = m / q;
s = n / q;%Photo_Tahlili_Dr.Valadan
%% Fx and Fy 
C = 153.692 / 1000; % mm -> m
Fx = x + C * r;
Fy = y + C * s;
%Now we need Jacobian of Fx & Fy based on Unknowns:

% F = [jacobian(Fx, [X0, Y0, Z0, Omega, Phi, Kappa]); jacobian(Fy, [X0, Y0, Z0, Omega, Phi, Kappa])]
% Fg_Tie = [jacobian(Fx, [X, Y, Z]); jacobian(Fy, [X, Y, Z])];
% Fg_Vertical = [jacobian(Fx, [X, Y]); jacobian(Fy, [X, Y])];
% Fg_Surface = [jacobian(Fx, Z); jacobian(Fy, Z)];

F = [ [ - (2768668935719305*cos(Kappa)*cos(Phi))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*sin(Phi)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), (2768668935719305*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*cos(Phi)*sin(Omega)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), (2768668935719305*cos(Omega)*cos(Phi)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2) - (2768668935719305*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))),   (2768668935719305*((Y - Y0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) + (Z - Z0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi))))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*(cos(Omega)*cos(Phi)*(Y - Y0) - cos(Phi)*sin(Omega)*(Z - Z0))*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), (2768668935719305*(cos(Kappa)*cos(Omega)*cos(Phi)*(Z - Z0) - cos(Kappa)*sin(Phi)*(X - X0) + cos(Kappa)*cos(Phi)*sin(Omega)*(Y - Y0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*(cos(Phi)*(X - X0) + cos(Omega)*sin(Phi)*(Z - Z0) + sin(Omega)*sin(Phi)*(Y - Y0))*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), -(2768668935719305*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0)))];
[ - (2768668935719305*cos(Phi)*sin(Kappa))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*sin(Phi)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), (2768668935719305*cos(Phi)*sin(Omega)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2) - (2768668935719305*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))), (2768668935719305*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*cos(Omega)*cos(Phi)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), - (2768668935719305*((Y - Y0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + (Z - Z0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi))))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*(cos(Omega)*cos(Phi)*(Y - Y0) - cos(Phi)*sin(Omega)*(Z - Z0))*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), (2768668935719305*(cos(Phi)*sin(Kappa)*sin(Omega)*(Y - Y0) - sin(Kappa)*sin(Phi)*(X - X0) + cos(Omega)*cos(Phi)*sin(Kappa)*(Z - Z0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*(cos(Phi)*(X - X0) + cos(Omega)*sin(Phi)*(Z - Z0) + sin(Omega)*sin(Phi)*(Y - Y0))*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2),  (2768668935719305*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0)))] ];
Fg_Tie = [[ (2768668935719305*cos(Kappa)*cos(Phi))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*sin(Phi)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), - (2768668935719305*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Phi)*sin(Omega)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2),   (2768668935719305*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Omega)*cos(Phi)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2)];...
[ (2768668935719305*cos(Phi)*sin(Kappa))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*sin(Phi)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2),   (2768668935719305*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Phi)*sin(Omega)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), - (2768668935719305*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Omega)*cos(Phi)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2)]];

Fg_Vertical = [[ (2768668935719305*cos(Kappa)*cos(Phi))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*sin(Phi)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), - (2768668935719305*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Phi)*sin(Omega)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2)];...
[ (2768668935719305*cos(Phi)*sin(Kappa))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*sin(Phi)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2),   (2768668935719305*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Phi)*sin(Omega)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2)]];

Fg_Surface = [  (2768668935719305*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Omega)*cos(Phi)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2);...
 - (2768668935719305*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Omega)*cos(Phi)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2)];
 

%% AS we've been Told Fx and Fy shouldn't be in the Loops! 
%now about A dimension
AEO = zeros(2 * len_data, 6 * 14);%EOP
Ag = zeros(2 * len_data, 50 * 3 + 1 * 2 + 2 * 1);%Unknown Tie and Controls coordinates
w = zeros(len_data * 2, 1);
delta_x = 1;
%% Calc X_cap
while norm(delta_x) > 10e-10
    for i = 1 : len_data
        image_number = data(i, 2);
        if data(i, 1) == 2
           image_number = image_number + 7;
        end% numbering Images we could do this out of while , but we should add another for like this.
        x = data(i, 4);
        y = data(i, 5);%Image coordinate
        X0 = result_images(image_number, 3);
        Y0 = result_images(image_number, 4);
        Z0 = result_images(image_number, 5);
        Omega = result_images(image_number, 6);
        Phi = result_images(image_number, 7);
        Kappa = result_images(image_number, 8);%Primary Values
        if data(i, 6) == 1
            num = data(i, 7) - 100;%Full Controls Starts at 100 
            X = Full_control_point(num, 2);
            Y = Full_control_point(num, 3);
            Z = Full_control_point(num, 4);
        elseif data(i, 6) == 2
            num  = data(i, 7) - 1000;%Surface Controls Starts at 1000 
            X = Surface_control_point(num, 2);
            Y = Surface_control_point(num, 3);
            Z = Surface_control_point(num, 4);
            Ag(2 * i - 1 : 2 * i, 50 * 3 + 2 + num) = eval(Fg_Surface);
        elseif data(i, 6) == 3
            num  = data(i, 7) - 10000;%Vertical Controls Starts at 10000 
            X = Vertical_control_point(num, 2);
            Y = Vertical_control_point(num, 3);
            Z = Vertical_control_point(num, 4);
            Ag(2 * i - 1 : 2 * i, 50 * 3 + 1 : 50 * 3 + 2) = eval(Fg_Vertical);
        elseif data(i, 6) == 0
            num  = data(i, 7);%Tie Points Starts at 1
            X = result_coordinates(num, 2);
            Y = result_coordinates(num, 3);
            Z = result_coordinates(num, 4);
            Ag(2 * i - 1 : 2 * i, num * 3 - 2 : num * 3) = eval(Fg_Tie);
        end
        AEO(2 * i - 1 : 2 * i, 6 * image_number - 5 : 6 * image_number) = eval(F);
        w(2 * i - 1) = eval(Fx);
        w(2 * i) = eval(Fy);
    end

    A = [AEO, Ag];
    delta_x = -inv(A' * A) * A' * w;%238->14*6 + 50*3 + 1*2(v) + 2*1(s) = 238

    for i = 1 : 14%Change EOP For All Images
       result_images(i, 3 : 8) = delta_x(i * 6 - 5 : i * 6, 1)' + result_images(i, 3 : 8); %C1 & C2 Are Ran And Photo Number ->const.
    end

    for i = 1 : 50
        result_coordinates(i, 2 : 4) = result_coordinates(i, 2 : 4) + delta_x(14 * 6 + 3 * i - 2 : 14 * 6 + 3 * i)';%14*6 for pass from EOPs.
    end

    Vertical_control_point(1, 2 : 3) = Vertical_control_point(1, 2 : 3) + delta_x(14 * 6 + 3 * 50 + 1 : 14 * 6 + 3 * 50 + 2)';%Just 1 V.Point with 2 Unknown
    Surface_control_point(1 : 2, 4) = Surface_control_point(1 : 2, 4) + delta_x(14 * 6 + 3 * 50 + 2 + 1 : 14 * 6 + 3 * 50 + 2 + 2);% 2 S.Points with 1 Unknown per each
end
disp('Norm delta_x: ');
disp(norm(delta_x));
%% Compute RMSE
data_control_points(4, 4) = 0;
data_control_points(10, 4) = 0;
data_control_points(6, 2) = 0;
data_control_points(6, 3) = 0;%for calc RMSE we should skip from GCPs that were involve in observation.

pair_points_matrix = zeros(20, 5);
temp = 1;
for i = 1 : len_con
    if data_control_points(i, 2) == 0
       continue;
    end%Jump from Elevational Points.
    temp2 = 1;
    for j = 1 : len_data
        if data(j, 3) == data_control_points(i, 1)
            if temp2 == 1
                pair_points_matrix(2 * temp - 1, :) = data(j, 1 : 5);
            elseif temp2 == 2
                pair_points_matrix(2 * temp, :) = data(j, 1 : 5);
            end
            temp2 = temp2 + 1;
        end
    end
    temp = temp + 1;
end
clear temp temp2;
%Now pair_point_matrix is full with 10 control points which are in 2 img.
for i = 1 : 10
   temp_image_number1 =  pair_points_matrix(2 * i - 1, 2);
   temp_image_number2 = pair_points_matrix(2 * i, 2);
   if pair_points_matrix(2 * i - 1, 1) == 2
      temp_image_number1 =  temp_image_number1 + 7;
   end
   if pair_points_matrix(2 * i, 1) == 2
      temp_image_number2 =  temp_image_number2 + 7;
   end
   
   omegal = result_images(temp_image_number1, 6);
   omegar = result_images(temp_image_number2, 6);
   phil = result_images(temp_image_number1, 7);
   phir = result_images(temp_image_number2, 7);
   kappal = result_images(temp_image_number1, 8);
   kappar = result_images(temp_image_number2, 8);
    Ml = [cos(kappal), -sin(kappal), 0;
        sin(kappal), cos(kappal), 0;
        0, 0, 1;] * [cos(phil), 0, sin(phil); 0, 1, 0; -sin(phil), 0, cos(phil)] * ...
        [1, 0, 0; 0, cos(omegal), -sin(omegal); 0, sin(omegal), cos(omegal)];
    Mr = [cos(kappar), -sin(kappar), 0;
        sin(kappar), cos(kappar), 0;
        0, 0, 1;] * [cos(phir), 0, sin(phir); 0, 1, 0; -sin(phir), 0, cos(phir)] * ...
        [1, 0, 0; 0, cos(omegar), -sin(omegar); 0, sin(omegar), cos(omegar)];
    L = Ml' * [pair_points_matrix(2 * i - 1, 4); pair_points_matrix(2 * i - 1, 5); -C];
    R = Mr' * [pair_points_matrix(2 * i, 4); pair_points_matrix(2 * i, 5);  -C];
    
    K = ((result_images(temp_image_number2, 3) - result_images(temp_image_number1, 3)) * L(2) - ...
    (result_images(temp_image_number2, 4) - result_images(temp_image_number1, 4)) * L(1)) / ...
    (R(2)*L(1) - L(2)*R(1));
    X(i) = K * R(1) + result_images(temp_image_number2, 3);
    Y(i) = K * R(2) + result_images(temp_image_number2, 4);
    Z(i) = K * R(3) + result_images(temp_image_number2, 5);
end
plot(X, Y, '*r', 'MarkerSize', 10);%Computed Control Points.
plot(data_control_points(:, 2), data_control_points(:, 3), '*g');%Origin Control Points.
plot(result_coordinates(:,2), result_coordinates(:,3) , 'ob')%Computed Tie points Coordinates.
legend('Tie Points','Control Points Computed', 'Control Points', 'Tie Points Computed');
clear temp_image_number1 temp_image_number2 kappal kappar Ml Mr L R K i j ;
X = X'; Y = Y'; Z = Z';

Ran_NUM = result_images(:,1) ; Photo_NUM = result_images(:,2) ; Omega =  result_images(:,3);
phi = result_images(:,4); Kappa = result_images(:,5); X0 = result_images(:,6);
Y0 = result_images(:,7); Z0 = result_images(:,8);
ID = result_coordinates(:,1); Xtie = result_coordinates(:,2);
Ytie = result_coordinates(:,3); Ztie = result_coordinates(:,4);

EOP = table(Ran_NUM,Photo_NUM,Omega,phi,Kappa,X0,Y0,Z0)
Tie_Coordinate = table(ID,Xtie,Ytie,Ztie)
fprintf('\n\n')

data_control_points(6, :) = [];
dX = data_control_points(:, 2) - X;
dY = data_control_points(:, 3) - Y;
dZ = data_control_points(:, 4) - Z;
dZ(4) = [];
dZ(8) = [];%remove surface and vertical from CALCs.
for i = 1 : length(dX)
    del(i, 1) = sqrt(dX(i) ^ 2 + dY(i) ^ 2);
end
RMSE = sqrt((del' * del) / (length(del) - 1))
%{
out put will be :
RMSE =

        0.0804346660992401

Elapsed time is 39.203486 seconds.
%}
toc