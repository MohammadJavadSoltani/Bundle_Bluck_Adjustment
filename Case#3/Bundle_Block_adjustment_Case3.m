%% Project 3 --> SK_CASE3
#Calc:EOP + Coods(Tie + GCPs) -> GCPs are wieghted.
tic; close all; clear; clc; format long g
%% First We Load Our Parameters Coumputed In Last Project And Our Data
load('parameters1.mat');
load('parameters2.mat');
load('data.mat');
load('data_control_points.mat');
plot(result_coordinates(:,2), result_coordinates(:,3), '*m', 'MarkerSize', 8);%Computed Control Points.
hold on;%we plot Tie coordinates this soon becuse it will change douring the App.
%% These Variables Are Primary Values
data_control_points(4, 4) = 641.272969739521;
data_control_points(10, 4) = 657.615318655257;
data_control_points(6, 2) = 1.817848254013564e+03;
data_control_points(6, 3) = 3.830066680064776e+03;
compelete_control_point = [data_control_points(1 : 3, :);data_control_points(5, :); data_control_points(7 : 9, :); data_control_points(11, :)];
compelete_control_point_copy = compelete_control_point;
height_control_point = [data_control_points(6, :)];
height_control_point_copy = height_control_point;
plan_control_point = [data_control_points(4, :); data_control_points(10, :)];
plan_control_point_copy = plan_control_point;
C = 153.692 / 1000; % mm -> m
%% Now We Will Create Our Main Matrix
result_images(:, 8) = result_images(:, 6);%addup 2 columns for Omega and phi
result_images(:, 6 : 7) = zeros(14, 2);
len_data = length(data);
% syms X0 Y0 Z0 Omega Phi Kappa X Y Z x y
% R = [cos(Kappa), -sin(Kappa), 0; sin(Kappa), cos(Kappa), 0; 0, 0, 1] * ...
%     [cos(Phi), 0, sin(Phi); 0, 1, 0; -sin(Phi), 0, cos(Phi)] * ...
%     [1, 0, 0; 0, cos(Omega), -sin(Omega); 0, sin(Omega), cos(Omega)];
%
% temp = R * [(X - X0); (Y - Y0); (Z - Z0)];
% m = temp(1);
% n = temp(2);
% q = temp(3);
% clear temp;
% r = m / q;
% s = n / q;
% Fx = x + C * r;
% Fy = y + C * s;
%
% F = [jacobian(Fx, [X0, Y0, Z0, Omega, Phi, Kappa]); jacobian(Fy, [X0, Y0, Z0, Omega, Phi, Kappa])];
% Fg = [jacobian(Fx, [X, Y, Z]); jacobian(Fy, [X, Y, Z])];
% Fg_V = [jacobian(Fx, [X, Y]); jacobian(Fy, [X, Y])];
% Fg_Plan = [jacobian(Fx, Z); jacobian(Fy, Z)];

%The computed values are added to loop.
VARIANCE = [(eye(len_data * 2) * 7 / 1000000) .^ 2, zeros(len_data * 2, 8 * 3 + 2 * 2 + 1 * 1); zeros(29, len_data * 2), (eye(29) * 15 / 100) .^ 2];
for i=1:8%sigma for Z is diffrent.
    VARIANCE(len_data * 2 + 3*i,len_data * 2 + 3*i) = (20 / 100) .^ 2;%Z full
    if i<2
        VARIANCE(len_data * 2 + 3*8 + 2*2 +1) = (20 / 100) .^ 2;%Z Vertical
    end
end
P = inv(VARIANCE);

AEO = zeros(2 * len_data, 6 * 14);
Ag = zeros(2 * len_data, 50 * 3 + 1 * 2 + 2 * 1);
w = zeros(len_data * 2, 1);
delta_x = 1;
while norm(delta_x) > 10e-12
    for i = 1 : len_data
        image_num = data(i, 2);
        if data(i, 1) == 2
            image_num = image_num + 7;
        end
        
        x = data(i, 4);
        y = data(i, 5);
        X0 = result_images(image_num, 3);
        Y0 = result_images(image_num, 4);
        Z0 = result_images(image_num, 5);
        Omega = result_images(image_num, 6);
        Phi = result_images(image_num, 7);
        Kappa = result_images(image_num, 8);
        if data(i, 6) == 1%FULL -> It's waigted so we have Ag here too
            num = data(i, 7) - 100;
            X = compelete_control_point(num, 2);
            Y = compelete_control_point(num, 3);
            Z = compelete_control_point(num, 4);
            Ag(2 * i - 1 : 2 * i, 154 + num * 3 - 2 : 154 + num * 3) = [ (2768668935719305*cos(Kappa)*cos(Phi))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*sin(Phi)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), - (2768668935719305*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Phi)*sin(Omega)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2),   (2768668935719305*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Omega)*cos(Phi)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2);
                (2768668935719305*cos(Phi)*sin(Kappa))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*sin(Phi)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2),   (2768668935719305*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Phi)*sin(Omega)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), - (2768668935719305*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Omega)*cos(Phi)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2)];
        elseif data(i, 6) == 2 % Ag1 -> Z  /  Ag2 -> XY
            num  = data(i, 7) - 1000;
            X = plan_control_point(num, 2);
            Y = plan_control_point(num, 3);
            Z = plan_control_point(num, 4);
            Ag(2 * i - 1 : 2 * i, 50 * 3 + 2 + num) = [(2768668935719305*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Omega)*cos(Phi)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2);
                - (2768668935719305*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Omega)*cos(Phi)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2)];
            Ag(2 * i - 1 : 2 * i, 154 + 8 * 3 + num * 2 - 1 : 154 + 8 * 3 + num * 2) = [ (2768668935719305*cos(Kappa)*cos(Phi))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*sin(Phi)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), - (2768668935719305*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Phi)*sin(Omega)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2);
                (2768668935719305*cos(Phi)*sin(Kappa))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*sin(Phi)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2),   (2768668935719305*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Phi)*sin(Omega)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2)];
        elseif data(i, 6) == 3% Ag1 -> XY  /  Ag2 -> Z
            num  = data(i, 7) - 10000;
            X = height_control_point(num, 2);
            Y = height_control_point(num, 3);
            Z = height_control_point(num, 4);
            Ag(2 * i - 1 : 2 * i, 50 * 3 + 1 : 50 * 3 + 2) = [ (2768668935719305*cos(Kappa)*cos(Phi))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*sin(Phi)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), - (2768668935719305*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Phi)*sin(Omega)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2);
                (2768668935719305*cos(Phi)*sin(Kappa))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*sin(Phi)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2),   (2768668935719305*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Phi)*sin(Omega)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2)];
            Ag(2 * i - 1 : 2 * i, 154 + 3 * 8 + 2 * 2 + 1) = [(2768668935719305*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Omega)*cos(Phi)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2);
                - (2768668935719305*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Omega)*cos(Phi)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2)];
        elseif data(i, 6) == 0
            num  = data(i, 7);
            X = result_coordinates(num, 2);
            Y = result_coordinates(num, 3);
            Z = result_coordinates(num, 4);
            Ag(2 * i - 1 : 2 * i, num * 3 - 2 : num * 3) = [ (2768668935719305*cos(Kappa)*cos(Phi))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*sin(Phi)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), - (2768668935719305*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Phi)*sin(Omega)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2),   (2768668935719305*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Omega)*cos(Phi)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2);
                (2768668935719305*cos(Phi)*sin(Kappa))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*sin(Phi)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2),   (2768668935719305*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Phi)*sin(Omega)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), - (2768668935719305*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*cos(Omega)*cos(Phi)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2)];
        end
        AEO(2 * i - 1 : 2 * i, 6 * image_num - 5 : 6 * image_num) = [ - (2768668935719305*cos(Kappa)*cos(Phi))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*sin(Phi)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), (2768668935719305*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*cos(Phi)*sin(Omega)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), (2768668935719305*cos(Omega)*cos(Phi)*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2) - (2768668935719305*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))),   (2768668935719305*((Y - Y0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) + (Z - Z0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi))))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*(cos(Omega)*cos(Phi)*(Y - Y0) - cos(Phi)*sin(Omega)*(Z - Z0))*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), (2768668935719305*(cos(Kappa)*cos(Omega)*cos(Phi)*(Z - Z0) - cos(Kappa)*sin(Phi)*(X - X0) + cos(Kappa)*cos(Phi)*sin(Omega)*(Y - Y0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*(cos(Phi)*(X - X0) + cos(Omega)*sin(Phi)*(Z - Z0) + sin(Omega)*sin(Phi)*(Y - Y0))*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), -(2768668935719305*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0)));
            - (2768668935719305*cos(Phi)*sin(Kappa))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*sin(Phi)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), (2768668935719305*cos(Phi)*sin(Omega)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2) - (2768668935719305*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))), (2768668935719305*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*cos(Omega)*cos(Phi)*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), - (2768668935719305*((Y - Y0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + (Z - Z0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi))))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) - (2768668935719305*(cos(Omega)*cos(Phi)*(Y - Y0) - cos(Phi)*sin(Omega)*(Z - Z0))*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2), (2768668935719305*(cos(Phi)*sin(Kappa)*sin(Omega)*(Y - Y0) - sin(Kappa)*sin(Phi)*(X - X0) + cos(Omega)*cos(Phi)*sin(Kappa)*(Z - Z0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))) + (2768668935719305*(cos(Phi)*(X - X0) + cos(Omega)*sin(Phi)*(Z - Z0) + sin(Omega)*sin(Phi)*(Y - Y0))*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0))^2),  (2768668935719305*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0)))];
        w(2 * i - 1) = x + (2768668935719305*((Z - Z0)*(sin(Kappa)*sin(Omega) + cos(Kappa)*cos(Omega)*sin(Phi)) - (Y - Y0)*(cos(Omega)*sin(Kappa) - cos(Kappa)*sin(Omega)*sin(Phi)) + cos(Kappa)*cos(Phi)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0)));
        w(2 * i) = y + (2768668935719305*((Y - Y0)*(cos(Kappa)*cos(Omega) + sin(Kappa)*sin(Omega)*sin(Phi)) - (Z - Z0)*(cos(Kappa)*sin(Omega) - cos(Omega)*sin(Kappa)*sin(Phi)) + cos(Phi)*sin(Kappa)*(X - X0)))/(18014398509481984*(cos(Omega)*cos(Phi)*(Z - Z0) - sin(Phi)*(X - X0) + cos(Phi)*sin(Omega)*(Y - Y0)));
    end
    
    for i = 1 : 8%-I3*3 for Full
        Ag(376 + i * 3 - 2 : 376 + i * 3, 154 + i * 3 - 2 : 154 + i * 3) =  -eye(3);
    end
    
    for i = 1 : 2%-I2*2 for Horizental
        Ag(376 + 24 + i * 2 - 1 : 376 + 24 + i * 2, 154 + 24 + i * 2 - 1 : 154 + 24 + i * 2) =  -eye(2);
    end
    
    Ag(376 + 24 + 5, 154 + 24 + 5) = -eye(1);%-I1*1 for vertical
    
    AEO(376 + 24 + 5, 1) = 0;%Dimantion
    for i = 1 : 8
        w(376 + 3 * i - 2 : 376 + 3 * i, 1) =  compelete_control_point_copy(i, 2 : 4) - compelete_control_point(i, 2 : 4);
    end
    for i = 1 : 2
        w(376 + 24 + 2 * i - 1 : 376 + 24 + 2 * i, 1) =  plan_control_point_copy(i, 2 : 3) - plan_control_point(i, 2 : 3);
    end
    w(376 + 24 + 4 + 1, 1) = height_control_point_copy(1, 4) - height_control_point(1, 4);
    
    A = [AEO, Ag];
    delta_x = -inv(A' * P * A) * A' * P * w;
    for i = 1 : 14% 1:EOP
        result_images(i, 3 : 8) = delta_x(i * 6 - 5 : i * 6, 1)' + result_images(i, 3 : 8);
    end
    
    for i = 1 : 50% 2:Tie_Coords ->X,Y for V & Z for H
        result_coordinates(i, 2 : 4) = result_coordinates(i, 2 : 4) + delta_x(14 * 6 + 3 * i - 2 : 14 * 6 + 3 * i)';
    end
    height_control_point(1, 4) = height_control_point(1, 4) + delta_x(end);
    height_control_point(1, 2 : 3) = height_control_point(1, 2 : 3) + delta_x(14 * 6 + 3 * 50 + 1 : 14 * 6 + 3 * 50 + 2)';
    plan_control_point(1 : 2, 4) = plan_control_point(1 : 2, 4) + delta_x(14 * 6 + 3 * 50 + 2 + 1 : 14 * 6 + 3 * 50 + 2 + 2);
    plan_control_point(1 : 2, 2 : 3) = [delta_x(154 + 84 + 24 + 1), delta_x(154 + 84 + 24 + 2); delta_x(154 + 84 + 24 + 3), delta_x(154 + 84 + 24 + 4)] + plan_control_point(1 : 2, 2 : 3);
    for i = 1 : 8% 3.Full_Coords
        compelete_control_point(i, 2 : 4) =  compelete_control_point(i, 2 : 4) + delta_x(154 + 84 + 3 * i - 2 : 154 + 84 + 3 * i)';
    end
end
GCP = table;
GCP.ID = [compelete_control_point(:, 1); plan_control_point(:, 1); height_control_point(:, 1)];
GCP.X = [compelete_control_point(:, 2); plan_control_point(:, 2); height_control_point(:, 2)];
GCP.Y = [compelete_control_point(:, 3); plan_control_point(:, 3); height_control_point(:, 3)];
GCP.Z = [compelete_control_point(:, 4); plan_control_point(:, 4); height_control_point(:, 4)];
GCP
result_coordinates

disp('Norm delta_x: ');
disp(norm(delta_x));
%% Intersection -> No Need :)
data_control_points(4, 4) = 0;
data_control_points(10, 4) = 0;
data_control_points(6, 2) = 0;
data_control_points(6, 3) = 0;

pair_points_matrix = zeros(20, 5);
temp = 1;
for i = 1 : 11%Number of Controls
    if data_control_points(i, 2) == 0
        continue;
    end
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
clear temp_image_number1 temp_image_number2 kappal kappar Ml Mr L R K;
toc