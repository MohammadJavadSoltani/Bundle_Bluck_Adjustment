%% PRJ1_Amaliat_Karboord_Photo_1401
% Mohammad Javad Soltani 9822663 
% Samaneh Karami 9824513
tic
clear all;close all ; clc ;format long;
%% Read Data
DATA = xlsread('data.xlsx');
DATA(:, 4 : 5) = DATA(:, 4 : 5) / 1000;% converting X & Y to meter mm->m
data_control = xlsread('data_control.xlsx'); %Read Control points
%% Encode (Tie, Control, Surface And Vertical)
len_data = length(DATA);%Length of Data
len_con = length(data_control);%Length of Control points
DATA(:, 6) = 0; %Tie :0
for i = 1 : len_con
   for j = 1 : len_data
      if DATA(j, 3) == data_control(i, 1) && ...
              data_control(i, 2) ~= 0 && ...
              data_control(i, 4) ~= 0             %No need to cheac y->Full_Control : 1
          DATA(j, 6) = 1;
      elseif DATA(j, 3) == data_control(i, 1) && ...
              data_control(i, 2) ~= 0 && ...
              data_control(i, 4) == 0             %Surface : 2
          DATA(j, 6) = 2;
      elseif DATA(j, 3) == data_control(i, 1) && ...
              data_control(i, 2) == 0 && ...
              data_control(i, 4) ~= 0             %Vertical : 3
          DATA(j, 6) = 3;
      end
   end
end
%% Numbering Points
%{
Tie Points      -> 1,    2,    3,    ...
Control Points  -> 101,  102,  103,  ...
Surface Points  -> 1001, 1002, 1003, ...
Vertical Points -> 10001,10002,10003,...
%}
control_id_temp = unique(DATA(:, 3));%POINT_IDs_UNIQUE
len_temp = length(control_id_temp);
T = 1; %T:Tie
C = 101; %C:Control
H = 1001; %H:Horizantal
V = 10001; %V:Vertical
DATA(:, 7) = 0;
checkT_temp = false;
checkC_temp = false;
checkH_temp = false;
checkV_temp = false;
for i = 1 : len_temp
    if checkT_temp == true
       T = T + 1;
    end
    if checkC_temp == true
       C = C + 1;
    end
    if checkH_temp == true
       H = H + 1;
    end
    if checkV_temp == true
       V = V + 1;
    end
    checkT_temp = false;
    checkC_temp = false;
    checkH_temp = false;
    checkV_temp = false;
    
    for j = 1 : len_data
        if DATA(j, 6) == 0 && control_id_temp(i) == DATA(j, 3)%Tie_Point
           DATA(j, 7) = T;
           checkT_temp = true;
        elseif DATA(j, 6) == 1 && control_id_temp(i) == DATA(j, 3)%Control_Point
            DATA(j, 7) = C;
            checkC_temp = true;
        elseif DATA(j, 6) == 2 && control_id_temp(i) == DATA(j, 3)%Surface_Point
            DATA(j, 7) = H;
            checkH_temp = true;
        elseif DATA(j, 6) == 3 && control_id_temp(i) == DATA(j, 3)%Vertical_Point
            DATA(j, 7) = V;
            checkV_temp = true;
        end
    end 
end

clear checkC_temp checkH_temp checkT_temp checkV_temp T C H V;
clear control_id_temp len_temp;
%% Which Point In Which Image
% First We Need Tie Points IDs.
temp = [];
j = 1;
for i = 1 : len_data
    if DATA(i, 6) == 0
       temp(j) = DATA(i, 3);
       j = j + 1;
    end
end
tie_points_ID = unique(temp)';%Now we have a Unique ID from Tie points.
clear temp;
for i=1:len_data
    DATA(i,8) = DATA(i,2);
    if DATA(i,1)==2
        DATA(i,8) = DATA(i,2)+7;%img8 = [1 2] 
    end
end%in C8 we have number of Images.
clear i;
% We Have 2x7 Images And 50 Tie Points.
len_tie_points = length(tie_points_ID);
len_image = 7;
check_matrix = zeros(len_tie_points , 15); % 15 = 14 Images + 1 Column For Points IDs
check_matrix(:, 1) = tie_points_ID;
for i = 1 : len_image
    for j = 1 : len_data
        for k = 1 : size(check_matrix, 1)
           if DATA(j, 1) == 1 && DATA(j, 2) == i && ...
                DATA(j, 3) == check_matrix(k, 1)
                check_matrix(k, i + 1) = 1;
           elseif DATA(j, 1) == 2 && DATA(j, 2) == i && ...%Now secound RAN
                   DATA(j, 3) == check_matrix(k, 1)
                check_matrix(k, i + 8) = 1;
           end
        end
    end
end
clear temp2 i j k%->?? ???? ???? ? ???? ?????? ? ??? ??? ?? ?? ??????? ?? ???? ???? ????.(2???? ?????)
%?? ?? ?? ?? ???? ????? ??????.????? ??? ?? ?? ??? ????? ????? ??????.
%???? ?? ????? ????? ?????? ?? ???? ?? ?? 11 ?? 21 ???? ?....
%% Create Our Main Matrix
A = zeros(len_data * 2, 14 * 4 + 51 * 2);% 51 : 50T + 1E
for i = 1 : len_data
    A(2 * i - 1 : 2 * i, (DATA(i, 8) - 1) * 4 + 1 : ...%(DATA(i, 8) -> Photo number [1-14]
        (DATA(i, 8) - 1) * 4 + 4) = [DATA(i, 4), DATA(i, 5), 1, 0; ...% 4->X0, Y0, Z0, Kappa
        DATA(i, 5), -DATA(i, 4), 0, 1];
    if DATA(i, 6) == 0
        A(2 * i - 1 : 2 * i, 14 * 4 + 2 * DATA(i, 7) - 1 : ...
            14 * 4 + 2 * DATA(i, 7)) = -eye(2, 2);
    elseif DATA(i, 6) == 3
        A(2 * i - 1 : 2 * i, 14 * 4 + 101 : 14 * 4 + 102) = -eye(2, 2);%Since the numbering of the...
                       %elevation points began from 10001, we had to separate them from tie points.
    end
end
%% Create Our L Matrix
l = zeros(len_data * 2, 1);
for i = 1 : len_data
    for j = 1 : len_con
        if (DATA(i, 6) == 2 || DATA(i, 6) == 1) && DATA(i, 3) == ...%full or surface
                data_control(j, 1)
           l(2 * i - 1 : 2 * i, 1) = ...
               [data_control(j, 2), data_control(j, 3)];
        end
    end
end
clear i j
%% Compute Parameters
x_cap = inv(A' * A) * A' * l;
%% Create Result Matrices
% Number Of Images = 14
% First We Need Mean Height.
H_mean = sum(data_control(:, 4)) / 9; % Ignore Surface Control Points
                                             % 11 - 2 = 9
% And Also We Need Focal Length.
C = 153.692 / 1000; % mm -> m
% Now We Can Compute Exterior Oreanation Parameters For Each Image.
result_EOPs = zeros(14, 2 + 4); %  |Ran Number|Image Number|X0|Y0|Z0|K0|
for i = 1 : 14% 2*7=14 Images
    temp_image = i;
    temp_ran = 1;
    if i > 7
       temp_image = temp_image - 7; 
       temp_ran = 2;
    end
    % These Variables Are Temps
    a = x_cap(4 * i - 3);
    b = x_cap(4 * i - 2);
    c = x_cap(4 * i - 1);
    d = x_cap(4 * i);
    landa = sqrt(a ^ 2 + b ^ 2);
    
    result_EOPs(i, 1 : 5) = [temp_ran, temp_image, c, d, landa * C + H_mean];
    if b > 0 && a > 0 %computing Theta
       result_EOPs(i, 6) = atan(b / a);
    elseif b > 0 && a < 0
        result_EOPs(i, 6) = pi - atan(abs(b / a));
    elseif b < 0 && a < 0
        result_EOPs(i, 6) = pi + atan(abs(b / a));
    elseif b < 0 && a > 0
        result_EOPs(i, 6) = 2 * pi - atan(abs(b / a));
    end
end
clear a b c d land temp_image temp_ran;
% Omega = 0 And Phi = 0; So We Will Not Save Them In Our Matrix. 14*8 -> 14*6
% result_EOPs Will be like:
%                             |Ran Number|Image Number|X0|Y0|Z0|K0|

% Now We Need To Compute Tie Points Coordinates.
% We Will Compute Z Elements Later on with space intersection.
result_coordinates = zeros(50, 2 + 3);%The last C is for Z.
for i = 1 : 50
    result_coordinates(i, 1 : 4) = [tie_points_ID(i), i, x_cap(14 * 4 + 2 * i - 1), ...
        x_cap(14 * 4 + 2 * i)];
end
clear i
%% intersection & Compute RMSE 
pair_points_matrix = zeros(20, 5);%all CPs are at least in 2 images.which ones?
temp = 1;
for i = 1 : len_con
    if data_control(i, 2) == 0 %egnore Elevational.C.Ps
       continue;
    end
    temp2 = 1;
    for j = 1 : len_data
        if DATA(j, 3) == data_control(i, 1)
            if temp2 == 1
                pair_points_matrix(2 * temp - 1, :) = DATA(j, 1 : 5);
            elseif temp2 == 2
                pair_points_matrix(2 * temp, :) = DATA(j, 1 : 5);
            end
            temp2 = temp2 + 1;
            if temp2 ==3%we dont need 3 points!
                break;
            end
        end
    end
    temp = temp + 1;
end
clear temp temp2;
%{
As we were told, a Point may be seen in several photos,
which the first two Are Imp for us.
So far, we know which checkpoints are in which photos 
and also we know their coordinates 
(except for the Elevational control points, of course).
%}
for i = 1 : 10
   temp_image_number1 =  pair_points_matrix(2 * i - 1, 2);
   temp_image_number2 = pair_points_matrix(2 * i, 2);%separating points in 2 new temp.
   if pair_points_matrix(2 * i - 1, 1) == 2
      temp_image_number1 =  temp_image_number1 + 7;
   end
   if pair_points_matrix(2 * i, 1) == 2
      temp_image_number2 =  temp_image_number2 + 7;
   end
   %we have image numbers in temp_image_number1 & 2
   kappal = result_EOPs(temp_image_number1, 6);%Kappa_Left
   kappar = result_EOPs(temp_image_number2, 6);%Kappa_Right
    Ml = [cos(kappal), -sin(kappal), 0;%Momega=0 & Mphi=0-> M = Mkappa
          sin(kappal),  cos(kappal), 0;
          0,            0          , 1;];
    Mr = [cos(kappar), -sin(kappar), 0;
          sin(kappar),  cos(kappar), 0;
          0,            0,           1;];
    L = Ml' * [pair_points_matrix(2 * i - 1, 4); pair_points_matrix(2 * i - 1, 5); -C];%L=Ml*(x0,y0,-c).:
    R = Mr' * [pair_points_matrix(2 * i, 4); pair_points_matrix(2 * i, 5);  -C];%R=Mr*(x0,y0,-c).:
    
    K = ((result_EOPs(temp_image_number2, 3) - result_EOPs(temp_image_number1, 3)) * L(2) - ...
        (result_EOPs(temp_image_number2, 4)  - result_EOPs(temp_image_number1, 4)) * L(1)) / ...
        (R(2)*L(1) - L(2)*R(1));%now we have Kappa
    X(i) = K * R(1) + result_EOPs(temp_image_number2, 3);%R(1)=U,result_images(temp_image_number2, 3)=X0
    Y(i) = K * R(2) + result_EOPs(temp_image_number2, 4);%R(2)=V,result_images(temp_image_number2, 4)=Y0
    Z(i) = K * R(3) + result_EOPs(temp_image_number2, 5);%R(3)=W,result_images(temp_image_number2, 5)=Z0
end
clear temp_image_number1 temp_image_number2 kappal kappar Ml Mr L R K;
plot(X, Y, '*r');
hold on;
plot(data_control(:, 2), data_control(:, 3), '*g');
X = X';
Y = Y';
Z = Z';
data_control(6, :) = [];%erase Elevational CP.
dX = data_control(:, 2) - X;
dY = data_control(:, 3) - Y;
dZ = data_control(:, 4) - Z;
dZ(4) = [];%P4 is Surface.
dZ(8) = [];%P8 is Surface.
for i = 1 : length(dX)
    r(i, 1) = sqrt(dX(i) ^ 2 + dY(i) ^ 2);
end
RMSE= sqrt((r' * r) / (length(r) - 1))
%% Compute Tie Points Coordinates ( will needed for prj2;) )
pair_points_matrix2 = zeros(100, 5);
temp = 1;
for i = 1 : 50
    temp2 = 1;
    for j = 1 : len_data
        if DATA(j, 3) == check_matrix(i, 1)
            if temp2 == 1
                pair_points_matrix2(2 * temp - 1, :) = DATA(j, 1 : 5);
            elseif temp2 == 2
                pair_points_matrix2(2 * temp, :) = DATA(j, 1 : 5);
            end
            temp2 = temp2 + 1;
        end
    end
    temp = temp + 1;
end
clear temp temp2;
for i = 1 : 50
   temp_image_number1 =  pair_points_matrix2(2 * i - 1, 2);
   temp_image_number2 = pair_points_matrix2(2 * i, 2);
   if pair_points_matrix2(2 * i - 1, 1) == 2
      temp_image_number1 =  temp_image_number1 + 7;
   end
   if pair_points_matrix2(2 * i, 1) == 2
      temp_image_number2 =  temp_image_number2 + 7;
   end
   
   kappal = result_EOPs(temp_image_number1, 6);
   kappar = result_EOPs(temp_image_number2, 6);
    Ml = [cos(kappal), -sin(kappal), 0;
          sin(kappal),  cos(kappal), 0;
          0,            0,           1;];
    Mr = [cos(kappar), -sin(kappar), 0;
         sin(kappar),   cos(kappar), 0;
         0,             0,           1;];
    L = Ml' * [pair_points_matrix2(2 * i - 1, 4); pair_points_matrix2(2 * i - 1, 5); -C];
    R = Mr' * [pair_points_matrix2(2 * i, 4); pair_points_matrix2(2 * i, 5);  -C];
    
    K = ((result_EOPs(temp_image_number2, 3) - result_EOPs(temp_image_number1, 3)) * L(2) - ...
    (result_EOPs(temp_image_number2, 4) - result_EOPs(temp_image_number1, 4)) * L(1)) / ...
    (R(2)*L(1) - L(2)*R(1));
    XT(i) = K * R(1) + result_EOPs(temp_image_number2, 3);
    YT(i) = K * R(2) + result_EOPs(temp_image_number2, 4);
    ZT(i) = K * R(3) + result_EOPs(temp_image_number2, 5);
end
clear temp_image_number1 temp_image_number2 kappal kappar Ml Mr L R K i j;
plot(XT, YT, '*b');
legend('Control Points Computed', 'Control Points', 'Tie Points');
toc
% The End Of Program ^^.
%{
out put will be :
RMSE =

   7.161393716152796

Elapsed time is 2.115668 seconds.
%}