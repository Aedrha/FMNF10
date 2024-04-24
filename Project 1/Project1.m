%Problem 1
%vectors t and y
t = [1940 1950 1960 1970 1980 1990 2000]';
y = [6371432
    7041829
    7497967
    8081229
    8317937
    8590630
    8882792];
% define basis functions for elementwise
phi_a = @(x) (x*ones(1,7)).^[0:6];
phi_b = @(x) ((x-1940)*ones(1,7)).^[0:6];
phi_c = @(x) ((x-1970)*ones(1,7)).^[0:6];
phi_d = @(x) (((x-1970)/30)*ones(1,7)).^[0:6];
%pre-allocate matricies
Aa=zeros(7,7);
Ab=zeros(7,7);
Ac=zeros(7,7);
Ad=zeros(7,7);
%fill out matricies
for i=1:7
    Aa(i,1:7) = phi_a(t(i));
    Ab(i,1:7) = phi_b(t(i));
    Ac(i,1:7) = phi_c(t(i));
    Ad(i,1:7) = phi_d(t(i));
end
%store result in a cell to simplyfy future use:
Acell = {Aa Ab Ac Ad};
% calculate the condition numbers
basisfs = {'a', 'b', 'c', 'd'};
%l1 norm (maximum absolute column sum of the matrix)
for i=1:4
    condnum = cond(Acell{i},1);
     fprintf("Condition number with l1 norm for basis %s is: %f\n", basisfs{i}, condnum);
     fprintf('\n');
end
%l2 norm square root of largest eigenvalue
for i=1:4
    condnum = cond(Acell{i},2);
     fprintf("Condition number with l2 norm for basis %s is: %f\n", basisfs{i}, condnum);
     fprintf('\n');
end
%linf norm (maximum absolut row number)
for i=1:4
    condnum = cond(Acell{i},"inf");
     fprintf("Condition number with linf norm for basis %s is: %f\n", basisfs{i}, condnum);
     fprintf('\n');
end
%% Problem 2
clc
close all
%First getting the coefficients d) gave best condition number so using Ad
c = Ad \ y;
%redefining phi_d;
phi_d2 = @(x) (((x-1970)./30));
% %Defining p_6(t) on nested form
%THIS COULD BE IMPLEMENTED ITERATIVELY WILL FIX EVENTUALLY
p_6 = @(t) c(1) + phi_d2(t).*(c(2) + phi_d2(t).*(c(3) + phi_d2(t).*(c(4) + phi_d2(t).*(c(5) + phi_d2(t).*(c(6) + phi_d2(t).*c(7))))));
% %one year interval
t_eval = (1940:2000)';
% %get p_6_values on 1 year intervals
p_6_values = p_6(t_eval);
% 
% %plotting together with interpolation points
figure;
plot(t_eval, p_6_values, '-', 'LineWidth', 2);
hold on
plot(t, y, 'o', 'MarkerSize', 10);
hold off
set(gca, 'FontSize', 20) 
xlabel('Year','FontSize',20);
ylabel('Population','FontSize',20);
title('Interpolated Population Data',FontSize=22);
legend('Interpolating Polynomial (p_6(t))', 'Interpolation Points','Fontsize', 20);
grid on;
%%
% Problem 3
clc 
close all
t_eval = (1940:2010)';
p_6_values = p_6(t_eval);
t_extp = [t
    2010];
y_extp = [y
    9345135];
% 
% %plotting together with interpolation points
figure;
plot(t_eval, p_6_values, '-', 'LineWidth', 2);
hold on
plot(t_extp, y_extp, 'o', 'MarkerSize', 10);
hold off
set(gca, 'FontSize', 20) 
xlabel('Year','FontSize',20);
ylabel('Population','FontSize',20);
title('Interpolated Population Data',FontSize=22);
legend('Interpolating Polynomial (p_6(t))', 'Interpolation Points','Fontsize', 20);
grid on;
%%
% Problem 4
clc
close all;
%adding on for the p7 version
y_p7 = [y;
    9345135];
%A function to get the coefficients definition at bottom of script
Cs = NewtonCoeff(t,y);

%time vectors 1 year intervals
t_eval = (1940:2000);
t_eval10 = (1940:2010);

%x is used in the phi part of the newton interpolation
x = t_eval-t;

%need an extra row of ones becuase of matlabs prod function
x = [ones(1,size(x,2));x];

%function to get interpolation, definition at end of script
P6 = newton_interpolation(x,Cs);

%new x matrix for 1940-2010
x10 = t_eval10-[t;2010];
x10 = [ones(1,size(x10,2));x10];
%have to modify size of x
P610 = newton_interpolation(x10(1:end-1,:),Cs);

% plotting the 6:th degree polynomials
figure;
plot(t_eval, P6, '-', 'LineWidth', 2);
hold on
plot(t, y, 'o', 'MarkerSize', 10);
hold off
set(gca, 'FontSize', 20) 
xlabel('Year','FontSize',20);
ylabel('Population','FontSize',20);
title('Interpolated Population Data up to 2000',FontSize=22);
legend('Interpolating Newton Polynomial (p_6(t))', 'Interpolation Points','Fontsize', 20);
grid on;

figure;
plot(t_eval10, P610, '-', 'LineWidth', 2);
hold on
plot([t;2010], y_p7, 'o', 'MarkerSize', 10);
hold off
set(gca, 'FontSize', 20) 
xlabel('Year','FontSize',20);
ylabel('Population','FontSize',20);
title('Interpolated Population Data upto 2010',FontSize=22);
legend('Interpolating Newton Polynomial (p_6(t))', 'Interpolation Points','Fontsize', 20);
grid on;

% making the 7th degree
[p7,c7] = P7(x10,P610);
%plotting the seventh degree
figure;
plot(t_eval10, p7, '-', 'LineWidth', 2);
hold on
plot([t;2010], y_p7, 'o', 'MarkerSize', 10);
hold off
set(gca, 'FontSize', 20) 
xlabel('Year','FontSize',20);
ylabel('Population','FontSize',20);
title('Interpolated Population Data upto 2010',FontSize=22);
legend('Interpolating Newton Polynomial (p_7(t))', 'Interpolation Points','Fontsize', 20);
grid on;
%%
%Problem 5
clc
close all
%using spline, uses cubic polynomials by default
pline = spline((1940:10:2000),y,(1940:2000));
%plotting the seventh degree
figure;
plot(t_eval, pline, '-', 'LineWidth', 2);
hold on
plot(t, y, 'o', 'MarkerSize', 10);
hold off
set(gca, 'FontSize', 20) 
xlabel('Year','FontSize',20);
ylabel('Population','FontSize',20);
title('Interpolated Population Data upto 2000',FontSize=22);
legend('Interpolation using spline', 'Interpolation Points','Fontsize', 20);
grid on;

 %%
 %In short this function worls as follows:
 % Iteratively builds a lower triagular matrix using divided differences
 % The c-coefficients will be located on the diagonal of the matrix.
 % Should work on datasets of arbitrary size.
 function Cs = NewtonCoeff(time,pop)
    Cmat=[];
    for i = 1:length(time)
        Cmat = [Cmat; pop(i), zeros(1, length(pop)-1)];
        for j = 2:size(Cmat,1)
            Cmat(size(Cmat,1),j) = (Cmat(size(Cmat,1),j-1)-Cmat(size(Cmat,1)-1,j-1))/...
            (time(length(time))-time(length(time)-j+1));
        end
    end
Cs = diag(Cmat);
return;
end

 
%Recursive function for newton interpolation, takes x as a matrix
%containing all possiple (x-x_k). Needs an extra rows of ones for the prod
%function to work. Returns a row vector containing final estimations.
 function P = newton_interpolation(x,c)
    if size(c,1)==1
        P=c*ones(1,size(x,2));
    return;
    end
    Pnminusone = newton_interpolation(x(1:end-1,:),c(1:end-1));
    c_phi = c(end).*prod(x(1:end-1,:));
    
    Pn = Pnminusone +c_phi;
    P=Pn;
    return;
 end
%adding the 7th degree
 function [p7,c7] = P7(x,p6)
    f7=9345135;
    phi7 = prod(x(1:end-1,end));
    c7=(f7-p6(end))/phi7;
    p7=p6+c7.*prod(x(1:end-1,:));

 end





