clc; clear; clf;
N = 10;

x = unique(rand(N,1));
y = -1 + 2*rand(length(x),1);

my_axes = axes;
plot_data = plot(x,y,'o','DisplayName',"Data"); hold on

% OPTION 1: two for loops
% SLOW!!!
tic
V1 = zeros(N);
for i = 1:N
    for j = 1:N
        V1(i,j) = x(i)^(j-1);
    end
end
toc

% OPTION 2: one for loop, clever use of vectorization
% a little faster when N is large
tic
V2 = ones(N,1);
for i = 2:N
    V2(:,i) = V2(:,i-1).*x;
end
toc

% OPTION 3: fully vectorized
% inbetween, not too fast when N is large, but compact
tic
V3 = x.^(0:(N-1));
toc

% all options are equivalent mathematically, but computations
% might differ slightly. Errors are ~1e-16 (machine precision)
abs(V1 - V2)
abs(V2 - V3)

% we pick an option
VM = V3;

% vector c will be such that p(x) = c(1) + c(2)*x + c(3)*x^2 ...
% but polyval function works with the reverse ordering
% Hence, we flip things upside down!
c = flipud(VM\y);
p = @(t) polyval(c,t);

xx = linspace(0,1);
plot_fit = plot(xx,p(xx),'r--','DisplayName',"Interpolating polynomial");
hold off

%my_axes.YLim = [-1 1];
%my_axes.XLim = [0 1];
xlabel(my_axes,"$x$",'interpreter','latex');
ylabel(my_axes,"$y$",'interpreter','latex');

legend(my_axes,'show','interpreter','latex')
exportgraphics(my_axes,"interpolating_polynomial.pdf",'ContentType','vector')

fID = fopen("table_info.tex",'w');
fprintf(fID,"\\begin{tabular}{c");
fprintf(fID,"%s",repmat("|c",N));
fprintf(fID,"}\n$x_i$");
fprintf(fID,"& $%.2f$",x);
fprintf(fID,"\\\\\n\\hline$y_i$");
fprintf(fID,"& $%.2f$",y);
fprintf(fID,"\n\\end{tabular}");
fclose(fID);


