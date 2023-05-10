clear; close all;
% Define the tau_1, tau_2, and K data
x = [0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75];
y = [0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75];
K = [53.69 61.94 82.38 65.49 49.71 17.85 42.83 14.71;
     98.42 46.87 109.41 99.40 7.01 16.71 20.70 1.88;
     41.81 6.32 20.75 31.51 6.11 26.88 33.71 13.48;
     149.19 11.47 0.63 14.88 8.84 73.17 40.83 29.96;
     140.93 30.31 1.04 0.92 2.81 34.85 3.31 0.24;
     105.74 1.27 10.58 0.21 0.04 0.57 2.92 7.09;
     99.05 12.11 0.12 0.97 5.09 6.90 0.65 1.29;
     164.42 7.38 13.35 10.88 8.53 2.22 3.26 0.73];
K = K.*10^-7;

% Create a grid of x and y values
[X,Y] = meshgrid(x,y);

% Plot the surface
figure(); surf(X,Y,K);

% Set the axis labels
xlabel('x');
ylabel('y');
zlabel('K');

%% Take the natural logarithm of z to transform it to log-normal scale
K_log = round(log(K),2);

% Assume the mean and standard deviation
mu_K_log = -13.86;
sigma2_K_log = 3.72;

% Plot the surface with the log-transformed z values
figure(); surf(X,Y,K_log-mu_K_log);

% Set the axis labels
xlabel('x');
ylabel('y');
zlabel('ln(K)-\mu_{ln(K)}');

%% Compute the directional sample correlation function
ny = length(y);
nx = length(x);
rho_y = zeros(ny,1);
rho_x = zeros(nx,1);

% Y direction
for j = 1:ny
    iter_j = j-1;
    prefactor = 1./(sigma2_K_log*(nx*(ny-iter_j )-1));
    for k = 1:nx
        for i = 1:ny-iter_j 
            l = i+iter_j ;
            X_ik = K_log(i,k)-mu_K_log;
            X_lk = K_log(l,k)-mu_K_log;
            rho_y(j) = rho_y(j)+prefactor*X_ik*X_lk;
        end
    end
end

% X direction
for j = 1:nx
    iter_j = j-1;
    prefactor = 1./(sigma2_K_log*(ny*(nx-iter_j )-1));
    for k = 1:ny
        for i = 1:nx-iter_j 
            l = i+iter_j ;
            X_ik = K_log(k,i)-mu_K_log;
            X_lk = K_log(k,l)-mu_K_log;
            rho_x(j) = rho_x(j)+prefactor*X_ik*X_lk;
        end
    end
end

% Isotropic sample correlation
n = max(ny,nx);
rho = zeros(n,1);
for j = 1:n
    iter_j = j-1;
    prefactor = 1./(sigma2_K_log*(ny*(nx-iter_j )+nx*(ny-iter_j )-1));
    for k = 1:nx
        for i = 1:ny-iter_j 
            l = i+iter_j ;
            X_ik = K_log(i,k)-mu_K_log;
            X_lk = K_log(l,k)-mu_K_log;
            rho(j) = rho(j)+prefactor*X_ik*X_lk;
        end
    end
    for k = 1:ny
        for i = 1:nx-iter_j 
            l = i+iter_j ;
            X_ik = K_log(k,i)-mu_K_log;
            X_lk = K_log(k,l)-mu_K_log;
            rho(j) = rho(j)+prefactor*X_ik*X_lk;
        end
    end
end

% Plot the sample correlation versus separation distance
figure();
hold on;
plot(x-x(1),rho_x);
plot(y-y(1),rho_y);
plot(y-y(1),rho);
legend('X direction','Y direction','Isotropic');

% Set the axis labels
xlabel('\tau');
ylabel('\rho(\tau)');

%% Define the exponential function
x = transpose(x-x(1));
exp_func = fittype('exp(-x/b)','independent','x','dependent','y');

% Fit the exponential function to the data
exp_fit = fit(x,rho,exp_func);

% Plot the data and the fits
figure();
plot(x,rho,'b');
hold on;
plot(exp_fit,'g');
legend('Data (isotropic)','Exponential Fit');
xlabel('\tau');
ylabel('\rho(\tau)');
saveas(gcf,'correlation_fitted.png')
fprintf('The estimated correlation length is %.3f m\n', exp_fit.b)