
clear
close all

%input
X = [-1.7,1.5,3];
Z = [-2,-0.1,2.25,4];
H = [-1.2653262,-0.9453262,-4.3203262];
H_prime = [1.7,-1.5,-3];
n = length(X);

x_grid = [-2:0.1:4]';
n_p = length(x_grid);

u = H + (x_grid-X).*H_prime;

%ploting functions
figid = figure;
plot(x_grid,u)
hold on
plot([Z;Z],repmat([12;-8],1,n+1))

%compute integrals
a_coeff = H-X.*H_prime;
b_coeff = H_prime;
for j = 1:n
    A(j) = exp(a_coeff(j))/b_coeff(j)*(exp(b_coeff(j)*Z(j+1))-exp(b_coeff(j)*Z(j)));
end
%normalize A
A = A/sum(A);
%vectorized function
A3 = exp(a_coeff)./b_coeff.*(exp(b_coeff.*Z(2:end))-exp(b_coeff.*Z(1:end-1)));

%verify integrals symbiotically
syms x
for j = 1:n
    f = exp(H(j)+(x - X(j))*H_prime(j));
    A2(j) = double(int(f,Z(j),Z(j+1)));
end

A = A/sum(A)
A_temp = [0.462418397 0.529588445 0.007993158]

