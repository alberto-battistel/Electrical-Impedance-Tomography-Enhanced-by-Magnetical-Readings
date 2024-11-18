close all

I = 0.10; % A

r = 0.1; % m
L = 2*r*sin(2*pi/16/2); % m

x1 = -L/2;
x2 = +L/2;

N = 1e6;

p = zeros(N,3);

p(:,1) = x1 + (0:N-1)'*L/N;
dl = zeros(size(p));
dl(:,1) = (x2-x1)/N;



r0 = [0,0.01,0];

r1 = r0-p;

Br_ = cross(dl,r1,2)./vecnorm(r1,2,2).^3;
mu0 = 4*pi*1e-7; % H m-1
Br = mu0/(4*pi)*I*sum(Br_) % T

% Br =
% 
%    1.0e-05 *
% 
%          0         0    0.1780

