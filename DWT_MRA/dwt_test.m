close all

wname = 'haar';
n = 32;
A = zeros(1,n);
D = zeros(1,n);
Avec = zeros(n,2*n);
Dvec = zeros(n,2*n);

figure(1)
tiledlayout
clf
nexttile
hold on
for ii = 1:n
    a = A;
    d = D;
    a(ii) = 1;
    Avec(ii,:) = idwt(a,d,wname);
    plot(Avec(ii,:))
end

nexttile
hold on
for ii = 1:n
    a = A;
    d = D;
    d(ii) = 1;
    Dvec(ii,:) = idwt(a,d,wname);
    plot(Dvec(ii,:))
end



fun = zeros(1,n*2);
fun(8:12) = -1;

[cA,cD] = dwt(fun,wname);
% cD = zeros(size(cD));
rec = idwt(cA,cD, wname);



%%
A = [Avec',Dvec'];

v = fun';
x = A\v;

vv = A*x;

%%
cc = [Avec;Dvec];

x = sum(fun.*cc,2);

vvv = cc'*x;

%%
figure(2)
clf
hold on
plot(fun,'o')
plot(rec)
plot(vv)
plot(vvv)

legend('Original', 'dwt-idwt', 'least squares', 'inner product')