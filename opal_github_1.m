clear all; close all; clc; load opal.mat; % clear workspace & load data
b = .4:.01:1; csr = .1:.01:.5; % parameter space to scan
j = find(Y<-30); % define region -- this one is SO, others are: j = find(Y>30 & X<-80 | Y>30 & X>30); // j = find(abs(Y)<30.001 & X<-80 | abs(Y)<30.001 & X>30); // j = find(Y>30 & X>-80 & X<30);
nboot = 1000; % bootstrap iterations

z = Z(j); s = S(j); i = I(j); o = O(j); % refine dataset to region
z = z(o>0); i = i(o>0); s = s(o>0); o = o(o>0);
z = z(s>0); i = i(s>0); o = o(s>0); s = s(s>0);
z = z(i>0); o = o(i>0); s = s(i>0); i = i(i>0);
o = o(z>0); i = i(z>0); s = s(z>0); z = z(z>0);

clearvars -EXCEPT z s o i csr b nboot;

xc = i;
xs = s;
lo = log(o);
n = length(lo);

for k = 1:nboot; % find bootstrap parameter distribution

btstrp = 1:n;
yb = lo(btstrp);
xcb = xc(btstrp);
xsb = xs(btstrp);
zb = z(btstrp);


for i = 1:length(b);
    for j = 1:length(csr);
        y = yb+b(i).*(log(zb)-log(1000));
        x = log(xcb+csr(j).*xsb);
        %x = log(xsb);
        [m(i,j),a(i,j),r(i,j),sm(i,j),sb(i,j)] = lsqfitma(x,y);
        B(i,j) = b(i);
        CSR(i,j) = csr(j);
    end
end

[R(k),indx] = max(r(:));
bb(k) = B(indx);
csrb(k) = CSR(indx);
M(k) = m(indx);
Mu(k) = sm(indx);
A(k) = a(indx);
Au(k) = sb(indx);
k
end

r2_ = median(R.^2)
b_ = median(bb)
bu_ = mad(bb)
beta_ = median(csrb)
betau_ = mad(csrb)
gamma_ = median(M)
gammau_ = median(Mu)
kappa_ = median(exp(A))
kappau_ = .5*median(exp(A+Au)-exp(A-Au))

zs = z; os = o; is = xc; ss = s; % plot
y = log(os)+b_.*(log(zs)-log(1000));
x = log(kappa_.^gamma_.*(is+beta_.*ss));
sp0 = scatter(exp(x),exp(y),50,log10(is./ss),'filled');
hold on;
set(gca,'xscale','log','yscale','log','ticklabelinterpreter','latex','fontsize',16)
box on
axis([.005 500 .005 500])
    title('Southern Ocean ($r^2 = 0.79$)','interpreter','latex')
    ylabel('$\mathcal{F}_{OC} \times z^{0.69(\pm0.03)}$','interpreter','latex') 
    xlabel('$1.02(\pm0.03)\mathcal{F}_{IC}+ 0.20(\pm0.02)\mathcal{F}_{Si}$','interpreter','latex') 
    caxis([-2 2])
hold on
    spp = plot(logspace(-5,5),kappa_.*logspace(-5,5).^gamma_,'-.k','linewidth',3)
    spl = legend([spp],'$y\sim x^{0.86(\pm0.02)}$','fontsize',16,'location','southeast')
set(spl,'interpreter','latex')
set(gca,'xtick',[.01 1 100],'ytick',[.01 1 100])