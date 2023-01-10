clear all
close all

% Fundamental matrix estimation

% keble images A and C
pl = [418, 134;495, 128;571, 119;573, 199;497, 203;416, 205;424, 292;496, 282];
pr = [137, 126;211, 125;289, 126;287, 206;210, 203;134, 202;134, 289;211, 287];

iml = imread("keble_a.jpg");
imr = imread("keble_b.jpg");

figure(1)
imshow(iml)
hold on
plot(pl(:,1),pl(:,2),'*r')
hold off

figure(2)
imshow(imr)
hold on
plot(pr(:,1),pr(:,2),'*r')
hold off

u = pl(:,1); v=pl(:,2); ub = pr(:,1);vb = pr(:,2);
A = [u.*ub, u.*vb, u, v.*ub, v.*vb, v, ub, vb, ones(8,1)];

[U,S,V] = svd(A);
F = V(:,end);

F = [F(1),F(2),F(3);F(4), F(5),F(6);F(7),F(8),F(9)];
rank(F)
det(F)

% epipolar lines
l = F*[ub';vb';ones(1,8)]; ml = -l(1,:)./l(2,:); cl = -l(3,:)./l(2,:);
lb = F'*[u';v';ones(1,8)]; mr = -lb(1,:)./lb(2,:); cr = -lb(3,:)./lb(2,:);

x = 1:720;
figure(3)
imshow(iml)
hold on
plot(pl(:,1),pl(:,2),'*r')
for i=1:8
    plot(x,ml(i)*x+cl(i),'LineWidth',2)
end
hold off

figure(4)
imshow(imr)
hold on
plot(pr(:,1),pr(:,2),'*r')
for i=1:8
    plot(x,mr(i)*x+cr(i),'LineWidth',2)
end
hold off

% force det(F)=0
[U,S,V]=svd(F);
S(3,3)=0;
F=U*S*V';
rank(F)
det(F)

% epipolar lines
l = F*[ub';vb';ones(1,8)]; ml = -l(1,:)./l(2,:); cl = -l(3,:)./l(2,:);
lb = F'*[u';v';ones(1,8)]; mr = -lb(1,:)./lb(2,:); cr = -lb(3,:)./lb(2,:);


figure(5)
imshow(iml)
hold on
plot(pl(:,1),pl(:,2),'*r')
for i=1:8
    plot(x,ml(i)*x+cl(i),'LineWidth',2)
end
hold off

figure(6)
imshow(imr)
hold on
plot(pr(:,1),pr(:,2),'*r')
for i=1:8
    plot(x,mr(i)*x+cr(i),'LineWidth',2)
end
hold off

% normalise data
p1 = pl; p2=pr;
p1o = p1; p2o = p2;
mp1 = mean(p1); mp2 = mean(p2);
s1 = mean(sqrt(sum((p1-mp1)'.^2))); s2 = mean(sqrt(sum((p2-mp2)'.^2)));
T1 = [sqrt(2)/s1,0,0;0,sqrt(2)/s1,0;0,0,1]*[1,0,-mp1(1);0,1,-mp1(2);0,0,1];
T2 = [sqrt(2)/s2,0,0;0,sqrt(2)/s2,0;0,0,1]*[1,0,-mp2(1);0,1,-mp2(2);0,0,1];
p1 = T1*[p1';ones(1,8)]; p2 = T2*[p2';ones(1,8)];
p1 = p1(1:2,:)'; p2 = p2(1:2,:)';

u = p1(:,1); v=p1(:,2); ub = p2(:,1);vb = p2(:,2);
A = [u.*ub, u.*vb, u, v.*ub, v.*vb, v, ub, vb, ones(8,1)];

[U,S,V] = svd(A);
F = V(:,end);

F = [F(1),F(2),F(3);F(4),F(5),F(6);F(7),F(8),F(9)];

% force det(F)=0
[U,S,V]=svd(F);
S(3,3)=0;
F=U*S*V';

F = T1'*F*T2;
u = pl(:,1); v=pl(:,2); ub = pr(:,1);vb = pr(:,2);
% epipolar lines
l = F*[ub';vb';ones(1,8)]; ml = -l(1,:)./l(2,:); cl = -l(3,:)./l(2,:);
lb = F'*[u';v';ones(1,8)]; mr = -lb(1,:)./lb(2,:); cr = -lb(3,:)./lb(2,:);

% finding epipole
el = round(l(1:2,:)'\(-l(3,:)'));
er = round(lb(1:2,:)'\(-lb(3,:)'));

x=1:el(1);
figure(7)
imshow([iml,zeros(size(iml));zeros(size(iml)),zeros(size(iml))])
hold on
plot(pl(:,1),pl(:,2),'*r')
plot(el(1),el(2),'*b')
for i=1:8
    plot(x,ml(i)*x+cl(i),'LineWidth',2)
end
hold off

x=1:er(1);
figure(8)
imshow([imr,zeros(size(iml));zeros(size(iml)),zeros(size(iml))])
hold on
plot(pr(:,1),pr(:,2),'*r')
plot(er(1),er(2),'*b')
for i=1:8
    plot(x,mr(i)*x+cr(i),'LineWidth',2)
end
hold off

% shift origin to center
T = [1,0,-size(iml,2)/2;0,1,-size(iml,1)/2;0, 0, 1];
elt = T*[el;1]; ert = T*[er;1];

% Shift epipoples to x-axis
if elt(1) >= 0 a = 1; else a = -1; end
r1 = a*elt(1)/norm(elt); r2 = a*elt(2)/norm(elt);
R1 = [r1,r2,0;-r2,r1,0;0,0,1];
els = round(R1*elt);

if ert(1) >= 0 a = 1; else a = -1; end
r1 = a*ert(1)/norm(ert); r2 = a*ert(2)/norm(ert);
R2 = [r1,r2,0;-r2,r1,0;0,0,1];
ers = round(R2*ert);

% shift epipoles to infinity
G1 = [1,0,0;0, 1, 0;-1/els(1), 0, 1];
G2 = [1,0,0;0, 1, 0;-1/ers(1), 0, 1];

H2 = inv(T)*G2*R2*T;
H1 = inv(T)*G1*R1*T;


% transform image left
a = H1*[1;1;1]; a = round(a/a(3));
b = H1*[1;size(imr,1);1]; b = round(b/b(3));
c = H1*[size(imr,2);1;1]; c = round(c/c(3));
d = H1*[size(imr,2);size(imr,1);1]; d = round(d/d(3));
p1t = H1*[pl';ones(1,8)];
p1t = round(p1t./p1t(3,:)); p1t = p1t(1:2,:); p1t(2,:)= p1t(2,:)+400;
ll = H1*l;mll = -ll(1,:)./ll(2,:); cll = -ll(3,:)./ll(2,:);
xmin = min([a(1),b(1),c(1),d(1)]); xmax = max([a(1),b(1),c(1),d(1)]);
ymin = min([a(2),b(2),c(2),d(2)]); ymax = max([a(2),b(2),c(2),d(2)]);
imls = zeros(ymax+400,xmax,3);
for j=1:size(iml,1)
    for i=1:size(iml,2)
        t = H1*[i;j;1];
        t = round(t/t(3));
        imls(t(2)+400,t(1),:)=iml(j,i,:);
    end
end

x=1:size(imls,2); y=ones(size(x));
figure(9)
imshow(uint8(imls))
hold on 
plot(p1t(1,:),p1t(2,:),'*r')
for i=1:8
    plot(x(:),p1t(2,i)*y(:),'LineWidth',2)
end
hold off

% transform image right
a = H2*[1;1;1]; a = round(a/a(3));
b = H2*[1;size(imr,1);1]; b = round(b/b(3));
c = H2*[size(imr,2);1;1]; c = round(c/c(3));
d = H2*[size(imr,2);size(imr,1);1]; d = round(d/d(3));
p2t = H2*[pr';ones(1,8)];
p2t = round(p2t./p2t(3,:)); p2t = p2t(1:2,:); p2t(2,:)= p2t(2,:)+400;
lbr = H2*lb;mrr = -lbr(1,:)./lbr(2,:); crr = -lbr(3,:)./lbr(2,:);
xmin = min([a(1),b(1),c(1),d(1)]); xmax = max([a(1),b(1),c(1),d(1)]);
ymin = min([a(2),b(2),c(2),d(2)]); ymax = max([a(2),b(2),c(2),d(2)]);
imrs = zeros(ymax+400,xmax,3);
for j=1:size(iml,1)
    for i=1:size(iml,2)
        t = H2*[i;j;1];
        t = round(t/t(3));
        imrs(t(2)+400,t(1),:)=imr(j,i,:);
    end
end
x=1:size(imrs,2); y=ones(size(x));
figure(10)
imshow(uint8(imrs))
hold on 
plot(p2t(1,:),p2t(2,:),'*r')
for i=1:8
    plot(x(:),p2t(2,i)*y(:),'LineWidth',2)
end
hold off




