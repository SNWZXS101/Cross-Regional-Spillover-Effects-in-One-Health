function dxdt = DWYode(tn,x,par1,par2)

dxdt=zeros(10,1);

d1 = par1(1) ;
d2 = par1(2) ;
m1 = par1(3) ;
m2 = par1(4) ;
beta = [par1(5),par1(5)] ; 
mu = par1(6) ;
gamma = par1(7) ; 

tn = floor(tn);

c1 = par2(tn,1);
c2 = par2(tn,2);
q1 = par2(tn,3);
q2 = par2(tn,4);
 
c = [c1 c2] ;
q = [q1 q2] ;
d = [d1 d2] ;
m = [m1 m2] ;



S = x(1:2)'; Sq = x(3:4)'; I = x(5:6)'; Iq = x(7:8)'; R = x(9:10)';
N = [sum(x(1:2:end)),sum(x(2:2:end))];
lambda = [c(1)*I(1)*S(1)/N(1), c(2)*I(2)*S(2)/N(2)] ;

dS = -lambda + (1-q).*(1-beta).*(1-m).*lambda + mu.*Sq + (1-flip(q)).*(1-flip(beta)).*flip(m).*flip(lambda)  ; 
dSq = q.*(1-beta).*lambda - mu*Sq ;
dI = (1-q).*beta.*(1-m).*lambda - gamma.*I + (1-flip(q)).*flip(beta).*flip(m).*flip(lambda) - d.*I + flip(d).*flip(I)   ;
dIq = q.*beta.*lambda - gamma.*Iq  ;
dR = gamma .* (I + Iq)  ;

% dS1 = -lambda(1) + (1-q(1)).*(1-beta(1)).*(1-m(1)).*lambda(1) + mu.*Sq(1) + (1-q(2)).*(1-beta(2)).*m(2).*lambda(2)  ; 
% dSq1 = q(1).*(1-beta(1)).*lambda(1) - mu*Sq(1) ;
% dI1 = (1-q(1)).*beta(1).*(1-m(1)).*lambda(1) - gamma.*I(1) + (1-q(2)).*beta(2).*m(2).*lambda(2) - d(1).*I(1) + d(2).*I(2)   ;
% dIq1 = q(1).*beta(1).*lambda(1) - gamma.*Iq(1)  ;
% dR1 = gamma .* (I(1) + Iq(1))  ;
% 
% dS2 = -lambda(2) + (1-q(2)).*(1-beta(2)).*(1-m(2)).*lambda(2) + mu.*Sq(2) + (1-q(1)).*(1-beta(1)).*m(1).*lambda(1)  ; 
% dSq2 = q(2).*(1-beta(2)).*lambda(2) - mu*Sq(2) ;
% dI2 = (1-q(2)).*beta(2).*(1-m(2)).*lambda(2) - gamma.*I(2) + (1-q(1)).*beta(1).*m(1).*lambda(1) - d(2).*I(2) + d(1).*I(1)   ;
% dIq2 = q(2).*beta(2).*lambda(2) - gamma.*Iq(2)  ;
% dR2 = gamma .* (I(2) + Iq(2))  ;

dxdt(1:2) = dS' ;
dxdt(3:4) = dSq' ;
dxdt(5:6) = dI' ;
dxdt(7:8) = dIq';
dxdt(9:10) = dR' ;


