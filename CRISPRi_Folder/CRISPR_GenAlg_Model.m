function dx= CRISPR_GenAlg_Model(t,x,p)

global par 

dx = zeros(10,1);

crRNA=x(1); 
trRNA=x(2);
gRNA=x(3);
dCas9=x(4);  %%% free/inactive dCas9
Complex=x(5); %%% active CRISPR regulator
Pyrep=x(6); %%% CRISPR bound GFP promoter
Y=x(7); %%% GFP mRNA

yi=x(8); 
G=x(9);
Gm=x(10);  %%% measurable GFP

Py=par.Pytot-Pyrep; %%% free GFP promoter

% 1 alpha_cr, 2 alpha_tr, 3 deg_cr, 4 deg_tr, 5 deg_m, 6 omega, 7 gamma1, 8 gamma2, 9 alpha_m, 10 KI, 11 KE, 12 alpha_gm

dx(1) = p(1)*par.Pcr - p(3)*crRNA - p(7)*crRNA*trRNA; %%% crRNA

dx(2) = p(2)*par.Ptr - p(4)*trRNA - p(7)*crRNA*trRNA; %%% trRNA using same transcription rate

dx(3) = p(7)*crRNA*trRNA - p(8)*gRNA*dCas9; %%gRNA

dx(4) = - p(8)*gRNA*dCas9; %%% inactive/free Cas9

dx(5) = p(8)*gRNA*dCas9 - p(6)*Complex*Py; 

dx(6) = p(6)*Complex*Py;

dx(7) = p(9)*Py - p(5)*Y-p(10)*Y + p(11)*yi;

dx(8) = p(10)*Y - p(11)*yi;

dx(9) = p(11)*yi - p(12)*G;

dx(10) = p(12)*G;


return 