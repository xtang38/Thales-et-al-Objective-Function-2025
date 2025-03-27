function dx= STAR_GenAlg_Model(t,x,p)

global par 

dx = zeros(5,1);

S=x(1);  %%% STAR 
Y=x(2);  %%% mRNA
Yi=x(3);  
G=x(4);
Gm=x(5);  %%% measurable GFP

Py=par.PYtot;  %%% free GFP plasmid

dx(1) = p(1)*par.Ps - p(2)*S - p(4)*S*Py;  %%% STAR

dx(2) = p(5)*p(4)*S*Py - p(3)*Y - p(6)*Y + p(7)*Yi; %%% mRFP production

dx(3) = p(6)*Y - p(7)*Yi; %%% maturation

dx(4) = p(7)*Yi - p(8)*G;   %%% elongation 

dx(5) = p(8)*G;  %%%measurement

return