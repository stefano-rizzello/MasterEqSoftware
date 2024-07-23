% Master Equations for River Networks:
% File example.red
% An example file for computing dispersion and diffusion
% of pollutants in river networks, using the Reduce library
% master_eq_river_network.red
% See https://reduce-algebra.sourceforge.io/
% for the Reduce Computer Algebra System.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
%    * Redistributions of source code must retain the relevant copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%      documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNERS OR
% CONTRIBUTORS
% BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% % *****************************************************************
% Authors and maintainers:
% S. Rizzello*, G. Napoli+, R. Vitolo**, S. De Bartolo*
% * Dipartimento di Ingegneria dell'Innovazione and
% ** Dipartimento di Matematica
% Universita' del Salento (Lecce, Italy)
% + and Dipartimento di Matematica
% Universita' Federico II, Napoli (Italy)
% email: samuele.debartolo@unisalento.it
% Date: 01 July 2024
% ===============================================================


operator rho;
for all t,s let rho(0,t,s)=0;
operator tprob;
for all m let tprob(m,0)=0;

in "master_eq_river_network.red";

net:={n1,{n2,{n3},{n4,{n5},{n6,{n7},{n8,{n9},{n10,{n11,{n12},{n13}},{n14,{n15},{n16,{n17},{n18,{n19,{n20},{n21}},{n22,{n23},{n24}}}}}}}}}}};

on rounded;

flatnet:=mkdepth_one(net);
validate_net1(net);
validate_net2(net);

for s:=1:length(flatnet) do rho(part(flatnet,s),t0,0):=0;
rho(n5,t0,0):=5;
rho(n15,t0,0):=4;
rho(n21,t0,0):=6;

tprob(n2,n1):=0.69;
tprob(n3,n2):=0.82;
tprob(n4,n2):=0.94;
tprob(n5,n4):=0.30;
tprob(n6,n4):=0.67;
tprob(n7,n6):=0.84;
tprob(n8,n6):=0.47;
tprob(n9,n8):=0.68;
tprob(n10,n8):=0.97;
tprob(n11,n10):=0.71;
tprob(n12,n11):=0.79;
tprob(n13,n11):=0.39;
tprob(n14,n10):=0.56;
tprob(n15,n14):=0.32;
tprob(n16,n14):=0.71;
tprob(n17,n16):=0.83;
tprob(n18,n16):=0.67;
tprob(n19,n18):=0.83;
tprob(n20,n19):=0.81;
tprob(n21,n19):=0.59;
tprob(n22,n18):=0.41;
tprob(n23,n22):=0.81;
tprob(n24,n22):=0.31;

tn:=15;
for i:=1:tn do all_master_equation(net,rho,t0,i,tprob);

rho(n1,t0,15);
results_n1:= for s:=0:tn collect {s,rho(n1,t0,s)};

end;