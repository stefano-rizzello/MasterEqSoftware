% Master Equations for River Networks:
% File master_eq_river_network.red
% a Reduce library for computing dispersion and diffusion
% of pollutants in river networks.
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
% Version: 2.0
% Date: 01 July 2024
% ===============================================================

% ALGORITHM

% 0 - Validation of the tree (abstact river network).
% 0.1 - First elements of any sublist of the tree must be nodes.
%       We use the index 0 to denote the successor of the last node
%       (in the sense of the lowermost node of the river basin).
%       That is the first element of the data structure.
%       The subsequent elements are the branches of the tree, each
%       connected in input by a given node.
% 0.2 - The list of all nodes shall not contain repeated nodes.
% 1 - Collect all nodes and all master equations at a given time t.
% 2 - Given master equations at time t, generate the equations at time t+dt.

% Loading the package assist, in order to manipulate lists
% in Reduce algebraic mode

load_package assist;

symbolic procedure h_alglistp(l);
% Check if l is an algebraic list
 (listp l) and ((car l) equal 'list);

symbolic operator h_alglistp;

%% algebraic procedure h_alglistp(l);
%%   if not(h_alglistp0(l)) then
%%     rederr "Error: not an algebraic list";

algebraic procedure validate_net1(net);
  % Validate the condition 0.1
  begin
    if net={} then return 0
    else
    if h_alglistp(first(net)) then
      rederr "The first item must be an atom!"
    else
      if length(net)=1 then return 0
    else
      for each el in rest(net) do validate_net1(el);
  end;

algebraic procedure validate_net2(net);
  % Validate the condition 0.2
  begin
    scalar flat_net,set_flat_net,temp_net;
    flat_net:=mkdepth_one(net);
    temp_flat_net:=flat_net;
    set_flat_net:=mkset(flat_net);
    if not(length(temp_flat_net)=length(set_flat_net))
      then rederr "The net contains repeated knots";
  end;

% 1 - Writing master equations at time t

algebraic procedure local_upknots(subnet);
  % 1.1 - Find all nodes that precede a given node
  begin
    if subnet={} then return subnet
    else if h_alglistp(first(subnet)) then
      rederr "Error: wrong hydrographic subnet";
    return first(subnet) . (for each el in rest(subnet) collect first(el))
  end;

algebraic procedure local_master_equation(downknot,subnet,rho_op,t0,n,tprob);
  % Here we compute the master equation in the first knot of subnet, f_subnet;
  % downknot is the (possible) knot after f_subnet (in hydrological sense);
  % rho_op is the density function
  % t0 is the initial time
  % dt is the time interval
  % tprob is the operator that yields the probability of transition,
  % tprob(a,b) gives the probability that there is a flow from a to b.
  begin
    scalar l_upknots;
    l_upknots:=local_upknots(subnet);
    rho_op(first(subnet),t0,n):=rho_op(first(subnet),t0,n-1) +
      (for each el in rest(subnet) sum
 	rho_op(first el,t0,n-1)*tprob(first el,first(subnet)))
	  - rho_op(first(subnet),t0,n-1)*tprob(first(subnet),downknot);
    rho_op(first(subnet),t0,n):=rho_op(first(subnet),t0,n);
  end;

algebraic procedure process_net(net);
% Here we 'open' a tree into the list of its branches
  if net={} then net
  else if h_alglistp(first(net)) then
    append(process_net(first(net)),process_net(rest(net)))
  else append(net,process_net(rest(net)));

algebraic procedure process_unfolded_net(unfolded_net,l_unet);
  % Here we read an 'open' list of branches from a tree and
  % return a list of branches
  begin
    scalar temp_branch,temppart,l_branches;
    l_branches:={};
    for i:=1:l_unet do
      if not(h_alglistp(temppart:=part(unfolded_net,i))) then
      <<
	temp_branch:=list(part(unfolded_net,i));
	j:=i+1;
        while (j<=l_unet and h_alglistp(part(unfolded_net,j))) do
	  <<
	    temp_branch:=part(unfolded_net,j) . temp_branch;
	    j:=j+1
	  >>;
	l_branches:=reverse(temp_branch) . l_branches;
      >>;
    return l_branches
  end;

algebraic procedure find_downknot(temp_bran,temp_iter_branches);
  begin
    scalar first_knot;
    first_knot:=first temp_bran;
    while not(temp_iter_branches={}) and
      not(first_knot member mkdepth_one(first(temp_iter_branches))) do
        temp_iter_branches:=rest(temp_iter_branches);
    return if temp_iter_branches={} then 0
      else first(first(temp_iter_branches))
  end;

% 2 - Given master equations at time t, find the master equations at time t+dt

algebraic procedure all_master_equation(net,rho_op,t0,n,tprob);
  begin
    scalar unfolded_net,l_branches,downknot,iter_branches,
      temp_bran,temp_iter_bran,cnt;
    cnt:=0;
    unfolded_net:=process_net(net);
    l_branches:=process_unfolded_net(unfolded_net,length(unfolded_net));
    iter_branches:=l_branches;
    while not(iter_branches={}) do
      <<
        temp_bran:=first(iter_branches);
        temp_iter_bran:=rest(iter_branches);
	downknot:=find_downknot(temp_bran,temp_iter_bran);
	local_master_equation(downknot,temp_bran,rho_op,t0,n,tprob);
	iter_branches:=temp_iter_bran
      >>
  end;


;end;

Local Variables:
mode:reduce
End: