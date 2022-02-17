function [p0,p1] = evolve_node_numbers(n0,n1,n2,dtSim,...
						koff,alpha,epsilon,offFlag)

	% use our kinetic scheme to decide the probability that a
	% zero-node will become a 1-node (p0) in time dtSim in the simulation
	% and a 1-node will become a 2-node (p1).
	
	% IMPORTANT NOTE: these probabilities are meant to be interpreted
	% sequentially. i.e. first convert zero nodes to one nodes using
	% the probability p0, then convert the one to two nodes using
	% probability p1. DO NOT DO THIS SIMULTANEOUSLY.
	
	% kinetic scheme:
	% n0dot = -k1*n0 + koff*(n1+n2)
	% n1dot = k1*n0 - k2*n1 - koff*n1
	% n2dot = k2*n1 - koff*n2
	
	% if offFlag is set to false, the koff
	% terms in the equation above are set to zero.
	% in simulation, new nodes coming in and getting removed
	% are handled separately from this part of the code, so set
	% offFlag to false. Otherwise set to true. If not, all
	% nodes will become 2 nodes eventually.
	
	% alpha and epsilon are defined as
	% k1 = koff * alpha
	% k2 = k1 * epsilon
	
	% dtSim is the timestep of the simulation
	% n0,n1,n2 are the current numbers of 0,1,and 2-nodes.
	
	% node-ring binding and removal is handled elsewhere
	% within the simulation
	
	% first, choose a timestep that ensures successful integration.
	% rate of 0->1 is the highest rate in the simulation, so choose
	% a timestep such that very little zero to one conversion happens
	% in one time step.
	
	k1 = koff * alpha;
	k2 = k1 * epsilon;
	
	dt = 0.01/(k1); % ensure time step is small compared to the fastest rate
	
	n0_0 = n0;
	n1_0 = n1;
	n2_0 = n2;
	ntot_0 = n0_0 + n1_0 + n2_0;
	
	if ~offFlag
		koff = 0;
	end
	
	for k=0:dt:dtSim
		n0 = n0 + dt*(koff*(n1+n2)-k1*n0);
		n1 = n1 + dt*(k1*n0-k2*n1-koff*n1);
		n2 = n2 + dt*(k2*n1-koff*n2);
	end
	
	p0 = (n0_0 - n0)/n0_0;
	p1 = (n2 - n2_0)/(n1_0 + p0 * n0_0);
	
end