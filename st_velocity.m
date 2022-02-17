function [vbead, vmyo, vmyp, fc_sept, mm] = ...
	st_velocity(pcpq,m,force,nbead,nmyo,cvec,cdot,tau_v,anc_flag)
	
	if nbead == 0
		vbead = double.empty(3,0);
		vmyo = double.empty(3,0);
		vmyp = double.empty(3,0);
		fc_sept = double.empty(3,0);
		return;
	end
	
    mm = [m, pcpq'];
    temp = sparse(size(pcpq,1),size(pcpq,1));
    temp2 = [pcpq, temp];
    mm = [mm; temp2];
    h = - cvec./tau_v - cdot;
    fvec = [force,h]';

	%{
	% gpu code
	gmm = gpuArray(mm);
	gfvec = gpuArray(fvec);
    gqvec = gmm \ gfvec;
	qvec = gather(gqvec);
	%}
	
	qvec = mm \ fvec;
	
    qdot = qvec(1:numel(force));
    lambda = qvec(numel(force)+1:end);
	lambda(~anc_flag) = 0;
    fc_sept = - pcpq' * lambda;
    
    % distribute the qdot vector into vbead and vmyo
    vcat = vec2mat(qdot,3)';
    vbead = vcat(:,1:nbead);
    vmyo = vcat(:,(nbead+1):(nbead+nmyo));
    vmyp = vcat(:,(nbead+nmyo+1):end);
	
end