function [vbead, vmyo, vmyp, fc, fanc] = velocity(pcpq,m,force,nbead,rmyo,cvec,cdot,tau,anc_flag)
% see platt_constraint_methods_thesis_caltech_89 p.97   
    mm = [m, pcpq'];
    temp = sparse(size(pcpq,1),size(pcpq,1));
    temp2 = [pcpq, temp];
    mm = [mm; temp2];
    h = - cvec/tau - cdot;
    fvec = [force,h]';
% disp('hello')
% size(mm)
% size(fvec)
    qvec = mm \ fvec;
    qdot = qvec(1:numel(force));
    lambda = qvec(numel(force)+1:end);
    fc = - pcpq' * lambda;
    lambda(~anc_flag) = 0;
    fanc = - pcpq' * lambda;
    
    % distribute the qdot vector into vbead and vmyo
    vcat = vec2mat(qdot,3)';
    vbead = vcat(:,1:nbead);
    vmyo = vcat(:,(nbead+1):(nbead+size(rmyo,2)));
    vmyp = vcat(:,(nbead+size(rmyo,2)+1):end);
end