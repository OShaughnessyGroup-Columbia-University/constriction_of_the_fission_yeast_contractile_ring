function [F1, s, t] = getSegForceGauss(k, l, rbead, kex, rex) 
% calculating excluded volume force on a pair of segments. Using scheme from
% geometrictools.com/Documentation/DistanceLine3Line3.pdf
% Code is at
% https://www.geometrictools.com/GTEngine/Include/Mathematics/GteDistSegmentSegment.h
% R(s, t) is distance between segments, where s and t are normalized distance
% along segments k and l
  rex2 = rex*rex;
  nbead = size(rbead,2);
  p0 = rbead(:, k-1);
  p1 = rbead(:, k);
  q0 = rbead(:, l-1);
  q1 = rbead(:, l);
  p10 = p1 - p0;
  q10 = q1 - q0;
  pq0 = p0 - q0;
  
  
  a = sum(p10.*p10); %= (p1 - p0)'*(p1 - p0) 
  b = sum(p10.*q10); %= (p1 - p0)'*(q1 - q0);
  c = sum(q10.*q10); %= (q1 - q0)'*(q1 - q0);
  d = sum(p10.*pq0);
  e = sum(q10.*pq0); %= (q1 - q0)'*(p0 - q0);
  f = sum(pq0.*pq0); %= (p0 - q0)'*(p0 - q0);

  f00 = d; f10 = d + a; f01 = d - b; f11 = d + a - b; % dR/ds at corners
  g00 = -e; g10 =-e - b; g01 =-e + c; g11 =-e - b + c; % dR/dt at corners

  % finds root of derivative of segment distance
  hatS = [getClampRoot(a, f00, f10);...
          getClampRoot(a, f01, f11)];
  classify = -1*(hatS<=0) + 1*(hatS>=1);

  if(all(classify==-1)) % mimimum distance is on s=0 edge
    params = [0; getClampRoot(c, g00, g01)];
  elseif(all(classify==1)) % minimum distance is on s=1 edge
    params = [1; getClampRoot(c, g10, g11)];
  else
    [edge, epts] = computeIntersect(hatS, classify);
    params = computeMinParams(edge, epts);
  end

  % pt of closest approach
  s = params(1);
  t = params(2);
  % squared distance at pt of closest approach
  R2 = a*s*s - 2*b*s*t + c*t*t + 2*d*s - 2*e*t + f;
  if(R2 < rex2 && R2 > 0)
      % unit vector from contact points
      n = (p0 + p10*s) - (q0 + q10*t);
      n = n/sqrt(sum(n.*n));
      F1 = n * kex * exp(-R2/rex2);
      %% mutual normal of two filaments
      %nx = cross(p10, q10);
      %nxmod = sqrt(sum(nx.*nx));
      %if(nxmod < 1e-6) % fils are parallel
      %    nx = pq0 - p10*sum(pq0.*p10)/sum(p10.*p10); % shortest vector cnncting segments
      %    nx = nx/sqrt(sum(nx.*nx));
      %else
      %    nx = nx/nxmod;
      %end
      %F1 = nx * (nx'*n) * kex * exp(-R2/rex2);
  else
      F1 = zeros(3, 1);
  end
  
%% Utility functions  
%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %
    function root = getClampRoot(sigma, h0, h1)
        if(h0 < 0)
            if(h1 > 0)
                root = -h0/sigma;
                if root > 1
                    root = 0.5;
                end
            else
                root = 1;
            end
        else
            root = 0;
        end
    end
%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %
  function [edge, epts] = computeIntersect(sval, classify)
      if(classify(1) < 0)
          edge(1) = 0;
          epts(1,:) = [0;  f00/b];
          if(epts(1,2) < 0 || 1 < epts(1,2))
              epts(1,2) = 0.5;
          end

          if classify(2) == 0
              edge(2) = 3;
              epts(2, :) = [sval(2), 1];
          else % classify(2) > 0
              edge(2) = 1;
              epts(2, :) = [1; f10/b];
              if(epts(2,2) < 0 || 1 < epts(2,2))
                  epts(2,2) = 0.5;
              end
          end
      elseif(classify(1) == 0)
          edge(1) = 2;
          epts(1, :) = [sval(1); 0];

          if(classify(2) < 0)
              edge(2) = 0;
              epts(2, :) = [0; f00/b];
              if(epts(2,2) < 0 || 1 < epts(2,2))
                  epts(2,2) = 0.5;
              end
          elseif(classify(2) == 0)
              edge(2) = 3;
              epts(2, :) = [sval(2); 1];
          else % classify(2) > 0
              edge(2) = 1;
              epts(2, :) = [1; f10/b];
              if(epts(2,2) < 0 || 1 < epts(2,2))
                  epts(2,2) = 0.5;
              end
          end
      else % classify(1) > 0
          edge(1) = 1;
          epts(1, :) = [1; f10/b];
          if(epts(1,2) < 0 || 1 < epts(1,2))
              epts(1,2) = 0.5;
          end

          if(classify(2) == 0)
              edge(2) = 3;
              epts(2, :) = [sval(2); 1];
          else
              edge(2) = 0;
              epts(2, :) = [0; f00/b];
              if(epts(2,2) < 0 || 1 < epts(2,2))
                  epts(2,2) = 0.5;
              end
          end
      end
  end
%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %
    function params = computeMinParams(edge, epts)
        delta = epts(2,2) - epts(1,2);
        h0 = delta * (-b*epts(1,1) + c*epts(1,2) - e);
        if(h0 >= 0)
            if(edge(1)==0)
                params = [0; getClampRoot(c, g00, g01)];
            elseif(edge(1)==1)
                params = [1; getClampRoot(c, g10, g11)];
            else
                params = epts(1, :);
            end
        else
            h1 = delta * (-b*epts(2,1) + c*epts(2,2) - e);
            if(h1 <= 0)
                if(edge(2)==0)
                    params = [0; getClampRoot(c, g00, g01)];
                elseif(edge(2)==1)
                    params = [1; getClampRoot(c, g10, g11)];
                else
                    params = epts(2, :);
                end
            else % h0 < 0 and h1 > 0
                z = min(max(h0/(h0-h1), 0), 1);
                params = (1-z)*epts(1,:) + z*epts(2,:);
            end
        end
    end
%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %%  %
end
