function [R, dRdx, e, dedx] = T_mat_OptimRotation(x,x0)
%
% This function determines the optimal rotation matrix to align a set of
% points x to a reference set x0.
% x and x0 should be [3 x n], with n>=3.

% Get optimal quaternions:
[e, dedx] = GetOptimalQuaternions(x,x0);

% Convert to rotation matrix:
[R]    = EvalQuatRot(e);
[dRdp] = EvalQuatRotDer(e);

% Rotation matrix derivatives
dRdx = zeros(3,3,numel(x));
for k = 1:numel(x)
   dRdx(:,:,k) = dRdp(:,:,1)*dedx(1,k) +dRdp(:,:,2)*dedx(2,k) +dRdp(:,:,3)*dedx(3,k) +dRdp(:,:,4)*dedx(4,k); 
end

end

%% Optimization function:
function [e,dedx] = GetOptimalQuaternions(x,x0)

    % Set up optimization matrices
    A = GetAMatrix(x, x0);
    
    % Solve eigenvalue problem
    [V,D] = eig(A);
    % the desired solution corresponds to the largest eigenvalue (why...)
    [lambda, iMax] = max(diag(D));
    
    % Select solution
    e = V(:,iMax);
    
    % Get derivatives wrt dofs q (derivative of eigenvalue problem)
    tempLeft  = [A-lambda*eye(4), -e; e', 0];
    dedx = zeros(4,numel(x)); 
    for k = 1:size(x,2)
        dAdx = GetAMatrix([1;0;0], x0(:,k));
        dAdy = GetAMatrix([0;1;0], x0(:,k));
        dAdz = GetAMatrix([0;0;1], x0(:,k));
        tempdx = tempLeft\[(lambda*eye(4)-dAdx)*e; 0];
        tempdy = tempLeft\[(lambda*eye(4)-dAdy)*e; 0];
        tempdz = tempLeft\[(lambda*eye(4)-dAdz)*e; 0];
        dedx(:,(k-1)*3+1) = tempdx(1:4);
        dedx(:,(k-1)*3+2) = tempdy(1:4);
        dedx(:,(k-1)*3+3) = tempdz(1:4);
    end
    
        

end

%% optimization matrix
function A = GetAMatrix(x,x0)
    
    ddRddp = EvalQuatRotDder(zeros(4,1)); % independent of p 

    A = zeros(4,4); 
    for i = 1:4
        for j = 1:4
            Atemp = sum(sum(x0.*(ddRddp(:,:,i,j)*x),1),2);
            A(i,j) = Atemp; 
        end
    end
    
end

%% Quaternion rotation matrices:
function [R] = EvalQuatRot(p)
a = p(1); b = p(2); c = p(3); d = p(4); 

R = [a^2+b^2-c^2-d^2, 2*b*c- 2*a*d, 2*b*d+2*a*c;
    2*b*c+2*a*d, a^2-b^2+c^2-d^2, 2*c*d-2*a*b;
    2*b*d-2*a*c, 2*c*d+2*a*b, a^2-b^2-c^2+d^2];
end

function dRdp = EvalQuatRotDer(p)
a = p(1); b = p(2); c = p(3); d = p(4); 

% Derivatives:
dRdp = zeros(3,3,4);

dRdp(:,:,1) = [2*a, -2*d, 2*c;
    2*d, 2*a, -2*b;
    -2*c, 2*b, 2*a];
dRdp(:,:,2) = [2*b, 2*c, 2*d;
    2*c, -2*b, -2*a;
    2*d, 2*a, -2*b];
dRdp(:,:,3) = [-2*c, 2*b, 2*a;
    2*b, 2*c, 2*d;
    -2*a, 2*d, -2*c];
dRdp(:,:,4) = [-2*d, -2*a, 2*b;
    2*a, -2*d, 2*c;
    2*b, 2*c, 2*d];
end

function ddRddp = EvalQuatRotDder(p)
a = p(1); b = p(2); c = p(3); d = p(4); 

% Derivatives:
ddRddp = zeros(3,3,4,4);

ddRddp(:,:,1,1) = [2, 0, 0;
    0, 2, 0;
    0, 0, 2];
ddRddp(:,:,1,2) = [0, 0, 0;
    0, 0, -2;
    0, 2, 0];
ddRddp(:,:,1,3) = [0, 0, 2;
    0, 0, 0;
    -2, 0, 0];
ddRddp(:,:,1,4) = [0, -2, 0;
    2, 0, 0;
    0, 0, 0];

ddRddp(:,:,2,1) = [0, 0, 0;
    0, 0, -2;
    0, 2, 0];
ddRddp(:,:,2,2) = [2, 0, 0;
    0, -2, 0;
    0, 0, -2];
ddRddp(:,:,2,3) = [0, 2, 0;
    2, 0, 0;
    0, 0, 0];
ddRddp(:,:,2,4) = [0, 0, 2;
    0, 0, 0;
    2, 0, 0];

ddRddp(:,:,3,1) = [0, 0, 2;
    0, 0, 0;
    -2, 0, 0];
ddRddp(:,:,3,2) = [0, 2, 0;
    2, 0, 0;
    0, 0, 0];
ddRddp(:,:,3,3) = [-2, 0, 0;
    0, 2, 0;
    0, 0, -2];
ddRddp(:,:,3,4) = [0, 0, 0;
    0, 0, 2;
    0, 2, 0];

ddRddp(:,:,4,1) = [0, -2, 0;
    2, 0, 0;
    0, 0, 0];
ddRddp(:,:,4,2) = [0, 0, 2;
    0, 0, 0;
    2, 0, 0];
ddRddp(:,:,4,3) = [0, 0, 0;
    0, 0, 2;
    0, 2, 0];
ddRddp(:,:,4,4) = [-2, 0, 0;
    0, -2, 0;
    0, 0, 2];

end



