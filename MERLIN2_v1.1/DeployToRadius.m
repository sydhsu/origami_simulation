% Stopping criterion for origami Flasher simulation
% ======================================================================= %
% Inputs: 
%    Node - initial nodal coordinates
%    U    - displacement of current time step icrm
%    icrm - time step
%
% Output:
%    flag - 1 if criterion met, 0 otherwise
% ======================================================================= %

function flag = DeployToRadius(Node,U,icrm)
    if icrm > 0
        Nodenw = Node;
        Nodenw(:,1) = Node(:,1)+U(1:3:end);
        Nodenw(:,2) = Node(:,2)+U(2:3:end);
        Nodenw(:,3) = Node(:,3)+U(3:3:end);
        norm2d = vecnorm(Nodenw(:,1:2),2,2);
        global R_planar;
        flag = any( norm2d >= R_planar*.99 );
        if(flag)
            fprintf('Stopping condition met. Radius: \n');
            fprintf(num2str(max(norm2d)));
            fprintf('\n')
        end
    else
        flag = 0;
    end
end
