function [w, v, vc, T] = moving_frames_algorithm(num_joints, DH, qdots, m,  rc, prismatic_indices,I)
    % This function performs the moving frame algorithm and returns the
    % vectors w, v, vc, and the values of T for each joint
    %
    %inputs:
    % num_joints= an int reapresenting the number of links o fthe robot
    %
    % DH = the DH matrix of the robot in the order [alpha, a, d, theta ]
    %
    % qdts = a vector containing qdots ex: [qd1; qd2; qd3]
    %
    % m = a vector of masses [m1 ; ... ; mn]
    %
    %
    % rc = a vector containing the values ^ir_{ci} that is the distance
    % between the i-th CoM from the i-th reference frame
    % rc1 = [-l1+dc1; 0; 0];
    % rc2 = [-l2+dc2; 0; 0];
    % rc3 = [-l3+dc3; 0; 0];
    % rc = [rc1, rc2, rc3]; 
    %
    % prismatic indices= a list containing the list of prismatic indices,
    % ex for a PPR= [1,2]
    %
    % I = a vector of 3 by 3 matrices [I1,..,In]
    %
    %output:

    syms Ixx Iyy Izz real
    


    % Check that the DH is consistent with the indicated number of joints 
    if size(DH, 1) ~= num_joints
        error('Number of rows in DH matrix should be equal to the number of joints.');
    end


    % Check that rc is consistent with the number of joints 
    if size(rc,2) ~= num_joints
        error('Length of rc vector should be equal to the number of joints.');
    end

    % Display the number of joints and prismatic indices
    disp(['Number of joints: ', num2str(num_joints)]);
    

    %START THE ALGORITHM:
    v00=[0;0;0]
    w00=[0;0;0]

    % Initialize outputs
    w  = sym(zeros(3, num_joints));
    v  = sym(zeros(3, num_joints));
    vc = sym(zeros(3, num_joints));
    T  = cell(1, num_joints);         % non cell(num_joints)
    z  = sym([0;0;1]);
    w_prev = sym([0;0;0]);
    v_prev = sym([0;0;0]);

    
    
    %Start the loop
    for i = 1:num_joints
        In=I(1:3,(i-1)*3+1:i*3);
      
        
        if ismember(i, prismatic_indices)
            sigma = 1;
        else
            sigma =0;
        end
        A = DHmatrix(DH(i,1), DH(i,2), DH(i,3), DH(i,4), true);
        R=A(1:3,1:3);
        p=A(1:3,4);

        r_i=R.'*p;
        fprintf("===JOINT %d===\n", i)
        % computing omega
        if i ==1
            w_prev= w00;
        else
            w_prev= wi;
        end
        
        wi=R.'*(w_prev + (1- sigma)*[0;0;qdots(i)]);
        fprintf("\nThe value of W_%d is:\n", i);
        fprintf_vec(simplify(wi));
        w(:,i)= wi; 
        
        % computing v
        if i ==1
            v_prev= v00;
        else
            v_prev= vi;
        end 
        
        vi=R.'*(v_prev + sigma*[0;0;qdots(i)])+ cross(wi,r_i,1);
        fprintf("\nThe value of V_%d is:\n", i);
        fprintf_vec(simplify(vi));
        v(:,i) = vi;

        %computing vc
        vci=vi+cross(wi,rc(:,i));
        fprintf("\nThe value of Vc_%d is:\n", i);
        fprintf_vec(simplify(vci));
        vc(:,i)=vci;

        %computing T
        Ti=simplify(0.5*m(i)*(vci'*vci)+0.5*wi'*In*wi);
        fprintf('\nThe value of T_%d is:\n', i);
        fprintf_vec(simplify(Ti));
        T{i}=Ti;
        disp("")
    
    
    end
end 


function fprintf_vec(vec)
    if isnumeric(vec)
        vec = vec(:); 
        for i = 1:length(vec)
            fprintf('[ %g ]\n', vec(i));
        end
        return
    end

    if isa(vec, 'sym')
        sz = size(vec);
        if isscalar(vec)
            fprintf('[ %s ]\n', string(vec));
        elseif isvector(vec)
            vec = vec(:); 
            for i = 1:length(vec)
                fprintf('[ %s ]\n', string(vec(i)));
            end
        else
            for r = 1:sz(1)
                rowStr = join(string(vec(r,:)), ' ');
                fprintf('[ %s ]\n', rowStr);
            end
        end
        return
    end

    error('Vec è di tipo non supportato: %s', class(vec));
end
