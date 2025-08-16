function M = inertia_matrix_from_kinetic_energy(T_in, qdot_vector)
    % Takes as inputs:
    %   - T = the kinetic energy of the robot 
    %   - qdot_vector = the vector of q_dot ex: [qd1,qd2,qd3]
    %
    % Output:
    %   - M = the inertia matrix 
    if iscell(T_in)
        T = sum([T_in{:}]);      % somma dei contributi
    else
        T = T_in;
    end
    qdot_vector = qdot_vector(:).';    % assicurati sia riga
    M = simplify(hessian(T, qdot_vector));
    % forza simmetria numerica
    M = simplify( (M + M.')/2 );
end