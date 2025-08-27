function [w_0, w_i] = compute_angular_velocity(DH_time, com_index)
%% The function compute the angular velocity w.r.t. base frame and local frame leveraging on 
% Kinematic relationships. 
%% IS FUNDAMENTAL THAT THE DH TABLE CONTAINS THE GENERALIZED COORDINATES EXPRESS W.R.T. TIME.
t = sym('t','real');

T = eye(4);
for i = 1:com_index
    alpha = DH_time(i,1);
    a     = DH_time(i,2);
    d     = DH_time(i,3);
    theta = DH_time(i,4);
    Ti = DHmatrix(alpha, a, d, theta, true);  % ok anche se theta = q_i(t)
    T  = simplify(T * Ti);
end


R = T(1:3, 1:3); 
R_dot = diff(R, t);
S_w = simplify(R_dot * R.'); 
w_0 = [S_w(3,2); S_w(1,3); S_w(2,1)];

fprintf("\n Angular velocity with respect to the base frame: "); 
disp(simplify(w_0)); 

w_i = R.'*w_0;
fprintf("\n Angular velocity with respect to the local frame: "); 
disp(simplify(w_i)); 

end 

