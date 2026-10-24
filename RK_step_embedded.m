%This function computes the value of X at the next time step
%for any arbitrary embedded RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%OUTPUTS:
%XB1: the approximate value for X(t+h) using the first row of the Tableau
%XB2: the approximate value for X(t+h) using the second row of the Tableau
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB1, XB2, num_evals] = RK_step_embedded(rate_func_in,t,XA,h,BT_struct)
    a = BT_struct.A; %current a value
    b = BT_struct.B; %current b value
    c = BT_struct.C; %current c value
    
    num_evals = 0;
    %Create K matrix
    m = length(XA); %should only equal 1
    s = length(c); %number of stages
    K = zeros(m, s);
    
    %Create for loop for K
    for i = 1:s
        ai = a(i, 1:i-1);
        ci = c(i);
        
        ti = t + ci*h;

        if i == 1
            Xi = XA;
        else
            sumi = K(:, 1:i-1)*ai';
            Xi = XA + h*sumi;
        end
        K(:, i) = rate_func_in(ti, Xi);
        num_evals = num_evals + 1;
    end

    sumb1 = K*b(1,:);
    sumb2 = K*b(2,:);

    XB1 = XA + h*sumb1;
    XB2 = XA + h*sumb2;
end