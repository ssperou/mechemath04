%Runs fixed step numerical integration
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
%              have the form: dXdt = rate_func_in(t,X) (t is before X)
%step_func: the function used to approximate X(t) at the next step
%           [XB,num_evals] = step_func(rate_func_in,t,XA,h)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%       num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals] = ...
        explicit_RK_fixed_step_integration(rate_func_in,tspan,X0,h_ref,BT_struct)
        

        % extract t0 and tf
        t0 = tspan(1);
        tf = tspan(2);

        % find smallest N that's smaller than h_ref
        % Find smallest N (number of steps rounded up to nearest whole number) so that h <= h_ref
        N = ceil((tf - t0)/h_ref);
        % step size 
        h_avg = (tf - t0)/N;
        t_list = linspace(t0, tf, N+1)';
        X_list = zeros(N+1,length(X0));
        X_list(1,:) = X0(:)';
        num_evals = 0;
        

        for j = 1:N
            t = t_list(j);
            XA = X_list(j,:);
            [XB,iter] = explicit_RK_step(rate_func_in,t,XA',h_avg,BT_struct);
           
            X_list(j+1,:) = XB;
            num_evals = num_evals + iter;
        end
  
end