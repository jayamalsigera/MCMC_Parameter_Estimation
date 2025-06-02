% function err = compute_error(params, x0, t, Vm, x_obs)
%     Br = params(1);
%     Bp = params(2);
%     kt = params(3);
%     km = params(4);
% 
%     % Simulate with proposed parameters
%     try
%         [~, x_sim] = ode45(@(t, x) rip_dynamics_est(t, x, Vm(t), Br, Bp, kt, km), t, x0);
%         err = sum((x_sim(:,1) - x_obs(:,1)).^2 + (x_sim(:,3) - x_obs(:,3)).^2); % theta and alpha only
%     catch
%         err = inf; % if simulation fails
%     end
% end


% function err = compute_error(params, x0, t, Vm, x_obs)
%     Br = params(1); Bp = params(2);
%     kt = params(3); km = params(4);
%     try
%         [~, x_sim] = ode45(@(t, x) rip_dynamics_est(t, x, Vm(t), Br, Bp, kt, km), t, x0);
%         err = sum((x_sim(:,1) - x_obs(:,1)).^2 + (x_sim(:,3) - x_obs(:,3)).^2);
%     catch
%         err = inf;
%     end
% end


function err = compute_error(params, x0, t, Vm, x_ref)
    % Simulate system using params
    x_sim = simulate_system(params, x0, t, Vm);  % You should define this
    alpha_sim = x_sim(:,3);
    alpha_ref = x_ref(:,3);

    err = mean((alpha_sim - alpha_ref).^2);  % MSE or use sqrt for RMSE
end