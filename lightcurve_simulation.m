function lightcurve_simulation()
    % 種子碼
    rng(5132004, 'twister');

    % 初始化參數
    params = init_parameters();
    escaped_initial_positions = [];
    absorbed_initial_positions = [];
    escaped_times = [];
    absorbed_times = [];
    absorbed_pos = [];
    arc_photons_t = [];
    all_escaped_times = [];

    % 光變曲線記錄
    t_records = [];
    t_absorbed_records = [];

    for i = 1:params.N_photon
        photon = generate_photon(params);
        [escaped, t_arrival, escaped_pos, t_thermal_obs] = propagate_photon(photon, params);
        
        if escaped % 處理逃逸的光子
            all_escaped_times = [all_escaped_times; t_arrival];
            
            if strcmp(params.dimension, '2D')
                angle = atan2(escaped_pos(2), escaped_pos(1));
                in_fov = (angle >= params.arc_angle_range(1) && angle <= params.arc_angle_range(2));
            else
                escaped_dir = escaped_pos / norm(escaped_pos);
                cos_angle = dot(escaped_dir, params.view_direction);
                in_fov = (cos_angle >= cos(params.view_cone_angle/2));
            end
            
            if in_fov
                escaped_initial_positions = [escaped_initial_positions; photon.pos];
                escaped_times = [escaped_times; t_arrival];
                arc_photons_t = [arc_photons_t; t_arrival];
            end
            
            t_records(end+1) = t_arrival;

        else % 處理被吸收的光子
            if strcmp(params.dimension, '2D')
                angle = atan2(escaped_pos(2), escaped_pos(1));
                in_fov = (angle >= params.arc_angle_range(1) && angle <= params.arc_angle_range(2));
            else
                absorbed_dir = escaped_pos / norm(escaped_pos);
                cos_angle = dot(absorbed_dir, params.view_direction);
                in_fov = (cos_angle >= cos(params.view_cone_angle/2));
            end
            
            if in_fov
                absorbed_pos = [absorbed_pos; escaped_pos];
                absorbed_times = [absorbed_times; t_arrival];
                absorbed_initial_positions = [absorbed_initial_positions; photon.pos];
            end
        
            t_absorbed_records(end+1) = t_thermal_obs;
        end
    end

    % 統計圓弧光子 t<1 和 t≥1 的數量
    if ~isempty(arc_photons_t)
        num_t_less_1 = sum(arc_photons_t < 1);
        num_t_greater_1 = sum(arc_photons_t >= 1);
        
        fprintf('===== Arc photon statistics =====\n');
        fprintf('Number of photons with t < 1: %d (%.2f%%)\n', num_t_less_1, 100*num_t_less_1/length(arc_photons_t));
        fprintf('Number of photons with t > 1: %d (%.2f%%)\n', num_t_greater_1, 100*num_t_greater_1/length(arc_photons_t));
        fprintf('Total arc photon number: %d\n', length(arc_photons_t));
    else
        fprintf('No photons escape!\n');
    end

    % 統計所有光子 t<1 和 t≥1 的數量
    if ~isempty(all_escaped_times)
        num_all_t_less_1 = sum(all_escaped_times < 1);
        num_all_t_greater_1 = sum(all_escaped_times >= 1);
        total_escaped = length(all_escaped_times);
        
        fprintf('\n===== All photon statistics =====\n');
        fprintf('Number of photons with t < 1: %d (%.2f%%)\n', num_all_t_less_1, 100*num_all_t_less_1/total_escaped);
        fprintf('Number of photons with t > 1: %d (%.2f%%)\n', num_all_t_greater_1, 100*num_all_t_greater_1/total_escaped);
        fprintf('Total photon number: %d\n', total_escaped);
    else
        fprintf('No photons escape!\n');
    end

    % 計算熱光變曲線
    edges = 0:params.dt:params.T_max;
    nb = numel(edges)-1;
    N_abs = record_lightcurve(t_absorbed_records, params.T_max, params.dt);
    
    sum_N_abs = 0;
    P = zeros(1, nb);
    N_th = zeros(1, nb);
    for k = 1:nb
        sum_N_abs = sum_N_abs + N_abs(k);
        if params.cooling_on
            P(k) = sum_N_abs;
            N_c = (P(k) * params.dt);
            if N_c > 0
                N_th(k) = N_c;
            end
            sum_N_abs = sum_N_abs - N_c / params.cooling_time;
            if sum_N_abs < 0
                sum_N_abs = 0;
            end
        else
            P(k) = sum_N_abs;
            N_c = (P(k) * params.dt);
            if N_c > 0
                N_th(k) = N_c;
            end
        end
    end

    % 繪製光變曲線
    plot_lightcurves(t_records, t_absorbed_records, N_th, params);

    % 光子初始位置可視化
    if params.plot_initial_positions && ~isempty(escaped_initial_positions)
        % 繪製光子初始位置的可視化圖
        figure;
        hold on;
        
        % 繪製 R_outer, R_inner
        theta = linspace(0, 2*pi, 100);
        plot(params.R_outer * cos(theta), params.R_outer * sin(theta), 'k-', 'LineWidth', 1.5);
        plot(params.R_inner * cos(theta), params.R_inner * sin(theta), 'k-', 'LineWidth', 1.5);
        
        % 標註觀察的圓弧
        arc_theta = linspace(params.arc_angle_range(1), params.arc_angle_range(2), 50);
        plot(params.R_outer * cos(arc_theta), params.R_outer * sin(arc_theta), 'r-', 'LineWidth', 2);
        
        % 繪製逃逸光子的初始位置
        scatter(escaped_initial_positions(:,1), escaped_initial_positions(:,2), 20, escaped_times, 'filled');
        colorbar;
        colormap(jet);
        clim([min(escaped_times), max(escaped_times)]);
        xlabel('X'); ylabel('Y');
        title('Initial positions of escaped photons');
        axis equal;
        hold off;

        % 繪製吸收光子初始位置可視化圖
        figure;
        hold on;
        
        theta = linspace(0, 2*pi, 100);
        plot(params.R_outer * cos(theta), params.R_outer * sin(theta), 'k-', 'LineWidth', 1.5);
        plot(params.R_inner * cos(theta), params.R_inner * sin(theta), 'k-', 'LineWidth', 1.5);

        arc_theta = linspace(params.arc_angle_range(1), params.arc_angle_range(2), 50);
        plot(params.R_outer * cos(arc_theta), params.R_outer * sin(arc_theta), 'r-', 'LineWidth', 2);

        if ~isempty(absorbed_initial_positions) && ~isempty(absorbed_times)
            scatter(absorbed_initial_positions(:,1), absorbed_initial_positions(:,2), 20, absorbed_times, 'filled');
            colorbar;
            colormap(jet);
            clim([min(absorbed_times), max(absorbed_times)]);
        end
        
        title('Initial positions of absorbed photons');
        xlabel('X'); ylabel('Y');
        axis equal;
        hold off;
    end
end

function params = init_parameters()
    % 初始化參數
    params.R_inner = 0.8;               % cocoon內半徑
    params.R_outer = 1;                 % cocoon外半徑
    params.R_inj = 0.2;                 % sphere的半徑（僅在 sphere 模式下使用）
    params.ring_inner = 0.75;           % ring的內半徑（僅在 ring 模式下使用）
    params.ring_outer = 0.75;           % ring的外半徑（僅在 ring 模式下使用）
    params.density = 5;                 % 環狀介質均勻密度
    params.kappa = 1;                   % 不透明度
    params.N_photon = 100000;           % 總模擬光子數
    params.dt = 0.01;                   % 時間 bin 寬度
    params.T_max = 2.0;                 % 最大觀測時間
    params.injection_mode = 'point';    % 注入模式: 'point', 'sphere', 'ring'
    params.flare_duration = 0.2;        % 耀班持續時間

    params.dimension = '2D';              % 2D or 3D
    params.arc_angle_range = [0, pi/6];   % for 2D, 想觀測的圓弧角度範圍 -pi to pi
    params.view_direction = [0, 0, 1];    % for 3D, 想觀測方向單位向量
    params.view_cone_angle = pi/6;        % for 3D, 想觀測圓錐角
    params.FOV2D = pi - abs(params.arc_angle_range(2) - params.arc_angle_range(1));   % 2D的FOV
    params.FOV3D = pi - params.view_cone_angle;                                       % 3D的FOV

    % cooling
    params.cooling_on = false;   % true有冷卻，false無冷卻
    params.cooling_time = 0.1;   % cooling time

    % thermal re-emit
    params.thermal_reemit_dir = false;  % false保持原來方向，true各向同性發射
    params.thermal_reemit_delay = 0;    % 再發射的額外延遲

    % visualization
    params.plot_initial_positions = false;   % 是否繪製光子初始位置可視化
end

function photon = generate_photon(params)

    switch params.injection_mode
        case 'point'
            r = 0; % 從中心注入
        case 'sphere'
            if strcmp(params.dimension, '2D')
                r = sqrt(rand()) * params.R_inj;
            else
                r = (rand())^(1/3) * params.R_inj;
            end
            if r > params.R_inner
                error('R_inj 必須小於等於 R_inner');
            end
        case 'ring'
            % 從有厚度的圓環注入
            if params.ring_inner > params.ring_outer
                error('ring_inner 必須小於等於 ring_outer');
            end
            if params.ring_outer > params.R_inner
                error('ring_outer 必須小於等於 R_inner');
            end
            % 使用平方根方法確保均勻分佈在圓環內
            r = sqrt(rand() * (params.ring_outer^2 - params.ring_inner^2) + params.ring_inner^2);
        otherwise
            error('未知的注入模式');
    end
    
    if strcmp(params.dimension, '2D')
        theta = rand() * 2 * pi;
        x = r * cos(theta);
        y = r * sin(theta);
        
        angle = rand() * 2 * pi;
        dx = cos(angle);
        dy = sin(angle);
        
        photon.pos = [x, y];
        photon.dir = [dx, dy];
        
    else
        theta = rand() * 2 * pi;
        cos_phi = 2 * rand() - 1;
        phi = acos(cos_phi);
        
        x = r * sin(phi) * cos(theta);
        y = r * sin(phi) * sin(theta);
        z = r * cos(phi);
        
        theta_dir = rand() * 2 * pi;
        cos_phi_dir = 2 * rand() - 1;
        phi_dir = acos(cos_phi_dir);
        
        dx = sin(phi_dir) * cos(theta_dir);
        dy = sin(phi_dir) * sin(theta_dir);
        dz = cos(phi_dir);
        
        photon.pos = [x, y, z];
        photon.dir = [dx, dy, dz];
    end

    photon.t_emit = rand() * params.flare_duration;
end

function [escaped, t_arrival, escaped_pos, t_thermal_obs] = propagate_photon(photon, params)

    pos = photon.pos;
    dir = photon.dir;

    % 計算到內外邊界的路徑長
    [s_inner, s_outer] = compute_path_lengths(pos, dir, params);

    t_thermal_obs = NaN;
    escaped_pos = [NaN, NaN];

    r = sqrt(dot(pos,pos));
    if r < params.R_inner
        if isnan(s_inner) || s_inner < 0
            escaped = false;
            t_arrival = NaN;
            return;
        end
        s_exit = s_outer - s_inner;
        if s_exit < 0 || isnan(s_outer)
            escaped = false;
            t_arrival = NaN;
            return;
        end
    else
        s_exit = s_outer;
        if isnan(s_outer) || s_outer < 0
            escaped = false;
            t_arrival = NaN;
            return;
        end
    end

    % 計算光學厚度
    tau = params.kappa * params.density * s_exit;

    if rand() < exp(-tau)
        escaped = true;
        c = 1;
        if r < params.R_inner
            t_arrival = photon.t_emit + (s_inner + s_exit) / c;
            escaped_pos = photon.pos + (s_inner + s_exit) * photon.dir;
        else
            t_arrival = photon.t_emit + (s_exit) / c;
            escaped_pos = photon.pos + s_exit * photon.dir;
        end
        t_thermal_obs = t_arrival;
    else
        escaped = false;
        s_absorb = rand() * s_exit;
        if r < params.R_inner
            s_total = s_inner + s_absorb;
            absorb_point = photon.pos + s_total * photon.dir;
        else
            s_total = s_absorb;
            absorb_point = photon.pos + s_total * photon.dir;
        end

        c = 1;
        t_arrival = photon.t_emit + s_total / c;

        % 計算 thermal photons 被觀測到的時間
        if ~params.thermal_reemit_dir
            s_remain = s_exit - s_absorb;
            t_thermal_obs = photon.t_emit + (s_total + s_remain) / c + params.thermal_reemit_delay;

        else
            new_dir_angle = rand() * 2*pi;
            new_dir = [cos(new_dir_angle), sin(new_dir_angle)];
            [~, s_out_new] = compute_path_lengths(absorb_point, new_dir, params);
            if isnan(s_out_new) || s_out_new < 0
                s_remain = max(s_exit - s_absorb, 0);
            else
                s_remain = s_out_new;
            end
            t_thermal_obs = photon.t_emit + s_total / c + s_remain / c + params.thermal_reemit_delay;
        end

        % 吸收點位置
        escaped_pos = absorb_point;
    end
end

function [s_inner, s_outer] = compute_path_lengths(pos, dir, params)

    dir = dir / norm(dir);
    if strcmp(params.dimension, '2D')
        a = dot(dir, dir);
        b = 2 * dot(pos, dir);
        
        % 到內邊界的距離
        c_inner = dot(pos, pos) - params.R_inner^2;
        discriminant_inner = b^2 - 4*a*c_inner;
        
        if discriminant_inner >= 0
            s_inner = (-b + sqrt(discriminant_inner)) / (2*a);
            if s_inner < 0
                s_inner = NaN;
            end
        else
            s_inner = NaN;
        end
        
        % 到外邊界的距離
        c_outer = dot(pos, pos) - params.R_outer^2;
        discriminant_outer = b^2 - 4*a*c_outer;
        
        if discriminant_outer >= 0
            s_outer = (-b + sqrt(discriminant_outer)) / (2*a);
            if s_outer < 0
                s_outer = NaN;
            end
        else
            s_outer = NaN;
        end
        
    else
        b = dot(pos, dir);

        % 計算到外球面的距離
        c_outer = dot(pos, pos) - params.R_outer^2;
        discriminant_outer = b^2 - c_outer;
        
        if discriminant_outer >= 0
            s_outer = -b + sqrt(discriminant_outer);
            if s_outer < 0
                s_outer = NaN;
            end
        else
            s_outer = NaN;
        end
        
        % 計算到內球面的距離
        c_inner = dot(pos, pos) - params.R_inner^2;
        discriminant_inner = b^2 - c_inner;
        
        if discriminant_inner >= 0
            s_inner = -b + sqrt(discriminant_inner);
            if s_inner < 0
                s_inner = NaN;
            end
        else
            s_inner = NaN;
        end
    end
end

function lightcurve = record_lightcurve(t_list, t_max, dt)
    % 記錄光變曲線
    t_list = round(t_list, 6);
    edges = 0:dt:t_max;
    lightcurve = histcounts(t_list, edges);
end

% function plot_lightcurves(t_records, t_absorbed_records, N_th, params) % 歸一化長條圖
%     % 繪製光變曲線
%     figure;
% 
%     % 第一個子圖: Non-thermal light curve
%     subplot(3, 1, 1);
%     lc = record_lightcurve(t_records, params.T_max, params.dt);
%     if max(lc) > 0 % 避免除以零
%         lc_normalized = lc / max(lc); % 最大值歸一化
%     else
%         lc_normalized = lc; % 如果數據全為零，保持原樣
%     end
%     t_axis = 0:params.dt:params.T_max-params.dt;
%     bar(t_axis, lc_normalized, 'FaceColor', 'b', 'BarWidth', 1);
%     title('Normalized Non-Thermal Light Curve');
%     xlabel('Time'); ylabel('Photon Counts');
%     grid on; ylim([0, 1]);
% 
%     % 第二個子圖: Absorbed light curve
%     subplot(3, 1, 2);
%     lc_absorbed = record_lightcurve(t_absorbed_records, params.T_max, params.dt);
%     if max(lc_absorbed) > 0 % 避免除以零
%         lc_absorbed_normalized = lc_absorbed / max(lc_absorbed); % 最大值歸一化
%     else
%         lc_absorbed_normalized = lc_absorbed; % 如果數據全為零，保持原樣
%     end
%     bar(t_axis, lc_absorbed_normalized, 'FaceColor', 'k', 'BarWidth', 1);
%     title('Normalized Absorbed Light Curve');
%     xlabel('Time'); ylabel('Photon Counts');
%     grid on; ylim([0, 1]);
% 
%     % 第三個子圖: Thermal light curve
%     subplot(3, 1, 3);
%     if max(N_th) > 0 % 避免除以零
%         N_th_normalized = N_th / max(N_th); % 最大值歸一化
%     else
%         N_th_normalized = N_th; % 如果數據全為零，保持原樣
%     end
%     bar(t_axis, N_th_normalized, 'FaceColor', 'r', 'BarWidth', 1);
%     title('Normalized Thermal Light Curve');
%     xlabel('Time'); ylabel('Photon Counts');
%     grid on; ylim([0, 1]);
% end

function plot_lightcurves(t_records, t_absorbed_records, N_th, params) % 歸一化曲線圖
    % 繪製光變曲線
    figure;

    % 第一個子圖: Non-thermal light curve
    subplot(3, 1, 1);
    lc = record_lightcurve(t_records, params.T_max, params.dt);
    if max(lc) > 0 % 避免除以零
        lc_normalized = lc / max(lc); % 最大值歸一化
    else
        lc_normalized = lc; % 如果數據全為零，保持原樣
    end
    plot(0:params.dt:params.T_max-params.dt, lc_normalized, 'color', 'b', 'LineWidth', 1);
    title('Normalized Non-Thermal Light Curve');
    xlabel('Time'); ylabel('Photon Counts');

    % 第二個子圖: Absorbed light curve
    subplot(3, 1, 2);
    lc_absorbed = record_lightcurve(t_absorbed_records, params.T_max, params.dt);
    if max(lc_absorbed) > 0 % 避免除以零
        lc_absorbed_normalized = lc_absorbed / max(lc_absorbed); % 最大值歸一化
    else
        lc_absorbed_normalized = lc_absorbed; % 如果數據全為零，保持原樣
    end
    plot(0:params.dt:params.T_max-params.dt, lc_absorbed_normalized, 'color', 'k', 'LineWidth', 1);
    title('Normalized Absorbed Light Curve');
    xlabel('Time'); ylabel('Photon Counts');

    % 第三個子圖: Thermal light curve
    subplot(3, 1, 3);
    if max(N_th) > 0 % 避免除以零
        N_th_normalized = N_th / max(N_th); % 最大值歸一化
    else
        N_th_normalized = N_th; % 如果數據全為零，保持原樣
    end
    t_axis = 0:params.dt:params.T_max-params.dt;
    plot(t_axis, N_th_normalized, 'color', 'r', 'LineWidth', 1);
    title('Normalized Thermal Light Curve');
    xlabel('Time'); ylabel('Photon Counts');
end
