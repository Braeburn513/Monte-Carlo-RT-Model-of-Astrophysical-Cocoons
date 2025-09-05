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

    % 初始化光變曲線記錄
    t_records = [];
    t_absorbed_records = [];

    % 模擬光子
    for i = 1:params.N_photon
        photon = generate_photon(params);
        [escaped, t_arrival, escaped_pos, t_thermal_obs] = propagate_photon(photon, params);
        
        if escaped % 處理逃逸的光子
            all_escaped_times = [all_escaped_times; t_arrival];
            angle = atan2(escaped_pos(2), escaped_pos(1));
            if angle >= params.arc_angle_range(1) && angle <= params.arc_angle_range(2)
                escaped_initial_positions = [escaped_initial_positions; photon.pos];
                escaped_times = [escaped_times; t_arrival];
                arc_photons_t = [arc_photons_t; t_arrival];
            end
            
            % 記錄光變曲線的程式碼
            t_records(end+1) = t_arrival;

        else % 處理被吸收的光子
            angle = atan2(photon.dir(2), photon.dir(1));
            if angle >= params.arc_angle_range(1) && angle <= params.arc_angle_range(2)
                absorbed_pos = [absorbed_pos; escaped_pos];  % escaped_pos 於被吸收時為吸收點
                absorbed_times = [absorbed_times; t_arrival]; % 真正被吸收的時間（可視化用）
                absorbed_initial_positions = [absorbed_initial_positions; photon.pos];
            end

            % 將「熱光子被觀測到的時間」加入熱光變曲線記錄（不是吸收瞬間）
            t_absorbed_records(end+1) = t_thermal_obs;
        end
    end

    % 統計圓弧光子的 t<1 和 t≥1 數量
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

    % 統計所有光子的 t<1 和 t≥1 數量
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
    if params.plot_thermal
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
    else
        N_th = [];
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
    params.N_photon = 50000;            % 總模擬光子數
    params.dt = 0.01;                   % 時間 bin 寬度
    params.T_max = 2.0;                 % 最大觀測時間
    params.injection_mode = 'point';     % 注入模式: 'point', 'sphere', 'ring'
    params.flare_duration = 0;          % 耀班持續時間
    params.arc_angle_range = [-pi, pi];     % 觀察的圓弧角度範圍 -pi to pi
    params.plot_initial_positions = false;   % 是否繪製光子初始位置可視化

    params.cooling_time = 0.1;   % cooling time
    params.cooling_on = false;    % true有冷卻，false無冷卻
    params.plot_thermal = true;  % 是否畫 thermal light curves

    params.thermal_reemit_delay = 0;   % 再發射的額外延遲
    params.thermal_isotropic = false;  % true再發射方向各向同性，false保持原來方向
end

function photon = generate_photon(params)
    % 生成光子，支援不同注入區域
    switch params.injection_mode
        case 'point'
            r = 0; % 從中心注入
        case 'sphere'
            r = sqrt(rand()) * params.R_inj;
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
    
    theta = rand() * 2 * pi;
    x = r * cos(theta);
    y = r * sin(theta);
    
    angle = rand() * 2 * pi;
    dx = cos(angle);
    dy = sin(angle);
    
    % 所有光子頻率都設為 1
    photon.pos = [x, y];
    photon.dir = [dx, dy];
    photon.t_emit = rand() * params.flare_duration;
end

function [escaped, t_arrival, escaped_pos, t_thermal_obs] = propagate_photon(photon, params)
    % simulate photon transport; return also the time when the thermal
    % re-emitted photon would be observed (t_thermal_obs).
    pos = photon.pos;
    dir = photon.dir;

    % 計算到內外邊界的路徑長
    [s_inner, s_outer] = compute_path_lengths(pos, dir, params);

    % 預設值
    t_thermal_obs = NaN;
    escaped_pos = [NaN, NaN];

    % 決定 photon 在 cocoon 中的 exit 距離 s_exit
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

    % 計算光學厚度（僅在環狀介質中）
    tau = params.kappa * params.density * s_exit;

    if rand() < exp(-tau)
        % photon 真正逃逸
        escaped = true;
        c = 1;
        if r < params.R_inner
            t_arrival = photon.t_emit + (s_inner + s_exit) / c;
            escaped_pos = photon.pos + (s_inner + s_exit) * photon.dir;
        else
            t_arrival = photon.t_emit + (s_exit) / c;
            escaped_pos = photon.pos + s_exit * photon.dir;
        end
        % 對於未被吸收的 photon，熱觀測時間同逃逸時間（若你想）
        t_thermal_obs = t_arrival;
    else
        % photon 在介質中被吸收，選一個 s_absorb
        escaped = false;
        s_absorb = rand() * s_exit;  % 在介質中的吸收距離
        if r < params.R_inner
            s_total = s_inner + s_absorb;  % 總飛行距離到吸收點
            absorb_point = photon.pos + s_total * photon.dir;
        else
            s_total = s_absorb;
            absorb_point = photon.pos + s_total * photon.dir;
        end

        c = 1;
        t_arrival = photon.t_emit + s_total / c; % 真正被吸收的時刻（吸收瞬間）

        % ------ 計算熱光子被觀測到的時間 ------
        if ~params.thermal_isotropic
            % 保持原方向再發射 -> 直接計算剩餘距離到外邊界
            s_remain = s_exit - s_absorb; % remaining distance inside cocoon along the same dir
            t_thermal_obs = photon.t_emit + (s_total + s_remain) / c + params.thermal_reemit_delay;
            % 注意：代數上 (s_total + s_remain) = (s_inner + s_exit) for r<R_inner，
            % 相當於 photon 若未被吸收本來會到達的時間（但我們把吸收時間和熱再發射間的延遲分開）
        else
            % 各向同性再發射 — 在吸收點抽一個新方向，從該點算到外邊界
            new_dir_angle = rand() * 2*pi;
            new_dir = [cos(new_dir_angle), sin(new_dir_angle)];
            % 使用 compute_path_lengths 從吸收點往 new_dir 算出到外邊界距離
            [~, s_out_new] = compute_path_lengths(absorb_point, new_dir, params);
            % 對於吸收點位於內圈或外圈，s_out_new 應該直接是我們要的剩餘距離
            if isnan(s_out_new) || s_out_new < 0
                % 若沒交點（理論上不該發生），退回保守估計：使用原方向的剩餘距離
                s_remain = max(s_exit - s_absorb, 0);
            else
                s_remain = s_out_new;
            end
            t_thermal_obs = photon.t_emit + s_total / c + s_remain / c + params.thermal_reemit_delay;
        end

        % 吸收點位置（作為被吸收 visual 用）
        escaped_pos = absorb_point;
    end
end

function [s_inner, s_outer] = compute_path_lengths(pos, dir, params)
    % 計算光子到環狀介質內外邊界的路徑長度
    a = dot(dir, dir);
    b = 2 * dot(pos, dir);
    
    % 到內邊界的距離
    c_inner = dot(pos, pos) - params.R_inner^2;
    discriminant_inner = b^2 - 4*a*c_inner;
    
    if discriminant_inner >= 0
        s_inner = (-b + sqrt(discriminant_inner)) / (2*a);
        if s_inner < 0
            s_inner = NaN; % 光子已在環狀介質內或無交點
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
end

function lightcurve = record_lightcurve(t_list, t_max, dt)
    % 記錄光變曲線
    t_list = round(t_list, 6);
    edges = 0:dt:t_max;
    lightcurve = histcounts(t_list, edges);
end

% function plot_lightcurves(t_records, t_absorbed_records, N_th, params) % 長條圖
%     % 繪製光變曲線 - 單一視窗內三個子圖（歸一化長條圖）
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
%     if params.plot_thermal && ~isempty(N_th)
%         subplot(3, 1, 3);
%         if max(N_th) > 0 % 避免除以零
%             N_th_normalized = N_th / max(N_th); % 最大值歸一化
%         else
%             N_th_normalized = N_th; % 如果數據全為零，保持原樣
%         end
%         bar(t_axis, N_th_normalized, 'FaceColor', 'r', 'BarWidth', 1);
%         title('Normalized Thermal Light Curve');
%         xlabel('Time'); ylabel('Photon Counts');
%         grid on; ylim([0, 1]);
%     end
% end

% function plot_lightcurves(t_records, t_absorbed_records, N_th, params) % 曲線圖
%     % 繪製光變曲線 - 單一視窗內三個子圖（歸一化）
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
%     plot(0:params.dt:params.T_max-params.dt, lc_normalized, 'color', 'b', 'LineWidth', 1);
%     title('Normalized Non-Thermal Light Curve');
%     xlabel('Time'); ylabel('Normalized Photon Counts');
% 
%     % 第二個子圖: Absorbed light curve
%     subplot(3, 1, 2);
%     lc_absorbed = record_lightcurve(t_absorbed_records, params.T_max, params.dt);
%     if max(lc_absorbed) > 0 % 避免除以零
%         lc_absorbed_normalized = lc_absorbed / max(lc_absorbed); % 最大值歸一化
%     else
%         lc_absorbed_normalized = lc_absorbed; % 如果數據全為零，保持原樣
%     end
%     plot(0:params.dt:params.T_max-params.dt, lc_absorbed_normalized, 'color', 'k', 'LineWidth', 1);
%     title('Normalized Absorbed Light Curve');
%     xlabel('Time'); ylabel('Normalized Photon Counts');
% 
%     % 第三個子圖: Thermal light curve
%     if params.plot_thermal && ~isempty(N_th)
%         subplot(3, 1, 3);
%         if max(N_th) > 0 % 避免除以零
%             N_th_normalized = N_th / max(N_th); % 最大值歸一化
%         else
%             N_th_normalized = N_th; % 如果數據全為零，保持原樣
%         end
%         t_axis = 0:params.dt:params.T_max-params.dt;
%         plot(t_axis, N_th_normalized, 'color', 'r', 'LineWidth', 1);
%         title('Normalized Thermal Light Curve');
%         xlabel('Time'); ylabel('Normalized Photon Counts');
%     end
% end
