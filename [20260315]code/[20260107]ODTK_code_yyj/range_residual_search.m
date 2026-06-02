function [best_rho, best_state] = range_residual_search(ra_deg, dec_deg, obs, station_per_obs, t, search_range_km, sta_alt_m)
%RANGE_RESIDUAL_SEARCH 在 rho 网格上搜索，使「匀速直线 + 名义地心距」模型与所有观测最一致。
%
%   obs : 参与拟合的观测在 ra_deg/dec_deg/t 中的下标（如 1:nObs 表示全部）
%   对每条 rho_test：先算各地心单位矢 u_geo(i)，再对中间历元 i=2..n-1 用时间插值系数 k_i，
%   堆叠方程 (1-k_i)*d1*u1 + k_i*dn*un ≈ rho_test*u_i，最小二乘求 d1,dn；残差为 ||G*d-b|| 的 RMSE (km)。
%   n=3 时与原先仅用三点一条中间约束等价。

    min_rho = search_range_km(1);
    max_rho = search_range_km(2);
    step_rho = search_range_km(3);

    obs = obs(:);
    n = numel(obs);
    if n < 3
        error('range_residual_search: 至少需要 3 条观测');
    end

    ra = ra_deg(obs);
    dec = dec_deg(obs);
    t_sub = t(obs);
    codes = station_per_obs(obs);


    [t_sorted, ord] = sort(t_sub(:));

    ra = ra(ord);
    dec = dec(ord);
    codes = codes(ord);

    if isdatetime(t_sorted)
        denom = seconds(t_sorted(end) - t_sorted(1));
        time_frac = @(idx) seconds(t_sorted(idx) - t_sorted(1)) / denom;
    else
        denom = t_sorted(end) - t_sorted(1);
        time_frac = @(idx) (t_sorted(idx) - t_sorted(1)) / denom;
    end

    if denom <= 0 || ~isfinite(denom)
        error('range_residual_search: 观测时间跨度无效');
    end

    best_res = inf;
    best_rho = NaN;
    best_state = [];

    nInt = n - 2;

    for rho_test = min_rho:step_rho:max_rho
        u_geo = zeros(n, 3);
        ok = true;
        for i = 1:n
            code = codes{i};
            if iscell(code)
                code = code{1};
            end
            code = upper(strtrim(char(string(code))));
            [Lon, Lat] = GetStationCoordinates(code);
            if isnan(Lon) || isnan(Lat)
                warning('range_residual_search: 站点 %s 无坐标', code);
                ok = false;
                break;
            end
            ti = t_sorted(i);
            u_geo(i, :) = topo_j2000_radec_to_geocentric_unit( ...
                ra(i), dec(i), Lat, Lon, sta_alt_m, ti, rho_test).';
        end
        if ~ok
            continue;
        end

        u1 = u_geo(1, :).';
        un = u_geo(n, :).';

        G = zeros(3 * nInt, 2);
        bb = zeros(3 * nInt, 1);
        row = 1;
        for i = 2:n - 1
            ki = time_frac(i);
            G(row:row + 2, :) = [(1 - ki) * u1, ki * un];
            bb(row:row + 2) = rho_test * u_geo(i, :).';
            row = row + 3;
        end

        d = G \ bb;
        if any(d <= 0) || any(~isfinite(d))
            continue;
        end

        res_vec = G * d - bb;
        res = sqrt(mean(res_vec.^2));

        if res < best_res
            best_res = res;
            best_rho = rho_test;
            best_state = d;
        end
    end

    if isnan(best_rho)
        warning('range_residual_search: 搜索范围内无合法解（检查 d1,dn>0 与测站坐标）');
    else
        fprintf('搜索完成，最优 rho = %.1f km，全观测向量残差 RMSE = %.6g km\n', best_rho, best_res);
    end
end
