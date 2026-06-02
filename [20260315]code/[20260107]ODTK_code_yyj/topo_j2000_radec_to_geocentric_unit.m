function u_geo = topo_j2000_radec_to_geocentric_unit(ra_deg, dec_deg, lat_deg, lon_deg, alt_m, utc_dt, rho_nom_km)
%TOPO_J2000_RADEC_TO_GEOCENTRIC_UNIT 测站 J2000 赤经赤纬 → 地心 J2000 单位视线矢量。
%
%   测站 J2000：原点在地心 J2000 中的测站位置，坐标轴与地心 MEME/J2000 平行；RA/Dec 给出
%   从测站指向目标的单位矢量 u_topo（与 ODTK 光学站心 MEME J2000 一致）。
%
%   地心 J2000：从地心指向目标的单位矢量 u_geo = (R + d*u_topo) / ||R + d*u_topo||，
%   其中 R 为测站地心矢径（本函数内由经纬度+UTC 近似到与 geodetic_to_eci_km 一致的 ECI），
%   d 为沿视线使 ||R + d*u_topo|| = rho_nom_km 的斜距（射线与地心球面的近交，取最小正根）。
%
%   rho_nom_km : 目标名义地心距（km），如 GEO≈42164、LEO≈6778。未知距离时可在主程序用
%   IOD/LS 得到的 |r| 迭代更新 rho（见 main_try 中 rho 迭代段）。

    u_topo = [cosd(dec_deg) * cosd(ra_deg); cosd(dec_deg) * sind(ra_deg); sind(dec_deg)];
    R = geodetic_to_eci_km(lat_deg, lon_deg, alt_m, utc_dt);

    b = dot(u_topo, R);
    c = dot(R, R) - rho_nom_km^2;
    disc = b * b - c;

    if disc < 0
        d_slant = max(rho_nom_km - norm(R), 1000);
    else
        s = sqrt(disc);
        t1 = -b - s;
        t2 = -b + s;
        cand = [t1; t2];
        pos = cand(cand > 1e-3);
        if isempty(pos)
            d_slant = max(rho_nom_km - norm(R), 1000);
        else
            d_slant = min(pos);
        end
    end

    r = R + d_slant * u_topo;
    nr = norm(r);
    if nr < 1e-12
        u_geo = u_topo;
    else
        u_geo = r / nr;
    end
end
