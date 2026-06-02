function r_eci_km = geodetic_to_eci_km(lat_deg, lon_deg, alt_m, utc_dt)
%GEODETIC_TO_ECI_KM 测站大地坐标 → 近似地心惯性系位置（km），用于与 MEME/J2000 赤经赤纬配套。
%   采用 WGS84 椭球 + UTC 近似 GMST 的 ECEF→ECI 绕 Z 旋转（与真 J2000 有小偏差，视差修正足够）。
%
%   lat_deg, lon_deg : 大地纬度、经度（度），东经为正
%   alt_m            : 椭球高（m）
%   utc_dt           : datetime，视为 UTC

    a_km = 6378.137;
    f = 1 / 298.257223563;
    e2 = 2 * f - f^2;

    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);
    h_km = alt_m / 1000;

    sinlat = sin(lat);
    N = a_km ./ sqrt(1 - e2 * sinlat.^2);
    x = (N + h_km) .* cos(lat) .* cos(lon);
    y = (N + h_km) .* cos(lat) .* sin(lon);
    z = (N * (1 - e2) + h_km) .* sinlat;
    r_ecef = [x; y; z];

    jd = juliandate(utc_dt);
    d = jd - 2451545.0;
    gmst_deg = mod(280.46061837 + 360.98564736629 * d, 360);
    gmst = deg2rad(gmst_deg);

    c = cos(gmst);
    s = sin(gmst);
    R = [c, s, 0; -s, c, 0; 0, 0, 1];
    r_eci_km = R * r_ecef;
end
