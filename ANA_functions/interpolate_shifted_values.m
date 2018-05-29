
function iC = interpolate_shifted_values(ts, vt)

iC.timestamps = ts;
iC.posx       = interp1(vt.timestamps, vt.posx, ts);
iC.posx2      = interp1(vt.timestamps, vt.posx2, ts);
iC.posy       = interp1(vt.timestamps, vt.posy, ts);
iC.posy2      = interp1(vt.timestamps, vt.posy2, ts);
iC.posx_c     = interp1(vt.timestamps, vt.posx_c, ts);
iC.posy_c     = interp1(vt.timestamps, vt.posy_c, ts);
iC.poshd      = interp1(vt.timestamps, vt.poshd, ts);
iC.vx         = interp1(vt.timestamps, vt.vx, ts);
iC.vy         = interp1(vt.timestamps, vt.vy, ts);
iC.speed      = interp1(vt.timestamps, vt.speed, ts);

end