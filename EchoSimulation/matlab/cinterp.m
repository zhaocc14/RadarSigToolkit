function y2 = cinterp(x1,y1,x2)
    % 对复数向量插值，通过对幅度和相位分别插值实现。
    % interp a complex array by its abs and angle, respectively.
    y1_abs = abs(y1);
    y1_ang = angle(y1);
    y1_ang = unwrap(y1_ang);
    
    y2_abs = interp1(x1, y1_abs, x2,'linear','extrap');
    y2_ang = interp1(x1, y1_ang, x2,'linear','extrap');
    y2 = y2_abs.*exp(1j*y2_ang);
    end