% script to test the quadratic spline interpolation function
% on a chosen function

clear; format compact
% interval of interpolation
a = -1; b = 1;
% number of times to interpolate
ntimes = 5;
for nn = 1:ntimes
    n = 2^(nn+1); nint(nn) = n;
    somefun = 'Arc';				% choose a function
    x = (linspace(a, b, n+1))';			% equidistant knots
    y = feval(somefun, x);
    t(1, 1) = x(1); t(n+2, 1) = x(n+1);		% define data points
    t(2:n+1, 1) = (x(1:n) + x(2:n+1))/2;	% for quadratic splines
    s = feval(somefun, t);
    xv = (linspace(a, b, 200))';		% 200 evaluation points
    yv = feval(somefun, xv);
    yv2= qspline(x, s, xv);			% construct and evaluate
    yv3=  spline(x, y, xv);			% two pp interpolants
    yx2= qspline(x, s, x);
    yt3=  spline(x, y, t);
    if n == 16
        subplot(2, 2, 1)
        plot(xv, yv, 'k:', xv, yv2, 'r-', xv, yv3, 'b--');
        xlabel('x'); ylabel('y');
        title(['pp interp. of the ' somefun ' funct. based on ' ...
               int2str(n+1) ' knots'])
        legend('exact', 'quadratic', 'cubic', 2)
        subplot(2, 2, 2)
        plot(xv, yv2-yv, 'r-', xv, yv3-yv, 'b--');
        xlabel('x'); ylabel('error');
        title(['pp interp. of the ' somefun ' funct. based on ' ...
               int2str(n+1) ' knots'])
        legend('quadratic', 'cubic', 2)
        subplot(2, 2, 3)
        plot(x, yx2-y, 'r-', t, yt3-s, 'b--');
        xlabel('x'); ylabel('error');
        title(['pp interp. of the ' somefun ' funct. based on ' ...
               int2str(n+1) ' knots'])
        legend('quadr@knots', 'cubic@midpts', 2)
    end
    err2(nn) = max(abs(yv2-yv));		% measure the max error
    err3(nn) = max(abs(yv3-yv));
    err2x(nn) = max(abs(yx2-y));
    err3t(nn) = max(abs(yt3-s));
end
% format the output
fprintf('pp interp. errors and converg. orders for the %s function\n', somefun);
fprintf('no. of     max error       max error       on knots        on midpts\n');
fprintf('grid       quadratic       cubic           quadratic       cubic\n');
fprintf('intervals  splines         splines         splines         splines\n');
nn = 1;
fprintf('%6d  %9.2e %15.2e %15.2e %15.2e\n', ...
        nint(nn), err2(nn), err3(nn), err2x(nn), err3t(nn));
for nn = 2:ntimes
    ord2(nn)  = log(err2(nn-1) /err2(nn)) /log(nint(nn)/nint(nn-1));
    ord3(nn)  = log(err3(nn-1) /err3(nn)) /log(nint(nn)/nint(nn-1));
    ord2x(nn) = log(err2x(nn-1)/err2x(nn))/log(nint(nn)/nint(nn-1));
    ord3t(nn) = log(err3t(nn-1)/err3t(nn))/log(nint(nn)/nint(nn-1));
    fprintf('%6d  %9.2e %5.2f %9.2e %5.2f %9.2e %5.2f %9.2e %5.2f\n', ...
        nint(nn), err2(nn), ord2(nn), err3(nn), ord3(nn), ...
                  err2x(nn), ord2x(nn), err3t(nn), ord3t(nn));
%   fprintf('%6d  %9.2e %15.2e %15.2e %15.2e\n', ...
%       nint(nn), err2(nn), err3(nn), err2x(nn), err3t(nn));
end
% plot error vs n
subplot(2, 2, 4)
loglog(nint, err2, 'r-', nint, err3, 'b--', ...
       nint, err2x, 'r-o', nint, err3t, 'b--*')
xlabel('number of intervals'); ylabel('error');
title(['convergence of pp interp. of the ' somefun ' funct.'])
legend('quadratic', 'cubic', 'quadr@knots', 'cubic@midpts', 3)
axis tight; hax = gca; set(hax, 'XTick', nint);
%print -dpsc -append testQ.m.ps
