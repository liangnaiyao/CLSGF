function [x ft] = EProjSimplex_new(v, k)

%
%% Problem
%
%  min  1/2 || x - v||^2  % 给定v求x，一般不给定k
%  s.t. x>=0, 1'x=1
%

if nargin < 2
    k = 1;
end;

ft=1;
n = length(v); % v的长度

v0 = v-mean(v) + k/n; % 和为 1
%vmax = max(v0);
vmin = min(v0);  
if vmin < 0
    f = 1;
    lambda_m = 0;
    while abs(f) > 10^-10
        v1 = v0 - lambda_m;
        posidx = v1>0;
        npos = sum(posidx);
        g = -npos;
        f = sum(v1(posidx)) - k;
        lambda_m = lambda_m - f/g;
        ft=ft+1;
        if ft > 100
            x = max(v1,0);
            break;
        end;
    end;
    x = max(v1,0);

else
    x = v0; % 和为 1的最小值为1直接输出结果
end;