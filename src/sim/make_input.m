function uFun = make_input(K, xeq, ueq)

% Build a state-feedback law u(t,x) = ueq + K*(x - xeq).

% Control law:  u = ueq + K*(x - xeq)
%
% K is designed on the linearised model in deviation coordinates
%   dx = x - xeq,  du = u - ueq,   d(dx)/dt = A*dx + B*du,
% where the LMI returns K such that du = K*dx stabilises it. Writing that
% relation back in absolute variables gives u = ueq + K*(x - xeq).

xeq = xeq(:); ueq = ueq(:);
if isempty(K)
    uFun = @(t, x) ueq;
else
    uFun = @(t, x) ueq + K * (x(:) - xeq);
end

end