function g = linearize(xss, u, funct)
% xss is an array containing the states at steady-state
% u is an array containing the inputs at steady-state
% funct is a function handle

% Number of states and inputs
nx = length(xss);
nu = length(u);

% Symbolic states and inputs
X = sym('X',[1 nx]).';
U = sym('U',[1 nu]).';

% Symbolic differential equation
dXdt = funct(0, X, U);

% Linearizing
Apos = jacobian(dXdt,X);
Bpos = jacobian(dXdt,U);

% Continuous state-space model
Apos = double(subs(Apos, [X; U], [xss'; u']));
Bpos = double(subs(Bpos, [X; U], [xss'; u']));

% Transfer function
s = tf('s');
g = (s*eye(size(Apos,1)) - Apos)\(Bpos);
end
