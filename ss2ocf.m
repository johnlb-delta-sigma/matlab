%% ss2ocf: transform state-space system into observable canonical form.
function [sysout,T] = ss2ocf(varargin)
	if nargin==1
		sysin = varargin{1};
		[A B C D] = ssdata(sysin);
	elseif nargin==4
		%% CT system
		sysin = ss(varargin{1},varargin{2},varargin{3},varargin{4});
		A = varargin{1};
		B = varargin{2};
		C = varargin{3};
		D = varargin{4};
	elseif nargin==5
		%% DT system
		sysin = ss(1,1,1,1,varargin{5});
		A = varargin{1};
		B = varargin{2};
		C = varargin{3};
		D = varargin{4};
	else 
		error('Wrong # of arguments')
	end


	n = size(A,1);
	m = size(B,2);
	p = size(C,2);
	O = obsv(A,C);
	Ao = O*A*inv(O);


	%% This transform exposes the numerator coefficients on B
	W = eye(n);
	for ii = 1:n-1
		W = W + diag(-ones(1,ii)*Ao(end,ii+1),-n+ii);
	end
	
	T = W*O;
	Aocf = T*A*inv(T);
	Bocf = T*B;
	Cocf = C*inv(T);
	Docf = D;


	%% Easiest way to transfer all other properties to output
	sysout = ss(sysin);
	sysout.A = Aocf;
	sysout.B = Bocf;
	sysout.C = Cocf;
	sysout.D = Docf;
