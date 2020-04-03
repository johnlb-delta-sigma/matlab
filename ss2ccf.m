%% ss2ccf: transform state-space system into controlable canonical form.
function [sysout] = ss2ccf(varargin)
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
	Ao = Ctr*A*inv(Ctr);


	%% This transform exposes the numerator coefficients on C
	W = eye(n);
	for ii = 1:n-1
		%% The order of the coefficients is reversed
		%% compared to ocf. Also, it is upper triangular
		%% instead of lower triangular.
		W = W + diag(-ones(1,ii)*Ao(ii+1,end),n-ii);
	end
	
	T = Ctr*W;
	Accf = inv(T)*A*T;
	Bccf = inv(T)*B;
	Cccf = C*T;
	Dccf = D;

	%% Easiest way to transfer all other properties to output
	sysout = ss(sysin);
	sysout.A = Accf;
	sysout.B = Bccf;
	sysout.C = Cccf;
	sysout.D = Dccf;
