%% c2d_general: Convert ct-system into equivalent dt-system 
%% 					based on any DAC shape.
%%
%% 				`dac_fn` should be function that takes one
%% 					scalar input (time) and returns a
%% 					scalar output.
%%
%% 					The DAC response may span up to 2
%% 					cycles (i.e. over the time interval
%% 					[0 2*Ts])
%%
%%				`c0` is the ELD path gain. 
%% 					c0 = 0 if dac_fn=0 for all t>Ts
function [sys_dt,Bd1,Bd2] = c2d_general(sys_ct,c0,Ts,dac_fn)

	[Ac Bc Cc Dc] = ssdata(sys_ct);

	n = size(Ac,1);
	m = size(Bc,2);
	p = size(Cc,2);


	%% Take care of the easy ones.
	Ad = expm(Ac*Ts);
	Cd = Cc;
	Dd = Dc;


	%% Calculate effects of the DAC
	integrand1 = @(t) expm(Ac*(Ts-t))*Bc*dac_fn(t);
	integrand2 = @(t) expm(Ac*(Ts-t))*Bc*dac_fn(t+Ts);
	integrand3 = @(t) expm(Ac*(Ts-t))*Bc*dac_fn(t+2*Ts);
	Bd1 = integral(integrand1, 0,Ts, 'ArrayValued',true);
	Bd2 = integral(integrand2, 0,Ts, 'ArrayValued',true);
	Bd3 = integral(integrand3, 0,Ts, 'ArrayValued',true);


	Aaug = [[Ad Bd2]; zeros(m,n+m)];
	Baug = [Bd1;eye(m,m)];
	Caug = [Cd c0];
	% Aaug = [[Ad Bd2 Bd3]; zeros(2*m,n+2*m)];
	% Baug = [Bd1; eye(m,m); eye(m,m)];
	% Caug = [Cd c0 0];
	Daug = Dd;


	% sys_dt = ss(Ad,Bd,Cd,Dd, Ts);
	sys_dt = ss(Aaug,Baug,Caug,Daug, Ts);
