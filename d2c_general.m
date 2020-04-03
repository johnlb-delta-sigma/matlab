%% d2c_general: Convert ct-system into equivalent dt-system 
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
%% 					`dac_fn` must return a scalar for
%% 					scalar input
%%
%%				`c0` is the ELD path gain. 
%% 					c0 = 0 if dac_fn=0 for all t>Ts
function [sys_ct, c0, phi] = d2c_general(sys_dt,dac_fn)

	%% Check that phi will be invertible
	if find( pole(sys_dt)<=0 & imag(pole(sys_dt))==0 )
		warning(['System has poles either on the negative real axis or ' ...
				'at the origin. D2C conversion is not well-defined.']);
	end


	%% Put into observable form first.
	sys_dt = ss2ocf(sys_dt);
	[Ad Bd Cd Dd] = ssdata(sys_dt);
	
	Ts = sys_dt.Ts;

	n = size(Ad,1);
	m = size(Bd,2);
	p = size(Cd,2);


	%% Take care of the easy ones.
	Ac = logm(Ad)/Ts;
	Dc = Dd;
	Cc = Cd;

	%% Calculate effects of the DAC
	integrand1 = @(t) expm(Ac*(Ts-t))*dac_fn(t);
	integrand2 = @(t) expm(Ac*(Ts-t))*dac_fn(t+Ts);
	Bd1_ = integral(integrand1, 0,Ts, 'ArrayValued',true);	
	Bd2_ = integral(integrand2, 0,Ts, 'ArrayValued',true);	


	%% Solve for c0 & Bc
	phi1 = [ Bd1_
			 zeros(1,n); ];
	phi2 = [ zeros(1,n);
			 Bd2_		 ];

	den_coeff = [1;-Ad(:,1)];
	phi = [den_coeff (phi1+phi2)];
	v = inv(phi)*[Bd; zeros(1,m)];

	c0 = v(1,:);
	Bc = v(2:end,:);

	%% Make final CT system object.
	sys_ct = ss(Ac,Bc,Cc,Dc);
end