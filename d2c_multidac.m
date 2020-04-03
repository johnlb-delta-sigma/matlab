%% d2c_multidac: Given multiple DACs, with connections restricted to a subspace
%%               defined by Ps, find a system, Gc(s), that is equivalent to 
%%               Gd(z) in the impulse-invariant sense.
%%
%%               Gc(s) is returned as a state-space system in a subspace that
%%               exposes the coefficients used to implement the filter. It
%%               does this by putting the output into the same subspace
%%               as the result of the function ``implementor_fn``, which takes
%%               an arbitrary A matrix and returns the A matrix of the form
%%               the output should be in. Note that Ps acts in this subspace,
%%               so that restricting one DAC to only, say, the inputs to 2 
%%               specific nodes in the filter behaves in the expected way.
function [sys_ct, Ad,Ac,Vd,lambda_d] = d2c_multidac(sys_dt, dac_fns, implementor_fn, Ps)
	Ts = sys_dt.Ts;

	Gd = ss2ocf(sys_dt);
	[Ad Bd Cd Dd] = ssdata(Gd);
	N = size(Ad,1);
	M = length(dac_fns);


	%% Vd:  Ad_ocf --> Ad_blkdiag
	% [Vd lambda_d] = jordan(Ad)
	% % [Vd lambda_d] = eig(Ad);
	% [Vd lambda_d] = cdf2rdf(Vd,lambda_d);

	%% Try schur form instead
	[Vd lambda_d] = schur(Ad);


	%% Desired output form
	Ac = implementor_fn(logm(Ad)/Ts);


	%% Calculate the change-of-basis between DT ocf
	%% and the desired CT form.
	[~, V_od] = ss2ocf(expm(Ac*Ts),zeros(size(Bd)),Cd,0);



	%% DACs
	Phi = @(t) Phi_fn(t,dac_fns,N);




	%% Calculate effects of the DAC

	Ac_ocf_d = logm(Ad)/Ts;
	integrand1 = @(t) expm(Ac_ocf_d*(Ts-t))*Phi(t);
	integrand2 = @(t) expm(Ac_ocf_d*(Ts-t))*Phi(t+Ts);
	Bd1_ = integral(integrand1, 0,Ts, 'ArrayValued',true);	
	Bd2_ = integral(integrand2, 0,Ts, 'ArrayValued',true);	
	Gamma1 = [ Bd1_;
			   zeros(1,2*N); ];
	Gamma2 = [ zeros(1,2*N);
			   Bd2_		     ];

	Gamma  = Gamma1 + Gamma2;


	%% Augment DT system with (z/z) s.t. we can equate.
	Ad_aug = [[Ad [zeros(N-1,1); 1]]; zeros(1,N+1)];
	Bd_aug = [Bd; 0];
	Cd_aug = [Cd 0];
	Dd_aug = Dd;



	Vstar_odes = blkdiag(V_od, V_od);
	Gamma_des = Gamma*Vstar_odes;

	Bc_star = Ps*pinv(Gamma_des*Ps)*Bd_aug;


	%% Check that the inversion is consistent
	%% (i.e. that a solution exists)
	Bd_err = [Bd;0] - Gamma*Vstar_odes*Bc_star;
	errnorm = norm(Bd_err);
	if errnorm > sqrt(eps)
		warning(['Failed to find exact solution. Norm of error is: ' ...
			      num2str(errnorm)]);
	end
		

	%% Adjust for desired CT subspace:  
	Bc = [];
	Dc = [];
	for ii=0:M-1
		Bc = [Bc Bc_star(ii*N+1:(ii+1)*N)];
		Dc = [Dc Dd];
	end
	Cc = Cd*inv(V_od);


	%% Output state-space system
	sys_ct = ss(real(Ac),real(Bc),real(Cc),real(Dc));

end



%% Phi_fn: calculate the function Phi(t) from provided dac_fns
%%         n is order of system.
function [x] = Phi_fn(t,dac_fns,n)
	x = [];
	for ii = 1:length(dac_fns)
		x = [x eye(n)*dac_fns{ii}(t)];
	end
end


%% fix_lambda_ordering: reorder Vc s.t. logm(lambda_d)/Ts == lambda_c
function [Vc lambda_c] = fix_lambda_ordering(lambda_d,Ac,Vc,Ts)

	tol = sqrt(eps);
	N = length(Vc);
	lambda_dc = logm(lambda_d)/Ts;

	ptrs = perms(1:N);
	success = false;
	minerr = Inf;
	%% Just brute force for lack of a better way.
	for ii=1:size(ptrs,1)
		thisptr = ptrs(ii,:);
		thisVc  = Vc(:,thisptr);

		thislambda = inv(thisVc)*Ac*thisVc;
		err = norm(norm(thislambda-lambda_dc));
		if err<minerr
			minerr=err;
			minVc =thisVc;
		end

		if err < tol
			success = true;
			break;
		end
	end

	if (~success)
		warning(['Failed to find correct ordering in d2c_multidac. ' ...
				 'Norm of minimum error = ' num2str(minerr)]);
	end
	Vc = minVc;
	lambda_c = thislambda;

end