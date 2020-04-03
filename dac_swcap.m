%% dac_swcap: Switched Cap DAC's impulse response.
function [out] = dac_swcap(t, t0,t1, vref,r,tau)
	window = (t>=t0 & t<=t1);
	out = ( vref*exp(-t/tau) ).*window;
end