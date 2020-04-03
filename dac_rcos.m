%% dac_rcos: Switched Cap DAC's impulse response.
function [out] = dac_rcos(t, t0,t1)
	window = (t>=t0 & t<=t1);
	f = pi/(t1-t0);
	out = ( cos(f*(t-t0) - pi/2) ).*window;
end