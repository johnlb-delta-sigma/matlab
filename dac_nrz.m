%% dac_nrz: Non-return-to-Zero DAC's impulse response.
function [out] = dac_nrz(t, Ts)
	out = dac_rz(t, 0,Ts);
end