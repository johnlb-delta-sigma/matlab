%% cheb_ntf: Create a Chebyshev NTF from specifications
%
% 	n: 			order
% 	rejection: 	difference b/w passband and stopband, in dB
% 	ftype: 		Filter type. Can be 'lowpass', 'highpass', or 'bandreject'
% 	fc: 		for 'lowpass' or 'highpass': fc = BW of NTF stopband
%				for 'bandreject': fc = [centerFreq, BW] of NTF stopband
% 	fs = sample rate of modulator
%
% 	Returns [zd pd kd] = [zeros poles NTFgain] in z-domain
function [zd, pd, k] = cheb_ntf(n,rejection,ftype,fc,fs)


	% Create Inverse Chebyshev Prototype
	% alpha = 10^(-(rejection+4.3)/10);
	% delta = 10^(-(rejection+4.3+1)/10);
	% alpha2 = alpha-delta;
	% e = sqrt(alpha2/(1-alpha2));

	alpha = 10^(-(rejection)/10);
	delta = 0;
	e = sqrt(alpha/(1-alpha));


	wc = 2*pi.*fc;
	wc_fs = wc./fs;

	pgen = @(a,b) 1/(-sinh(a)*sin(b) + j*cosh(a)*cos(b));
	% zgen = @(a,b) j*sec(b);
	zgen = @(a,b,c) j*sec(b+j*c);


	% Create frequency transform maps
	% NOTE: Both HP and LP have same xforms.
	% 		This is because a LP NTF can only
	% 		be realized by a HP NTF mirrored
	% 		about IMAG axis. This is done later.
	if 	( strcmp(ftype,'bandreject') )

		freqmap = @(rt,wc) [ 	(wc(2) + sqrt( wc(2)^2 - 4*rt^2*wc(1)^2 ))/(2*rt);
								(wc(2) - sqrt( wc(2)^2 - 4*rt^2*wc(1)^2 ))/(2*rt) 	];

	else 

		freqmap = @(rt,wc) wc/rt;

	end



	% Generate p,z,k
	
	% s-domain + required frequency transform
	a = (1/n)*asinh(1/e);
	k = 1;
	p = [];
	z = [];
	for i=1:n
		b = (2*i-1)*pi/(2*n);
		c = sqrt(delta)/(e*n);
		% c = 0;

		p(i,:) = freqmap( pgen(a,b), wc );
		z(i,:) = freqmap( zgen(a,b,c), wc );

	end

	% make sure all poles/zeros are in a single vector
	% (bandreject and bandpass come out as 2-D arrays)
	rts_shape = size(p);
	p = reshape(p, 1, rts_shape(1)*rts_shape(2));
	z = reshape(z, 1, rts_shape(1)*rts_shape(2));


	% z-domain (bilinear xform)
	% Note: using a higher order exponential
	% 		approximation than is standard.
	% 		We need this for accuracy when
	% 		generating high OSR modulators.
	x=50;
	% pd = (p + 2*fs*x).^x./(2*fs*x - p).^x;
	% zd = (z + 2*fs*x).^x./(2*fs*x - z).^x;
	pd = exp(p/fs);
	zd = exp(z/fs);

	% calculate bilinear xform gain
	kp = (2*fs*x - p);
	kz = (2*fs*x - z);
	k = (prod(kz)/prod(kp))^x;


	% Only practical lowpass NTF is bandpass
	% NTF mirrored about IMAG axis.
	if ( strcmp(ftype,'lowpass') )
		pd = -pd;
		zd = -zd;
	end


	% Make sure TF has purely real coeffs.
	% (i.e. fix rounding errors in p/z locs)
	zd = cplxpair(zd);
	pd = cplxpair(pd);


	if (imag(k)>1e-3)
		error('Imaginary part of Gain is too large in cheb_ntf');
	else
		k = real(k);
	end

end