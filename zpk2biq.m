%% zpk2biq: transforms a zpk filter into a set of biquads
function [biqs] = zpk2biq(X)

	z = X.z{1};
	p = X.p{1};
	k = X.k;


	zelts = construct_rt_elts(z);
	pelts = construct_rt_elts(p);

	% assuming the system is causal, if there are more
	% zero elements than pole elements, we will need to
	% pair 2 real zeros with one pole pair.
	realpair_flag = length(zelts)>length(pelts);

	% sort poles by highest freq
	for ii=1:length(pelts)
		pimag(ii) = imag(pelts{ii}(1));
	end
	parray = flipud(sortrows([pimag' (1:length(pelts))']));
	pelts = {pelts{parray(:,2)}};


	% pair up elements by distance
	for ii = 1:min(length(zelts),length(pelts))
		zarray = [];
		for jj = 1:length(zelts)
			zarray(jj) = zelts{jj}(1);
		end
		
		d{ii} = abs( zarray.' - pelts{ii}(1) );
		sorted{ii} = sortrows([d{ii} (1:length(zelts))' imag(zarray.')]);
		

		choice_index = sorted{ii}(1,2);
		pz{ii} = {pelts{ii} zelts{choice_index}};
		if (imag(zelts{choice_index})==0 & realpair_flag)
			sorted_index = find((sorted{ii}(:,3)==0) & (sorted{ii}(:,2)~=choice_index));
			choice_index_ = sorted{ii}(sorted_index,2);
			pz{ii} = {pelts{ii} [zelts{choice_index} zelts{choice_index_}]};

		end

		zelts(choice_index) = [];
		realpair_flag = length(zelts)>length(pelts)-ii;
	end

	% any remaining poles left over will be synthesized by themselves.
	for ii = [ii+1:length(pelts)]
		% zarray = [];
		% for jj = 1:length(zelts)
		% 	zarray(jj) = zelts{jj}(1);
		% end
		
		% d{ii} = abs( zarray.' - pelts{ii}(1) );
		% sorted{ii} = sortrows([d{ii} (1:length(zelts))' imag(zarray.')]);
		

		% choice_index = sorted{ii}(1,2);
		% pz{ii} = {pelts{ii} zelts{choice_index}};
		% if (imag(zelts{choice_index})==0 & realpair_flag)
		% 	sorted_index = find((sorted{ii}(:,3)==0) & (sorted{ii}(:,2)~=choice_index))
		% 	choice_index_ = sorted{ii}(sorted_index,2);
		% 	pz{ii} = {pelts{ii} [zelts{choice_index} zelts{choice_index_}]};
		% end
		
		pz{ii} = {pelts(ii) []};
	end



	% synthesize all biquads from poles/zeros.
	for ii = 1:length(pz)
		p_ = pz{ii}{1};
		z_ = pz{ii}{2};
		if ii==1
			k_ = k;
		else
			k_ = 1;
		end

		biqs{ii} = zpk(z_,p_,k_);
	end


	% reorder biqs so lowest Q filter is first.
	Qs = calc_Qs(biqs)'
	sorted = sortrows([Qs (1:length(Qs))']);
	biqs = {biqs{sorted(:,2)}};

end


%% construct_rt_elts: construct datastructure to hold either a single real or piar of cconj roots.
function [rt_elts] = construct_rt_elts(rts)

	[real_index] = find(imag(rts)==0);
	[imag_index] = find(imag(rts)>0);
	for ii = 1:length(imag_index)
		conj_index(ii,1) = find(imag(rts)==-imag(rts(imag_index(ii))));
	end

	rt_elts = {};
	for ii = 1:(length(real_index))
		rt_elts{ii} = [rts(real_index(ii))];
	end
	if length(ii)<1
		ii=0;
	end
	for jj = 1:(length(imag_index))
		rt_elts{ii+jj} = [rts(imag_index(jj)) rts(conj_index(jj))];
	end


end


function [Q] = calc_Qs(biqs)

	for ii = 1:length(biqs)
		G = tf(biqs{ii});

		if (length(G.den{1}) < 3)
			Q(ii) = 0;
		else
			Q(ii) = abs( sqrt(G.den{1}(3))/G.den{1}(2) );
		end
	end

end



% %% construct_rt_elts: 
% function [rt_elts] = construct_rt_elts(rts)

% 	[real_index] = find(imag(rts)==0);
% 	[imag_index] = find(imag(rts)>0);
% 	for ii = 1:length(imag_index)
% 		conj_index(ii,1) = find(imag(rts)==-imag(rts(imag_index(ii))));
% 	end

% 	rt_elts = {};
% 	for ii = 1:(length(real_index))
% 		rt_elts{ii} = [rts(real_index(ii))];
% 	end
% 	if length(ii)<1
% 		ii=0;
% 	end
% 	for jj = 1:(length(imag_index))
% 		rt_elts{ii+jj} = [rts(imag_index(jj)) rts(conj_index(jj))];
% 	end


% end