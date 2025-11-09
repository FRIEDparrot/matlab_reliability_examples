function [Pf, msc] = time_MCS(mu_, sigma_, g, t, num_MCS)
    if (size(t, 1) == 1)
        t_new = t';
    else
        t_new = t;
    end
	xp = lhsnorm(mu_, sigma_, num_MCS , "on");
	If = zeros(num_MCS, 1);
	for i = 1:num_MCS
		yp = g(xp(i,:), t_new);  % 投影到整个t部分
		if ~isempty(find(yp < 0, 1))
			If(i,:) = 1;
		end
		if mod(i, 5000) == 0
			fprintf("epoch: %d\n", i)
		end
	end
	Pf = sum(If) ./ num_MCS;

    % ?
    msc.Pf_var = (Pf - Pf.^2)./(num_MCS -1);
    msc.Pf_cov = sqrt(msc.Pf_var)./Pf;
end
