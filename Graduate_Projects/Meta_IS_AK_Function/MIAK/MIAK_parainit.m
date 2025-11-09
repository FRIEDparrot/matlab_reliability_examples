function params = MIAK_parainit(n)
%% intitialize a default params struct for MIAK_modelbuild function 
%% n: number of the variables 
	
	% params for first model to build and optimize 
	params.num_MKIS_init = 30;   	     % 初始样本点 m0
	params.max_epoch1 = 10;              % 第一次代理模型优化的最大数量
	params.nKMeans = 10;			     % 优化过程每次通过KMeans添加的点数
	params.num_MCMC = 1e4;               % 更新初始代理模型时每次抽取的点
	params.a_LOSS_min = 0.1;             % 初始模型接受阈值a_LOSS下界
	params.a_LOSS_max = 10;			     % 初始模型接收阈值a_LOSS上界

	% params for submodel build (dace.dacefit)
	params.regr = @regpoly0;             % regr function 
	params.corr = @corrgauss; 		     % corr function 
	params.theta0 = ones(1,n);		     % theta definition 
	params.lob  =  1e-5 .* ones(1,n);
	params.upb  = 100 .* ones(1,n);
	
	% params for optimization of submodel in iteration process 
	params.num_ISAK   = 2e3;         % ********* 注意: 这个会影响收敛性, 所以不要调特别高*****
	params.max_epoch2 = 1000;        % ****** 从ISAK样本中最大添加样本数量(< num_ISAK) ****** 
end
