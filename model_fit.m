%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% <- Author: Junmin Wang -> %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% read data
full_data = readtable('yfp_as_fn_of_binned_mKate_BFP_tidied_20x20.csv');

%% make matrix
% column1: mKate; column2: EBFP2; column3: gRNA; 
% column4: uORF; column5: EYFP; column6: cellCount
full_data_mat = [full_data{:, :}];

%% partition into training and testing datasets
full_data_mat(:, [1 2 5]) = log10(full_data_mat(:, [1 2 5]));
x_train = full_data_mat((full_data_mat(:, 3) ~= 3 & full_data_mat(:, 3) ~= 12), 1:4);
x_test = full_data_mat((full_data_mat(:, 3) == 3 | full_data_mat(:, 3) == 12), 1:4);
y_train = full_data_mat((full_data_mat(:, 3) ~= 3 & full_data_mat(:, 3) ~= 12), 5);
y_test = full_data_mat((full_data_mat(:, 3) == 3 | full_data_mat(:, 3) == 12), 5);

%% training the model
% X --- 
% column1: mKate; column2: EBFP2;
% column3: gRNA; column4: uORF;
minfn = @(param) sum((y_train - (param(1) + (...
    param(8) * x_train(:, 3) ./ (x_train(:, 3) + param(9)) + ...
    param(10) * x_train(:, 4)) ./ ...
    (1 + exp(param(2) + x_train(:, 1) * param(3) + ...
    x_train(:, 2) * param(4) + x_train(:, 3) * param(5) + ...
    x_train(:, 4) * param(6) + ...
    x_train(:, 1) .* x_train(:, 2) * param(7) )))).^2);
start = [5 0 0 0 0 0 0 2.5 1 -2.5];
ub = [7 10 1 1 1 1 1 5 2 0];    
lb = [3 -10 -1 -1 -1 -1 -1 0 0 -5]; 
opts = optimoptions('fmincon','Algorithm','sqp');
problem = createOptimProblem('fmincon', 'objective', minfn, ...
    'x0', start, 'ub', ub, 'lb', lb, 'options', opts);
gs = GlobalSearch;

params_tbl = zeros(100, length(start));
err_lst = zeros(100, 1);
for i=1:100,
    rng(i, 'twister');
    [params, err] = run(gs, problem);
    params_tbl(i, :) = params;
    err_lst(i) = err;
end

% retrieve minimal error and best parameters
min_err = min(err_lst);
best_params = params_tbl(err_lst == min(err_lst), :);

%% simulate the output using the best parameters
simfn = @(param, x) param(1) + (...
    param(8) * x(:, 3) ./ (x(:, 3) + param(9)) + ...
    param(10) * x(:, 4)) ./ ...
    (1 + exp(param(2) + x(:, 1) * param(3) + ...
    x(:, 2) * param(4) + x(:, 3) * param(5) + ...
    x(:, 4) * param(6) + ...
    x(:, 1) .* x(:, 2) * param(7)));

model_out_train = simfn(best_params, x_train);
model_out_test = simfn(best_params, x_test);

% save parameters
best_params_min_err_tbl = array2table([best_params min_err]);
best_params_min_err_tbl.Properties.VariableNames = {'k0', 'a0', 'a1', ...
    'a2', 'a3', 'a4', ...
    'a12', 'c', 'm', ...
    'd', 'minErr'};
writetable(best_params_min_err_tbl, 'best_params_min_err_tbl.csv');

% save simulated output
train_res_mat = [x_train y_train model_out_train];
test_res_mat = [x_test y_test model_out_test];
train_res_tbl = array2table(train_res_mat);
train_res_tbl.Properties.VariableNames = {'log10mKate', 'log10EBFP2', 'sgRNA', ...
    'uORF', 'log10obsEYFP', 'log10predEYFP'};
test_res_tbl = array2table(test_res_mat);
test_res_tbl.Properties.VariableNames = {'log10mKate', 'log10EBFP2', 'sgRNA', ...
    'uORF', 'log10obsEYFP', 'log10predEYFP'};
writetable(train_res_tbl, 'train_res_tbl.csv');
writetable(test_res_tbl, 'test_res_tbl.csv');

