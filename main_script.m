function main_script()
%code for running the experiments in the paper
%`Quasi Branch and bound for smooth global optimization'
%

%add subfolders to path
current_place=pwd;
addpath(genpath(current_place));
%% choose paramters
time_max=30;% time limit (in seconds)
tol=10^-8;  %error tolerance

experiment_name='Dixon'; %options are Dixon, Rastrigin, or Rastrigin_constrained
switch experiment_name
    case 'Dixon'
    %Dixon Szego experiment
    %problem can be one of {'branin','camel','gold','shubert','hartman3','shekel5','shekel7','shekel10','hartman6'};
    problem='branin';
    solvers={'2','3','alpha'};
    case 'Rastrigin_constrained'
    %Rastrigin  constrained experiment
    problem='cosine';
    solvers={'2lb','2','constrained'};
    d=3;
    alpha=rand(d,1);
    alphaNormalized=2*pi^2*mean(abs(alpha));    
    alpha=alpha./alphaNormalized;
    delta=-1;
    case 'Rastrigin'
    %Rastrigin
    d=2;
    alpha=[10 10]';
    delta=1;
    problem='cosine';
    %solvers={'1','2','3','32','alpha'}; all solvers used for this experiment
    solvers={'2','32','3'};
    otherwise
        error('did not recognize experiment name');
end



    %% choose problem and compute bounds on derivatives
    switch problem
        case 'cosine'
            [c0,f,g,H,L1,L2,L3,gt]=get_cosine_exp_bounds(d,alpha,delta);
        case 'branin'
            [c0,f,g,H,L1,L2,L3,gt]=get_branin_bounds();
        case 'camel'
            [c0,f,g,H,L1,L2,L3,gt]=get_camel_bounds();
        case 'gold'
            [c0,f,g,H,L1,L2,L3,gt]=get_gold_bounds();
        case 'shubert'
            [c0,f,g,H,L1,L2,L3,gt]=get_shubert_bounds();
        case 'hartman3'
            [c0,f,g,H,L1,L2,L3,gt]=get_hartman3_bounds();
        case 'hartman6'
            [c0,f,g,H,L1,L2,L3,gt]=get_hartman6_bounds();
        case 'shekel5'
            [c0,f,g,H,L1,L2,L3,gt]=get_shekel5_bounds();
        case 'shekel7'
            [c0,f,g,H,L1,L2,L3,gt]=get_shekel7_bounds();
        case 'shekel10'
            [c0,f,g,H,L1,L2,L3,gt]=get_shekel10_bounds();
        otherwise
            error('did not recongnize problem name');
    end



%% solve
solver_num=length(solvers);

ub_vec=inf(solver_num,1);
time_vec=inf(solver_num,1);
computed_error_vec=inf(solver_num,1);
generation_cell=cell(solver_num,1);
survivingCubes=cell(solver_num,1);
x_opt_cell=cell(solver_num,1);
fun_eval_num_cell=cell(solver_num,1);
lin_solve_num_cell=cell(solver_num,1);
fun_eval_num_vec=inf(solver_num,1);
lin_solve_num_vec=inf(solver_num,1);
more_data_cell=cell(solver_num,1);

for ii=1:solver_num
    [ub_vec(ii),x_opt_cell{ii},computed_error_temp,generation_cell{ii},time_vec(ii),survivingCubes{ii},fun_eval_num_cell{ii},lin_solve_num_cell{ii},more_data_cell{ii}]=qBnB(c0,f,g,H,L1,L2,L3,tol,solvers{ii},time_max);
    computed_error_vec(ii)=computed_error_temp(end);
    fun_eval_num_vec(ii)=sum(fun_eval_num_cell{ii});
    lin_solve_num_vec(ii)=sum(lin_solve_num_cell{ii});
end

%% plot results
%iteration comparison
figure;
hold on
for ii=1:solver_num
    plot(log2(cumsum(generation_cell{ii})));
    title('cumulative iterations per search level (log scale) ');
end
legend(solvers);

%time comparison
figure;
hold on
for ii=1:solver_num
    plot(log2(more_data_cell{ii}.elapsed_time));
    title('cumulative time per search level ');
end
legend(solvers);

end