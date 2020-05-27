function [minima,problem_names]=get_dixon_optimal()
problem_names={'branin','camel','gold','shubert','hartman3','shekel5','shekel7','shekel10','hartman6'};
minima=inf(9,1);
for jj=1:9
    %% choose problem and compute bounds on derivatives
    switch problem_names{jj}
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
    end
    if ~isnan(gt)
        minima(jj)=f(gt(:,1));
    end


end
minima(4)=-186.7309;
results_folder='/home/postdoc/nadavd/.sshfs/corwin/ytmp/cloud/nadavd/dropbox/Dropbox/smoothQuasiBnB/code/results';
cd(results_folder);
save('dixon_minima','minima','problem_names');
end