function [ub,x_opt,computedError,generationSizeVec,totalTime,survivingCubes,fun_eval_num,lin_solve_num,more_data]=qBnB(c0,f,g,H,L1,L2,L3,tol,bound_type,time_max)
%global optimization using quasibnb of the function f in the given cube,
%using the bound on the hessian L2, up to an error tolerance tol


%doDisplay=getoptions(params,'doDisplay',true);
%iterMax=getoptions(params,'iterMax',10^6); %maximal number of evaluations we allow
doDisplay=false;
iterMax=inf;
list={c0};
d=length(c0.x);
r0=norm(c0.h);

%start qBnB
runTime=tic;

depth=0;
ub=inf;
local_ub=inf;
generationSize=length(list);
doStop=false;
totalRuns=0;
generationSizeVec=[];
keptIndNum=[];
allIterInd=0;
%ubVec=inf(iterMax,1);
ubVec=[];
fun_eval_num=[];
lin_solve_num=[];
elapsed_time=[]; %we only need this for one experiment
while ~doStop
    depth=depth+1;
    generationSizeVec(depth)=generationSize;
    objective=nan(generationSize,1);
    qlb=nan(generationSize,1);
    
    %solve on all generations
    %x_cell=cell(generationSize,1);
    for ii=1:generationSize
        %check timing condition
        if mod(ii,10)==0
            time_now=toc(runTime);
            if time_now>time_max
                fprintf('stopped qBnB because time limit was exceeded');
                doStop=true;
            end
        end
        if doStop
            break;
        end
        
        allIterInd=allIterInd+1;
        cube=list{ii};
        %compute local_ub,x,, qlb
        switch bound_type
            case '1'
                [x,local_ub,qlb(ii),fun_eval_num(allIterInd),lin_solve_num(allIterInd)]=compute_lip_bounds(cube,f,L1);
                
            case '2'
                [x,local_ub,qlb(ii),fun_eval_num(allIterInd),lin_solve_num(allIterInd)]=compute_2bounds(cube,f,L2);
% seems to be superflous            case '2constrainedB'
%                 [x,local_ub,qlb(ii),fun_eval_num(allIterInd),lin_solve_num(allIterInd)]=compute_2bounds_constrained_b(cube,c_0,f,L2);
            
            case '2lb'
                [x,local_ub,qlb(ii),fun_eval_num(allIterInd),lin_solve_num(allIterInd)]=compute_2lb_bounds(cube,f,g,L2);
            case '3'
                %[x,local_ub,qlb(ii)]=compute_3bounds_old_version(cube,f,g,H,L3);
                [x,local_ub,qlb(ii),fun_eval_num(allIterInd),lin_solve_num(allIterInd)]=compute_3bounds(cube,f,g,H,L3,tol/100);
            case '32'
                r=norm(cube.h);
                if (0.5*L2*r^2<3*L3*r^3) || (2*r+norm(cube.x-c0.x)>r0)
                    [x,local_ub,qlb(ii),fun_eval_num(allIterInd),lin_solve_num(allIterInd)]=compute_2bounds(cube,f,L2);
                else
                    %[x,local_ub,qlb(ii)]=compute_3bounds(cube,f,g,H,L3);
                    [x,local_ub,qlb(ii),fun_eval_num(allIterInd),lin_solve_num(allIterInd)]=compute_3bounds(cube,f,g,H,L3,tol/100);
                end
            case 'alpha'
                %min_eig_est=L2;
                min_eig_est=min(eig(H(cube.x)))-L3*norm(cube.h);
                alpha_vec=max(0,-0.5*min_eig_est)*ones(d,1);
                [x,local_ub,qlb(ii),fun_eval_num(allIterInd),lin_solve_num(allIterInd)]=compute_alpha_bounds(cube,f,g,H,alpha_vec);
            case 'constrained'
                [x,local_ub,qlb(ii),fun_eval_num(allIterInd),lin_solve_num(allIterInd)]=compute_2bounds_constrained(cube,c0,f,L2);
            otherwise
                error('did not recongnize solver name');
                
        end
        if local_ub<ub
            ub=local_ub;
            x_opt=x;
        end
        ubVec(allIterInd)=ub;
    end
    if doStop %stopped due to time overflow
        totalRuns=totalRuns+ii;
        UBall=ub;
    else
        totalRuns=totalRuns+generationSize;
        %prepare for next generation
        keepInd=qlb<=ub;
        keptIndNum(depth)=sum(keepInd);
        prevList=list(keepInd,1);
        
        UBall=ub;
        LBall=min(qlb);
        
        
        constrained_split=false; %this option currently disabled
        if constrained_split
            preBoundaryList=list(~keepInd,1);
            boundaryList=getBoundaryList(preBoundaryList,c0);
        else
            boundaryList={};
        end
        list=getNewListBinarySplit(prevList,d); %FILL IN
        list={list{:},boundaryList{:}};
        list=list';
        generationSize=length(list);
        if isempty(list)
            doStop=true;
            fprintf('stopped qBnB because list is empty');
        end
        computedError(depth)=UBall-LBall;
        if UBall-LBall<tol
            doStop=true;
            fprintf('stopped qBnB because requested energy tolerance was reached \n');
        elseif totalRuns+generationSize>iterMax
            doStop=true;
            warning('stopped BnB before reaching requested tolerance \n');
        end
        
    end
    
    
    
    elapsed_time(depth)=toc(runTime);
end

more_data.elapsed_time=elapsed_time;
more_data.computed_error=computedError;
survivingCubes=prevList;

totalTime=toc(runTime);

if doDisplay
    figure;
    plot(generationSizeVec);
    title('number of BnB iterations per generation');
    figure;
    plot(log2(generationSizeVec));
    title('log2 number of BnB iterations per generation');
    figure;
    plot(ubVec);
    title('upper bound vs iteration num');
end

end
