% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Modification to parts of the dcegm function to include the closed for
% solution for c1 in the case of no taste shocks
% Marta Oliva Riera
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [fun] = 3Afun_closed(grid,fun,param)

    disp('get closed form solutions for comparison');

    fun.c0_closed=zeros(param.na,1);
    fun.c1_closed=zeros(param.na,1);
        
    if (param.sigma<=sqrt(eps))  % option w/o taste shocks
        
        Gam_tilde=param.Gam/(param.Gam-1.0);
        temp=Gam_tilde^param.phi;
        param.temp=temp;
    
        func=@fun_thrs;
        w0=1.0;
        thres_w0=fzero(func,w0,[],param);
        thres_x0=thres_w0-1.0;

        for ac=1:param.na
            
            if (param.phi<eps)
                fun.c0_closed(ac)=2/3+grid.coh_anal(ac)/3;
            else
            
                if (grid.coh_anal(ac)<thres_x0)     % HH decides to work
                    fun.c0_closed(ac)=2/3+grid.coh_anal(ac)/3;
                else                                % HH decides to not work
                    fun.c0_closed(ac)=1/3+grid.coh_anal(ac)/3;
                end
            end
            fun.c1_closed = 2*fun.c0_closed - 1;    % This is what I added for question A.3
            fun.a1_closed=grid.coh_anal-fun.c0_closed;
            fun.a2_closed=grid.coh_anal+1.0-2.0*fun.c0_closed;
        end
    
        
    elseif (param.sigma == 1.0) % option w/ taste shocks, assuming sigma=1
        
        fun.a1_closed = 1/3 * grid.coh_anal - 4/9;
        fun.a2_closed = 1/3 * grid.coh_anal + 1/9;
        fun.c0_closed = 1/3 * grid.coh_anal + 4/9;
        fun.c1_closed = 1/3 * grid.coh_anal + 4/9;
        fun.prob_closed = (10 + 3.0*grid.coh_anal)./(12 + 9.0*grid.coh_anal);
        
    end    

end   % end function fun_closed
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Also changed the graph in this part:

% Closed form solution of dynamic program (fully deterministic + var of taste shock = 1)
if ( sig2e<=sqrt(eps) && ( param.sigma<=sqrt(eps) || (param.sigma==1.0) ) ),
    
    fun = fun_closed(grid,fun,param);
    
    % one cannot do this comparison bec the objects do not live on the same grid,
    % first one would have to interpolate
    c1_intp=zeros(param.na,1);
    for ac=1:param.na
        c1_intp(ac) = interp1(grid.coh(:,2,opt.meth),fun.cons(:,2,opt.meth),grid.coh_anal(ac,2),'linear');
    end
    dist_c1=abs(fun.c1_closed./c1_intp-1.0);
    max_dist_c1=max(dist_c1);
    disp(['maximum distance from closed form solution: ', num2str(max_dist_c1)]);
    
    if (opt.plt)
        xla='$x$';
        yla='$c$';
        leg={'$c_{1}$ - analytical','$c_{1}$ - numerical'};
        tit=['comparision of c(1): analytical and numerical solution with $\sigma$ = ',num2str(param.sigma)] ;
        outpath=[];
        grname='c1_anal_num';
        ax=[];
        onesideplt_marker([grid.coh_anal,grid.coh(:,2,opt.meth)],[fun.c1_closed,fun.cons(:,2,opt.meth)],xla,yla,leg,tit,outpath,grname,ax,opt.pr)
    end
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
