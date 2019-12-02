% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Modification to fun_aggr to improve speed
% Marta Oliva Riera
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for jc=2:nj
    TT = zeros(ny,nx,ny,nx);    % transfer function
    
    for xc=1:nx
            for yc=1:ny
                for ycc=1:ny

                % income (wages and pensions) in current period/age:
                inc=epsi(jc)*netw*gridy(ycc)+(1-epsi(jc))*pens;
                
                % cash on hand: x=a*(1+r)+y = s(-1)*(1+r)+y;
                cah=inc+(1.0+r)*gridsav(xc);
                
                [vals,inds]=basefun(gridx(jc,ycc,:),cah,nx);
                
                TT(ycc,inds(1),yc,xc)=vals(1)*pi(yc,ycc);
                TT(ycc,inds(2),yc,xc)=vals(2)*pi(yc,ycc);
                
                for xcc=1:nx
     
                    % transfer distribution:
                    Phi(jc,ycc,xcc)=Phi(jc,ycc,xcc)+Phi(jc-1,yc,xc)*TT(ycc,xcc,yc,xc)*sr(jc-1);
                
                    if Phi(jc-1,yc,xc) == 0;
                        break
                    end
                    
                 end;    
            end;    
         end; 
    end;
end;

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [Phi,PhiAss,ass]=func_aggr(gridx,gridsav,cfun,gridass)

global r nj nx ny pi gridy netw pens sr epsi pini frac totpop

disp('aggregation and cross-sectional measure');

% Compute Cross sectional distributions and aggregate variables
Phi = zeros(nj,ny,nx);          % distribution of assets conditional by age and shock
PhiAss = zeros(nx,1);             % distribution of assets

% Distribution of newborns over cash at hand
for yc=1:ny
    
    % income (wages and pensions) in current period/age:
    inc=epsi(1)*netw*gridy(yc)+(1-epsi(1))*pens;
    
    % initial cash-on-hand:
    cahini=inc;
    
    [vals,inds]=basefun(gridx(1,yc,:),cahini,nx);
    Phi(1,yc,inds(1))=vals(1)*pini(yc)*frac(1);
    Phi(1,yc,inds(2))=vals(2)*pini(yc)*frac(1);
end;

for jc=2:nj     % this is the modified part: 
    TT = zeros(ny,nx,ny,nx);    % transfer function
    
    for xc=1:nx
            for yc=1:ny
                for ycc=1:ny

                % income (wages and pensions) in current period/age:
                inc=epsi(jc)*netw*gridy(ycc)+(1-epsi(jc))*pens;
                
                % cash on hand: x=a*(1+r)+y = s(-1)*(1+r)+y;
                cah=inc+(1.0+r)*gridsav(xc);
                
                [vals,inds]=basefun(gridx(jc,ycc,:),cah,nx);
                
                TT(ycc,inds(1),yc,xc)=vals(1)*pi(yc,ycc);
                TT(ycc,inds(2),yc,xc)=vals(2)*pi(yc,ycc);
                
                for xcc=1:nx
     
                    % transfer distribution:
                    Phi(jc,ycc,xcc)=Phi(jc,ycc,xcc)+Phi(jc-1,yc,xc)*TT(ycc,xcc,yc,xc)*sr(jc-1);
                
                    if Phi(jc-1,yc,xc) == 0;
                        break
                    end
                    
                 end;    
            end;    
         end; 
    end;
end;


% Check that for each country distribution sums to 1
sumprob=sum(sum(sum(Phi(:,:,:))));
if ( ( sumprob < 0.999 ) || ( sumprob > 1.001) ),
    beep; beep; beep;
    warning('distribution fucked up');
end;

% Check if Grid is Big enough
sumprob=sum(sum(Phi(:,:,nx)));
if (sumprob > 0.001 ),
    beep; beep; beep;
    warning('grid too small -- increase your grid');
    pause
end;

ass=0.0;
cons=0.0;
lab=0.0;
ret=0.0;

% aggregation
for jc=1:nj
    for yc=1:ny
        for xc=1:nx,
            PhiAss(xc)=PhiAss(xc)+Phi(jc,yc,xc);
            
            % asset holdings = capital stock in general equilibrium
            ass=ass+totpop*Phi(jc,yc,xc)*gridsav(xc);
            
            cons=cons+totpop*Phi(jc,yc,xc)*cfun(jc,yc,xc);
            
            lab=lab+totpop*Phi(jc,yc,xc)*gridy(yc)*epsi(jc);
            ret=ret+totpop*Phi(jc,yc,xc)*gridy(yc)*(1.0-epsi(jc));
        end;
    end;
end;


% ---------------------------------------------------------------------
    function [vals,inds]=basefun(grid_x,x,nx)
        % this subroutine returns the values and the indices of the two basis
        % functions that are positive on a given x in the grid_x
        
        % MF function to lookup the current position
        i=lookup(grid_x,x,0);
        
        if ( (i+1)>nx),
            inds(1)=nx;
            inds(2)=nx;
            vals(2)=0.0;
            vals(1)=1.0;
        elseif (i==0),
            inds(1)=1;
            inds(2)=1;
            vals(1)=1.0;
            vals(2)=0.0;
        else
            inds(1)=i;
            inds(2)=i+1;
            dist = grid_x(i+1)-grid_x(i);
            vals(2)=( x-grid_x(i) )/dist;
            vals(1)=( grid_x(i+1)-x )/dist;
        end;
        
    end 	% end function basefun
% ---------------------------------------------------------------------




end     % end function func_aggr
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++