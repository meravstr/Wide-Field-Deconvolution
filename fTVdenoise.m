function [ xsol ] = fTVdenoise( lambda,y )
%TVDENOISE implement A direct Algorithm for 1D Total Variation, Laurent
%Condat IEEE 2013

N = length(y);
xsol = zeros(N,1);

k = 1;
k0 = 1;
kmin = 1;
kplus = 1;
numin = y(1)-lambda;
numax = y(1)+lambda;
umin = lambda;
umax = -lambda;
inloop = true;
while inloop
    if k==N
        xsol(N) = numin+umin;
        inloop = false;
    else
        try
            while k<N 
                   if y(k+1)+umin<numin-lambda
                       xsol(k0:kmin) = numin;
                       k = kmin+1;
                       k0 = kmin+1;
                       kmin = kmin+1;
                       kplus = kmin+1;
                       numin = y(k);
                       numax = y(k)+2*lambda;
                       umin = lambda;
                       umax = -lambda;
                   else
                        if y(k+1)+umax>numax+lambda
                           xsol(k0:kplus) = numax;
                           k = kplus+1;
                           k0 = kplus+1;
                           kmin = kplus+1;
                           kplus = kplus+1;
                           numin = y(k)-2*lambda;
                           numax = y(k);
                           umin = lambda;
                           umax = -lambda;
                        else
                            k = k+1;
                            umin = umin+y(k)-numin;
                            umax = umax+y(k)-numax;
                            if umin>=lambda
                                numin = numin+(umin-lambda)/(k-k0+1);
                                umin = lambda;
                                kmin = k;
                            end
                            if umax<=-lambda
                                numax = numax+(umax+lambda)/(k-k0+1);
                                umax = -lambda;
                                kplus = k;
                            end
                        end
                   end
            end
            if umin<0
                xsol(k0:kmin) = numin;
                k = kmin+1;
                k0 = kmin+1;
                kmin = kmin+1;
                numin = y(k);
                umin = lambda;
                umax = y(k)+lambda-numax;
            else
                if umax>0
                    xsol(k0:kplus) = numax;
                    k = kplus+1;
                    k0 = kplus+1;
                    kplus = kplus+1;
                    numax = y(k);
                    umax = -lambda;
                    umin = y(k)-lambda-numin;
                else
                    xsol(k0:N) = numin+umin/(k-k0+1);
                    inloop = false;
                end
            end
        catch
           %warning('am i stuck?');
           inloop = false;   
        end 
    end
end

end

