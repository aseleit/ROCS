function [ dx ] = nae_su(F, B, x, x0, method, parameter, iter)
% update the solution to NAE

if(strcmp(method,'NEWTON'))
    if(norm(F)<1e5*eps)
        dx = zeros(size(B,2),1);
    else
        dx = -B\F';
    end
elseif(strcmp(method,'SFHM1'))
    v = parameter(1);
    m = parameter(2);
    if(norm(F)<1e5*eps)
        dx = zeros(size(B,2),1);
    else
        t = iter;
        QdQ = v/(1+t)^m;
        Q = exp(v/(1-m)*((1+t)^(1-m)-1));
        u = B'*F+(x-x0)/Q;
        dx = -0.5*QdQ*norm(F)^2/norm(u)^2*u;
    end
elseif(strcmp(method,'SFHM2'))
    if(norm(F)<1e5*eps)
        dx = zeros(size(B,2),1);
    else
        t = iter;
        QdQ = 1;
        Q = exp(t);
        u = B'*F+(x-x0)/Q;
        dx = -0.5*QdQ*norm(F)^2/norm(u)^2*u;
    end
elseif(strcmp(method,'SNHM1'))
    v = parameter(1);
    m = parameter(2);
    if(norm(F)<1e5*eps)
        dx = zeros(size(B,2),1);
    else
        t = iter;
        QdQ = v/(1+t)^m;
        u = B'*F;
        dx = -0.5*QdQ*norm(F)^2/norm(u)^2*u;
    end
elseif(strcmp(method,'SNHM2'))
    if(norm(F)<1e5*eps)
        dx = zeros(size(B,2),1);
    else
        QdQ = 1;
        u = B'*F;
        dx = -0.5*QdQ*norm(F)^2/norm(u)^2*u;
    end
elseif(strcmp(method,'GOIA'))
    gamma = parameter(1);
    if(norm(F)<1e5*eps)
        dx = zeros(size(B,2),1);
    else
        A = B*B';
        v1 = A*F;
        v2 = B*F;
        ac = (norm(v1)^2*norm(v2)^2 - (v1'*v2)^2)/norm((v1'*F)*v2 - (v2'*F)*v1)^2;
        alphac = (ac*(F'*v1)*(F'*v2) - (v1'*v2))/(norm(v2)^2-ac*(F'*v2)^2);
        v = v1 + alphac*v2;
        if (abs(F'*v) <= 1e-6)
            alphac = 0;
            v = v1 + alphac*v2;
        end
        u = alphac*F + B'*F;
        dx = -(1-gamma)*(F'*v)/norm(v)^2*u;
    end
elseif (strcmpi(method,'LSQ'))
    if(norm(F)<1e5*eps)
        dx = zeros(size(B,2),1);
    else
        dx  = -inv(B'*B)*B'*F;
    end
    
elseif (strcmpi(method,'FTIM'))
    nu = parameter(1);
    if(norm(F)<1e5*eps)
        dx = zeros(size(B,2),1);
    else
        t = iter;
        qt = exp(t);
        dx  = -nu/qt*F;
    end
    
elseif (strcmpi(method,'OIA-ODV'))
    gamma = parameter(1);
    if(norm(F)<1e5*eps)
        dx = zeros(size(B,2),1);
    else
        A = B*B';
        v1 = A*F;
        v2 = B*F;
        alphao = dot(TripleVec(v1,F,v2),v1)/dot(TripleVec(v2,F,v1),v2);
        v = v1 + alphao*v2;
        u = alphao*F + B'*F;
        dx = -(1-gamma)*dot(F,v)/norm(v)^2*u;
    end
else
    error('Unsupport NAE Solver')
end

end
