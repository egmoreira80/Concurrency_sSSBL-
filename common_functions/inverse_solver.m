function [T,J,stat,indms] = inverse_solver(L,parcellation,W,Svv,Nsamples,ismethod,isfield)
%% MEG Inverse solution filter
if ismethod == 1 % sSSBL
    [miu,sigma_post,T]    = sSSBLpp(Svv,L,Nsamples,W,parcellation);
    if (isfield == 2) || (isfield == 3) %3D Lead Field
        miu               = sum(reshape(abs(miu),3,length(L)/3),1)';
        sigma_post        = sum(reshape(abs(sigma_post),3,length(L)/3),1)';
    end
    stat                  = sqrt(miu./sigma_post);
    indms                 = find(stat > 1);
    J                     = miu;
elseif ismethod > 1
    indms                 = [1:length(L)]';
    p                     = size(L,1);
    Ip                    = eye(p);
    scaleL                = sqrt(sum(abs(diag(L*L')))/p);
    L                     = L/scaleL;
    scaleV                = (sum(abs(diag(Svv)))/p);
    Svv                   = Svv/scaleV;
    switch ismethod
        case 2
            gamma1                = 0.01;
            gamma2                = 0.1;
            delta_gamma           = 0.01;
            gamma_grid            = gamma1:delta_gamma:gamma2;
            gcv                   = zeros(length(gamma_grid),1);
            count                 = 1;
            for gamma = gamma_grid
                [T,Wout]          = mkfilt_eloreta(L,gamma);
                T                 = T';
                Txiv              = Ip - L*T;
                gcv(count)        = (1/p)*sum(abs(diag(Txiv*Svv*Txiv')))/((1/p)*sum(abs(diag(Txiv))))^2;
                disp(['eloreta gcv param_' , num2str(gamma)])
                count             = count + 1;
            end
            [gcv_opt,idx_gamma]   = min(gcv);
            gamma                 = gamma_grid(idx_gamma);
            [T,Wout]              = mkfilt_eloreta(L,gamma);
        case 3
            gamma                 = sum(abs(diag(Svv)))/(length(Svv)*100);
            [T,T1,Wout]           = mkfilt_lcmv(L,Svv,gamma);
    end
    T                     = transpose(T);
    J                     = abs(diag(T*Svv*T'));
    if isfield == 3 %3D Lead Field
        J                 = sum(reshape(abs(J),3,length(L)/3),1)';
    end
    stat                  = J;
end
end