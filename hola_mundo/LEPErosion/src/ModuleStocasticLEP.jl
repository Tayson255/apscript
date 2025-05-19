module ModuleStocasticLEP

using Distributions
using PDMats


export logRegression, simplify_randomness


function logRegression(V,N)

    X=log10.(V)
    Y=log10.(N)
    n=size(X)

    Z=X.-mean(X)
    ZZ=[ones(n[1],1) Z]

    #Procedure 2
    SS_zy=Z'*Y
    SS_z=Z'*Z
    beta_0=mean(Y)
    beta_1=SS_zy/SS_z
    beta=[beta_0;beta_1]
    SS_R=Y'*Y-beta'*(ZZ'*ZZ)*beta

    #Procedure 3
    MS_R=SS_R/(n[1]-2);


    mu_beta = [beta_0, beta_1]
    sigma_beta=[[MS_R/n[1], 0 ] [0, MS_R/SS_z]]

    tM=[[1.0, 0] [-mean(X), 1]];
    mu_alpha = tM*mu_beta;
    sigma_alpha= tM*sigma_beta * tM';

    D=MvNormal(mu_alpha,PDMat(sigma_alpha))

    return D

end

function simplify_randomness(D::MvNormal)
    mu = mean(D)
    sigma = cov(D)
    sigma[2,1] = sign(sigma[2,1]) * sqrt(sigma[1,1] * sigma[2,2])
    sigma[1,2] = sigma[2,1]
    return MvNormal(mu, PDMat(sigma))
end



end