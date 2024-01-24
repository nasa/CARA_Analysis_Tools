
% standard normal distribution of residuals; covariances undersized by factor of 3
function CreateTestInput03
    Covariances=NaN(3,3,1000);
    Residuals=normrnd(0,1,1000,3);
    for i=1:1:1000
        TempResids=normrnd(0,.333,1000,3);
        Covariances(:,:,i)=cov(TempResids);
    end
    save('TestInput03','Residuals','Covariances');
end