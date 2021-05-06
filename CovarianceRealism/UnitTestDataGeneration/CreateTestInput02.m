
% standard normal distribution of residuals; covariances oversized by factor of 4
function CreateTestInput02
    Covariances=NaN(3,3,1000);
    Residuals=normrnd(0,1,1000,3);
    for i=1:1:1000
        TempResids=normrnd(0,4,1000,3);
        Covariances(:,:,i)=cov(TempResids);
    end
    save('TestInput02','Residuals','Covariances');
end