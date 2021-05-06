
% standard normal distribution of residuals; covariances to match
function CreateTestInput01
    Covariances=NaN(3,3,1000);
    Residuals=normrnd(0,1,1000,3);
    for i=1:1:1000
        TempResids=normrnd(0,1,1000,3);
        Covariances(:,:,i)=cov(TempResids);
    end
    save('TestInput01','Residuals','Covariances');
end