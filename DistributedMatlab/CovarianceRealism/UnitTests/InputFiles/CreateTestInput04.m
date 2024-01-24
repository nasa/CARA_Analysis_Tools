
% student's t-distribution used instead of normal distribution
% symmetric distribution but not Gaussian, so might pass some of the tests with a scale factor applied
function CreateTestInput04
    Covariances=NaN(3,3,1000);
    Residuals=trnd(2,1000,3);
    for i=1:1:1000
        TempResids=trnd(2,1000,3);
        Covariances(:,:,i)=cov(TempResids);
    end
    save('TestInput04','Residuals','Covariances');
end