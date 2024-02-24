function  c = mahalTune_func_cross_temp(data,theta,angspace,bin_width)

% computes the mahalanobis distances between the data of single test-trial of a particular orientation and the rest of the data,
% for all time-point combinations.

trl_ind=1:size(data,1);
sigma=nan(size(data,1),size(data,3),size(data,2),size(data,2));
% prepare sigmas (slightly faster when done here)
for trl=1:size(data,1);
    trn_dat = data(setdiff(trl_ind,trl),:,:);
    for t1=1:size(data,3);
        sigma(trl,t1,:,:) = covdiag(trn_dat(:,:,t1));
    end
end
c_temp=nan(size(data,1),size(data,3),size(data,3));
d_temp=nan(length(angspace),size(data,3),size(data,3));
reverseStr='';
for trl=1:size(data,1)
    trn_ind = setdiff(trl_ind,trl);
    trn_dat = data(trn_ind,:,:);
    trn_angle = theta(trn_ind);
    msg=sprintf('%d percent\n',round((trl/size(data,1))*100));
    fprintf([reverseStr,msg]);
    reverseStr=repmat(sprintf('\b'),1,length(msg));
    test_dat=squeeze(data(trl,:,:))';
    for b=1:length(angspace)
        m=squeeze(mean(trn_dat(abs(angle(exp(1i*trn_angle)./exp(1i*(theta(trl)-angspace(b)))))<bin_width,:,:),1));
        for t1=1:size(data,3)
            d_temp(b,t1,:) = pdist2(squeeze(m(:,t1))', test_dat,'mahalanobis',squeeze(sigma(trl,t1,:,:)));
        end
    end
    c_temp(trl,:,:)=-mean(bsxfun(@times,cos(angspace),d_temp),1);
end
c=squeeze(mean(c_temp,1));
%%
    function sigma=covdiag(x)
        
        % x (t*n): t iid observations on n random variables
        % sigma (n*n): invertible covariance matrix estimator
        %
        % Shrinks towards diagonal matrix
        % as described in Ledoit and Wolf, 2004
        
        % de-mean returns
        [t,n]=size(x);
        meanx=mean(x);
        x=x-meanx(ones(t,1),:);
        
        % compute sample covariance matrix
        sample=(1/t).*(x'*x);
        
        % compute prior
        prior=diag(diag(sample));
        
        % compute shrinkage parameters
        d=1/n*norm(sample-prior,'fro')^2;
        y=x.^2;
        r2=1/n/t^2*sum(sum(y'*y))-1/n/t*sum(sum(sample.^2));
        
        % compute the estimator
        shrinkage=max(0,min(1,r2/d));
        sigma=shrinkage*prior+(1-shrinkage)*sample;
    end
end

