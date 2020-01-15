function [bias,nrmse,rmse,si,symr,cc]=wave_stats(x,y)

N=length(find(~isnan(x.*y)));
ii=find(~isnan(x.*y));
bias=1/N*sum(y(ii)-x(ii));
nrmse=sqrt(sum((y(ii)-x(ii)).^2)/sum(x(ii).^2));
rmse=sqrt(1/(N-1)*sum((y(ii)-x(ii)).^2));
si=rmse/mean(x(ii));
symr=sqrt(sum(y(ii).^2)./sum(x(ii).^2));
%cc=sum((x(ii)-mean(x(ii))).*(y(ii)-mean(y(ii))))./(sqrt(sum((x(ii)-mean(x(ii))).^2).*sum((y(ii)-mean(y(ii))).^2)));
c=corrcoef(x(ii),y(ii)); %same as above
cc=c(2);

%disp(['Bias: ',num2str(bias)]);
%disp(['RMSE: ',num2str(rmse)]);
%disp(['NRMSE: ',num2str(nrmse)]);
%disp(['SI: ',num2str(si)]);
%disp(['symr: ',num2str(symr)]);
%disp(['cc: ',num2str(cc)]);

disp([num2str(bias),' ', num2str(rmse),' ', num2str(nrmse),' ', num2str(si), ' ',num2str(symr), ' ',num2str(cc)]);
