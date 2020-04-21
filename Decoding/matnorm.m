%function m=matnorm(M,normaxe)
%
%normaxe==1 normalize the columns
%normaxe==2 normalize the rows

function m=matnorm(M,normaxe)

m=M;

if normaxe==1

	for ii=1:length(m(1,:))
		m(:,ii)=m(:,ii)/(max(m(:,ii))+eps);
	end
	
elseif normaxe==2

	for ii=1:length(m(:,1))
		m(ii,:)=m(ii,:)/(max(m(ii,:))+eps);
	end
	
end
