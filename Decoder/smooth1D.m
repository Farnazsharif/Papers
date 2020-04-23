%function Sdata=smooth1D(data,halfwidth,dim)
%
%smooth data by convolving it with a gaussian of
%halwidth 'halwidth'. halfwidth in number of points.
%if dim=1, smooth column, else, smooth rows


function Sdata=smooth1D(data,halfwidth,dim)

Fullwin=(-halfwidth:halfwidth)';
Svector=exp(-Fullwin.^2/(halfwidth/2)^2);

[a,b]=size(data);
if a>1&b>1

	if dim==1
		for ii=1:length(data(1,:))
		smoothed=conv(data(:,ii),Svector)/sum(Svector);
		Sdata(:,ii)=smoothed((halfwidth+1):(end-halfwidth));
		end
	else
		for ii=1:length(data(:,1))
		smoothed=conv(data(ii,:),Svector)/sum(Svector);
		Sdata(ii,:)=smoothed((halfwidth+1):(end-halfwidth));
		end
	end
else

	smoothed=conv(data,Svector)/sum(Svector);
	Sdata=smoothed((halfwidth+1):(end-halfwidth));

end








