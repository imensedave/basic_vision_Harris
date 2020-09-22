function [lim] = read_rlm(fname)
%
% fiction to read a run length map file .
%
% run length encoding is good for image region label images.
%
%
%

fid=fopen(fname);  
if( fid < 0 )
	['bad',fname];
	lim = [];
	return;
end
tline = fgetl(fid);
[nrc] = sscanf(tline,'%d %d');
lim = zeros(nrc(1), nrc(2) );
nr = nrc(1);, nc = nrc(2);

for c=1:nrc(2) % read in collumns of run length data.
    tline = fgetl(fid);
	 df = sscanf(tline,'%d ');
	 numd = size( df,1);
	 lbls = df(1:2:end);
	 lgs = df(2:2:end);
	 legs = [0;cumsum(lgs)];
	 
	 for k=1:numd/2  % index of labels
		 lim(legs(k)+1:legs(k+1),c) = lbls(k);
	 end
		 
end

 fclose(fid);


return;
