function m = read_ascii(fname)
%
% pts=read_mat('corners.mat');
%

m=-1;
fid = fopen(fname, 'r');
if( fid < 0 )
   return
end

tline = fgetl(fid);
[nrc] = sscanf(tline,'%d %d ');

nr = nrc(1);
nc = nrc(2);
m = zeros(nr,nc);

for r=1:nr % read in collumns of run length data.
    tline = fgetl(fid);

	 df = sscanf(tline,'%g ');
     if( size(df,1) >= nc)
         m(r,:)  = df(1:nc);
     else
         disp('duff line length ');
     end

end 

fclose(fid);

return;




