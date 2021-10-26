%compresses the data due to issue with the MS data
%michael whipped this up

function [mzc,xicc] = mzcompress(mz,xic)

%round the mz vector
mzs = round(mz);

qq = 1;

for ii = min(mzs):max(mzs)
   
    indx = find(mzs==ii);
    
    xicc(:,qq) = sum(xic(:,indx),2); %#ok
    mzc(qq) = ii; %#ok
    qq = qq + 1;
    
end

end
