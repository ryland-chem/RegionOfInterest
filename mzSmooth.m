function [smoothedMZ] = mzSmooth(xic)

sz = size(xic);

numbChan = sz(2);

smoothedMZ = [];

for i = 1:numbChan
   
    smoothedMZ(:,i) = smooth(xic(:,i), 10);
    
end

end