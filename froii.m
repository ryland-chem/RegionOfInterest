%%Region of interest selection for 1D GC-MS data using pseudo-fisher ratios
%
%(c) 2021 Michael Sorochan Armstrong & Ryland T. Giebelhaus
%
%Utilisation of a moving window function to calculate the f ratios for a
%particular region of interest. Smaller windows are more sensitive to
%smaller features, but may split larger regions of interest. Larger windows
%are less sensitive to small regions and create larger regions of interest.
%Currently this returns a vector of probabilities for one-dimensional data,
%and it is up to the user to determine the significance level.
%
%Currently not optimised for speed.
%
%v1.1

%main branch, code does the following
%takes xic data, scan window (between 5-40 ideally) and p value cutoff from
%user
%returns the pvals for each scan, the modPVans which is the pvals less the
%pvalues below the cutoff, the ticData which is just the XIC data converted
%to TIC, noiseDropped which is the TIC without regions of non-interest, and
%boolCutOff which is a binary yes or no if a p-value is above (1) or below
%(0) the cutoff input by the user. This program also outputs graphs.

%to normalize the x axis to time need scan speed
%will assume 200 Hz for now
function [pv, modPVans, ticData, noiseDropped, boolCutOff] = froii(data, wndw, CutOff)

%bool to print graph
%at the start so user can input then let run
prompt = 'Output graph (y/n)';
choicePrint = input(prompt, 's');

%Initialisation
sz = size(data);

%Start at 1:wndw
indx(1) = 1;
indx(2) = wndw;

%number of scans in dataset
numbScans = sz(1);

secTime = [];

for i = 1:sz(1)
   
    %converts to seconds
    %i - 1 as the first scan is 0 sec
    secTime(i) = (i - 1)/200;
    
end

mat = [];

iter = 1;

%Using a while loop instead of doing arithmetic; change this later
while indx(2) <= sz(1)
    
    %C is the autoscaled matrix
    C = (data(indx(1):indx(2),:) - mean(data(indx(1):indx(2),:)))./std(data(indx(1):indx(2),:));
    
    %C
    C(isnan(C)) = 0;
    
    %s are the singular values
    s = svds(C,2);
    
    %f is the pseudo fisher ratio
    f(iter) = s(1)^2/s(2)^2; %#ok
    
    %matrix of probabilities
    mat(indx(1):indx(2),iter) = fcdf(f(iter),wndw,wndw); %#ok
    
    %increase the iteration number, change the region where the window is
    %active
    iter = iter + 1;
    indx(1) = indx(1) + 1;
    indx(2) = indx(2) + 1;
    
end

for ii = 1:sz(1)
    
   %ii is the observation number we are looping through
   %the degrees of freedom is the number of observations, since there is no
   %class information
   dof(ii) = sum(mat(ii,:) > 0); %#ok
   
   %retrieve the nonzero elements. mat is an off-diagonal matrix.
   nzx2 = nonzeros(mat(ii,:));
   
   for qq = 1:max(size(nzx2))
       %every probability is translated to a chi2 value with a dof of 1.
       x2mat(ii,qq) = chi2inv(nzx2(qq),1); %#ok
   end
   
   %sum the chi2 values
   x2vec(ii) = sum(x2mat(ii,:),2); %#ok
   
   %translate the chi2 values back to
   pv(ii) = chi2cdf(x2vec(ii),dof(ii)); %#ok
   
end

%Change the orientation to a column
pv = pv';

%copy pv's so we dont ruin column
modPVans = pv; 
 
for i = 1:length(pv) 
     
    if pv(i) < CutOff 
         
        modPVans(i) = 0; 
       
    elseif pv(i) > CutOff 
         
        modPVans(i) = pv(i); 
     
    end 
         
end 

%need to drop the noise from the TIC
%first generate the TIC so the user doesnt have to input it
%first calculate the TIC
ticData = zeros(numbScans,1);

for i = 1:numbScans
    
    ticData(i) = sum(data(i,:));

end

%transpose to make it go the right direction
ticData = ticData';

%column of zeros for noise dropped from TIC
%for speed
noiseDropped = zeros(numbScans, 1);

%column of bools, 1 for above p cutoff and 0 for below
boolCutOff = zeros(numbScans, 1);

for i = 1:numbScans
   
    %if the p value wasnt cutoff carry the tic value over
    if modPVans(i) > 0
       
        noiseDropped(i) = ticData(i);
        boolCutOff(i) = 1;
    
    %if the p value was cutoff then drop tic value
    elseif modPVans(i) == 0
        
        noiseDropped(i) = 0;
        boolCutOff(i) = 0;
        
    end
    
end 

%conditional whether to print or not

if choicePrint == 'y'
    
    yyaxis right; hold on; plot(ticData); ylabel('Intensity');
    
        %plots boxes around the ROI
        yyaxis left; hold on; area(boolCutOff);

        %area color
        newcolors = [0.7 0.7 0.7]; %grey
        colororder(newcolors);

        %set transparent
        alpha(0.4);
    
    xlabel('time (s)');
    
else
    
end


end

