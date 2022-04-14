%%Region of interest selection for 1D GC-MS data using pseudo-fisher ratios
%
%(c) 2022 Ryland T. Giebelhaus, Michael D.S. Armstrong, A. Paulina de la
%Mata, James J. Harynuk
%
%Utilisation of a moving window function to calculate the f ratios for a
%particular region of interest. Smaller windows are more sensitive to
%smaller features, but may split larger regions of interest. Larger windows
%are less sensitive to small regions and create larger regions of interest.
%Currently this returns a vector of probabilities for one-dimensional data,
%and it is up to the user to determine the significance level.
%
%
%
%v1.2

%main branch, code does the following
%takes xic data, scan window (between 5-40 ideally) and p value cutoff from
%user
%returns the pvals for each scan, the modPVans which is the pvals less the
%pvalues below the cutoff, the ticData which is just the XIC data converted
%to TIC, noiseDropped which is the TIC without regions of non-interest, and
%boolCutOff which is a binary yes or no if a p-value is above (1) or below
%(0) the cutoff input by the user. This program also outputs graphs.

%%%inputs
%%data: M x N array of ion intensities, where M is scans and N is ion m/z
%%wndw: moving window size. Ideally this should be about the approximate
%width of the average peak in your separation. 10 is a good starting point
%if unsure
%%CutOff: probability cutoff. 0.7 is a good place to start.

%%%outputs
%%pv: probability values for each scan. Is overlaid on output plot by
%default
%%modPVans: modified probability values, above the CutOff threshold.
%%ticData: Total Ion Chromatogram of the data. Plotted by default.
%%noiseDroppedTIC: TIC with non ROIs dropped.
%%noiseDropped: all ions and scans with noise (non ROIs) dropped. This is
%intended for use in further analysis.
%%boolCutOff: 0 = scan not in ROI, 1 = scan in ROI.

function [pv, modPVans, ticData, noiseDroppedTIC, noiseDropped, boolCutOff] = froiispeedop(data, wndw, CutOff)

%bool to print graph
%at the start so user can input then let run
prompt = 'Output graph (y/n)';
choicePrint = input(prompt, 's');

%Initialisation
sz = size(data);

%Start at 1:wndw
indx(1) = 1;
indx(2) = wndw;

mat = zeros(sz(1),1);

iter = 1;

chisqVec = [];

sumChiSq = zeros(sz(1),1);

pv = [];

%for loop 
for i = 1:(sz(1)-wndw)
    
    %C is the autoscaled matrix
    C = (data(indx(1):indx(2),:) - mean(data(indx(1):indx(2),:)))./std(data(indx(1):indx(2),:));
    
    %C converts NA values to 0
    C(isnan(C)) = 0;
    
    %s are the singular values
    s = svds(C,2);
    
    %f is the fisher ratio
    f(i) = ((s(1)^2)/(wndw-1))/((s(2)^2)/(wndw-2)); %#ok
    
    %vector of probabilities
    mat(indx(1):indx(2),1) = fcdf(f(i),wndw - 1, wndw - 2);
    
    parfor j = 1:sz(1)
       
        chisqVect(1,j) = chi2inv(mat(j),1);
        
    end
    
    chisqVect = chisqVect';
    
    sumChiSq = sumChiSq + chisqVect;
    
    chisqVect = [];
    
    mat = zeros(sz(1),1);
    
    %increase the iteration number, change the region where the window is
    %active

    indx(1) = indx(1) + 1;
    indx(2) = indx(2) + 1;
    
end


for ii = 1:sz(1)
    
   %ii is the observation number we are looping through
   %the degrees of freedom is the number of observations, since there is no
   %class information
   %getting DOF through relationship between the number of observations and
   %the window size
    for k = 1:sz(1)

        if k <= wndw

            dof(k) = k;

        elseif k > sz(1) - wndw

            dof(k) = wndw - (k - (sz(1) + 1 - wndw));

        elseif k >= wndw

            dof(k) = wndw;

        end

    end
   
   parfor jj = 1:sz(1)
       
       %every chi sq is translated to a p value with a dof of 1.
       pv(jj) = chi2cdf(sumChiSq(jj), dof(jj));
       
   end
   

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
ticData = zeros(sz(1),1);

for i = 1:sz(1)
    
    ticData(i) = sum(data(i,:));

end

%transpose to make it go the right direction
ticData = ticData';

%column of zeros for noise dropped from TIC
%for speed
noiseDroppedTIC = zeros(sz(1), 1);

%empty matrix for entire chromatogram with just regions
%to preserve the spectral data and export it
noiseDropped = zeros(sz(1), sz(2));

%column of bools, 1 for above p cutoff and 0 for below
boolCutOff = zeros(sz(1), 1);

for i = 1:sz(1)
   
    %if the p value wasnt cutoff carry the tic value over
    if modPVans(i) > 0
        
        %for TIC data
        noiseDroppedTIC(i) = ticData(i);
        
        %for spectral data
        noiseDropped(i,:) = data(i,:);
        
        %for the cutoff
        boolCutOff(i) = 1;
    
    %if the p value was cutoff then drop tic value
    elseif modPVans(i) == 0
        
        noiseDroppedTIC(i) = 0;
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
    
else
    
end


end