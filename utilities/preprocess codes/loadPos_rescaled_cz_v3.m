function [post,posx,posy] = loadPos_rescaled_cz_v3(file,scalex,scaley,vfs,track_center)
% this code is for the vnt files created in Circular track experiment
% the NVT file should have several segments, with big gaps
% find out each segment, and then do linear interpolation in each segment

% updated on 03/14/2017
% this code is for circular track
% if there's a missing tracking data, interp it in angle but not linear
% position

% file is the video file
% scale is the ratio that transfer the location from pixels to cm
% vfs is the sample frequency of video

fieldSelection(1) = 1; % Timestamps
fieldSelection(2) = 1; % Extracted X
fieldSelection(3) = 1; % Extracted Y
fieldSelection(4) = 0; % Extracted Angel
fieldSelection(5) = 0;  % Targets
fieldSelection(6) = 0; % Points

% Do we return header 1 = Yes, 0 = No.
extractHeader = 0;
% 5 different extraction modes, see help file for Nlx2MatVT
extractMode = 1; % Extract all data
[t, x, y] = Nlx2MatVT(file,fieldSelection,extractHeader,extractMode);

%tracklength = 193; %track length in cm
% scale = 0.275;% Conversion from pixels to cm, frame is 640 by 480
x = x * scalex; % == transfer to cm ==
y = y * scaley;
t = t ./ 1000000;

if max(x) == 0 & max(y) == 0
    % no video tracking data
    post=t;
    posx=x;
    posy=y;
    return
end

% find the segments
ind = find(diff(t)>1);
ind_seg = [[1,ind+1];[ind,length(t)]]';
post=nan(size(t));
posx=nan(size(x));
posy=nan(size(y));
for nseg = 1:size(ind_seg,1)
    x0 = x(ind_seg(nseg,1):ind_seg(nseg,2));
    y0 = y(ind_seg(nseg,1):ind_seg(nseg,2));
    t0 = t(ind_seg(nseg,1):ind_seg(nseg,2));
    
    %Clear zero values
    ind = find((y0>0 & x0>0)==1 & ~isnan(x0));
    
    if ~isempty(ind)
        x0_center = x0-track_center(1);
        y0_center = y0-track_center(2);
        ang0 = nan(size(x0));
        ind0=find(x0_center>0 & x0~=0);
        ang0(ind0)=mod(atan(-y0_center(ind0)./x0_center(ind0)),2*pi);
        ind0=find(x0_center<0 & x0~=0);
        ang0(ind0)=mod(atan(-y0_center(ind0)./x0_center(ind0))+pi,2*pi);
        
        xx0 = x0_center(ind);
        yy0 = y0_center(ind);
        radius00 = sqrt(xx0.^2 +yy0.^2);
        tt0 = t0(ind);
        ang00 = ang0(ind);
        
        %Fill in the gaps that the camera missed by circular interpolation
        radius0 = interp1(tt0,radius00,t0,'linear');
        
        ang00_unwrap = unwrap(ang00);
        ang0_unwrap_interp = interp1(tt0,ang00_unwrap,t0,'linear');
        ang0_interp = mod(ang0_unwrap_interp,2*pi);
        posx_nseg = radius0.*cos(ang0_interp);
        posy_nseg = radius0.*(-1).*sin(ang0_interp);
        posx_nseg(ind) = xx0;
        posy_nseg(ind) = yy0;
        ind1 = ind(1);
        posx_nseg(1:ind1-1) = posx_nseg(ind1);
        posy_nseg(1:ind1-1) = posy_nseg(ind1);
        ind2 = ind(end);
        posx_nseg(ind2+1:end) = posx_nseg(ind2);
        posy_nseg(ind2+1:end) = posy_nseg(ind2);
        posx_nseg = posx_nseg+track_center(1);
        posy_nseg = posy_nseg+track_center(2);
        
        post_nseg=t0;
        
        % Treshold for how far a rat can move (100cm/s), in one sample (sampFreq =
        % 30 Hz)
        threshold = 120/vfs;
        [posx0,posy0] = removeBadTracking_cz(posx_nseg,posy_nseg,threshold);
        %Clear zero values
        ind = find((posy0>0 & posx0>0)==1 & ~isnan(posx0));
        x0_center = posx0-track_center(1);
        y0_center = posy0-track_center(2);
        ang0 = nan(size(posx0));
        ind0=find(x0_center>0 & posx0~=0);
        ang0(ind0)=mod(atan(-y0_center(ind0)./x0_center(ind0)),2*pi);
        ind0=find(x0_center<0 & posx0~=0);
        ang0(ind0)=mod(atan(-y0_center(ind0)./x0_center(ind0))+pi,2*pi);
        
        xx0 = x0_center(ind);
        yy0 = y0_center(ind);
        radius00 = sqrt(xx0.^2 +yy0.^2);
        tt0 = t0(ind);
        ang00 = ang0(ind);
        
        %Fill in the gaps that the camera missed by circular interpolation
        radius0 = interp1(tt0,radius00,t0,'linear');
        
        ang00_unwrap = unwrap(ang00);
        ang0_unwrap_interp = interp1(tt0,ang00_unwrap,t0,'linear');
        ang0_interp = mod(ang0_unwrap_interp,2*pi);
        posx_nseg = radius0.*cos(ang0_interp);
        posy_nseg = radius0.*(-1).*sin(ang0_interp);
        posx_nseg(ind) = xx0;
        posy_nseg(ind) = yy0;
        ind1 = ind(1);
        posx_nseg(1:ind1-1) = posx_nseg(ind1);
        posy_nseg(1:ind1-1) = posy_nseg(ind1);
        ind2 = ind(end);
        posx_nseg(ind2+1:end) = posx_nseg(ind2);
        posy_nseg(ind2+1:end) = posy_nseg(ind2);
        posx_nseg = posx_nseg+track_center(1);
        posy_nseg = posy_nseg+track_center(2);
        
        post_nseg=t0;
        
        % Moving window mean filter
        for cc = 8:length(posx_nseg)-7
            posx_nseg(cc) = nanmean(posx_nseg(cc-7:cc+7));
            posy_nseg(cc) = nanmean(posy_nseg(cc-7:cc+7));
        end
        
        posx(ind_seg(nseg,1):ind_seg(nseg,2)) = posx_nseg;
        posy(ind_seg(nseg,1):ind_seg(nseg,2)) = posy_nseg;
        post(ind_seg(nseg,1):ind_seg(nseg,2)) = post_nseg;
    end
end

ind = find(~isnan(post));
posx = posx(ind);
posy = posy(ind);
post = post(ind);
	





% Removes position "jumps", i.e position samples that imply that the rat is
% moving quicker than physical possible. Samples in the "jump" parts are
% set to NaN
%
% (c) Raymond Skjerpeng, 2008.
function [x,y] = removeBadTracking_cz(x,y,threshold)

% Indexes to position samples that are to be removed
remInd = [];

diffX = diff(x);
diffY = diff(y);
diffR = sqrt(diffX.^2 + diffY.^2);
ind = find(diffR > threshold);
if isempty(ind)
    return
end

if ind(end) == length(x)
    offset = 2;
else
    offset = 1;
end

for ii = 1:length(ind)-offset
    if ind(ii+1) == ind(ii)+1
        % A single sample position jump, tracker jumps out one sample and
        % then jumps back to path on the next sample. Remove bad sample.
        remInd = [remInd; ind(ii)+1];
        ii = ii+1;
        continue
    else
        % Not a single jump. 2 possibilities:
        % 1. Tracker jumps out, and stay out at the same place for several
        % samples and then jumps back.
        % 2. Tracker just has a small jump before path continues as normal,
        % unknown reason for this. In latter case the samples are left
        % untouched.
        
        idx = find(x(ind(ii)+1:ind(ii+1))==x(ind(ii)+1));
        if length(idx) == length(x(ind(ii)+1:ind(ii+1)));
            remInd = [remInd; (ind(ii)+1:ind(ii+1))'];
        end
    end
end
% Remove the samples
x(remInd) = 0;
y(remInd) = 0;
