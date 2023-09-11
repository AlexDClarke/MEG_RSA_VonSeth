function D = spm_eeg_remove_evoked(S)
% Removes evoked signal from single trial M/EEG data
% FORMAT D = spm_eeg_zscore(S)
%
% S        - optional input struct
%      fields of S:
%   S.D       - MEEG object or filename of M/EEG mat-file with epoched data
%   S.save    - save the baseline corrected data in a separate file [default: true]
%   S.updatehistory - update history information [default: true]
%   S.prefix     - prefix for the output file (default - 'b')
%
% D        - MEEG object (also saved on disk if requested)
%__________________________________________________________________________
%
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% AC based on spm_eeg_bc.m 5212 2013-01-26 13:16:36Z vladimir $

SVNrev = '$Rev: 5212 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG remove evoked'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
D = spm_eeg_load(S.D);

%-Get input parameters
%--------------------------------------------------------------------------
if ~isfield(S, 'prefix'),          S.prefix = 'n';           end
if ~isfield(S, 'save'),            S.save = 1;               end
if ~isfield(S, 'updatehistory'),   S.updatehistory = 1;      end

%-Remove evoked from trials
indchannels = [D.indchantype('Filtered') D.indchantype('MEGCOMB')];
if S.save
    D = D.copy([S.prefix D.fname]);
end

tmp = mean(D(indchannels, :, :),3);  % evoked signals
D(indchannels,:,:) = (D(indchannels,:,:) - repmat(tmp,[1,1,D.ntrials])); % remove evoked

%-Save data
%--------------------------------------------------------------------------
if ~isfield(S, 'updatehistory') || S.updatehistory   
    D = D.history(mfilename, S);
    save(D);
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG remove evoked: done'); spm('Pointer','Arrow');
