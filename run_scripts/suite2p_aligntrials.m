%% Folder

[basedata baseanalysis] = DataFolders;

datafolder = uigetdir(basedata, 'Enter raw data folder.');

analysisfolder = uigetdir(baseanalysis, 'Enter analysis folder.');

[data_path, data_name, ~] = fileparts(datafolder);
[analysis_path, analysis_name, ~] = fileparts(analysisfolder);

if ~strcmp(data_name, analysis_name)
   error 'Raw data and analysis folders don''t match'
end

finfo = fullfile(analysisfolder,'ExperimentInfo.mat');

% if exist('DataPreprocessScript_restart_','var')paid
%     % DataPreprocessScript was restarted by command line inside
%     % DataPreprocessScript itself, folder is already defined
%     clear DataPreprocessScript_restart_
% elseif exist('','var') && any(strfind(finfo,baseanalysis))
%     fullfolder = uigetdir(fileparts(fileparts(finfo)));
% else
%     fullfolder = uigetdir(baseanalysis);
% end
% 
% if ~strfind(baseanalysis,fullfolder), error 'data folder must be inside ''ANALYSIS''', end
% folder = fullfolder(length(baseanalysis)+1:end);
% basefolder = fn_fileparts(folder,'base');
% if isempty(basefolder), return, end
% date = basefolder(1:6);
% if isempty(regexp(date,'^1\d[01]\d[0-3]\d$')) %#ok<RGXP1>
%     error 'folder name must start with date in format ''yymmdd''.'
% end
% 
% disp(folder)
% datafolder = fullfile(basedata,folder);
% analysisfolder = fullfile(baseanalysis,folder);
% if ~exist(datafolder,'dir')
%     error('no corresponding folder ''%s'' in ''RAW DATA''',folder)
% end



%% Check if any registration result already exists
% 
% douseprev = false;
% 
% finfo = fullfile(analysisfolder,'ExperimentInfo.mat');
% if exist(finfo,'file')
%         hmesc = fn_loadvar(finfo,'hmesc');
%         nunits = size(hmesc);
%         nunits = nunits(2);
%         answer = questdlg(sprintf('Manual alignment already done with %d selected trials, what do you want to do?', nunits),...
%             'Restart from scratch?',...
%             'Keep them and run suite2p', ...
%             'Manually align trials again',...
%             'Cancel');
%         if isempty(answer), disp 'interrupted', return, end
%         if strcmp(answer, 'Keep them and run suite2p');
%             douseprev = true;
%         end
% end

%% Preprocessing settings


%% Scan for MESc files - get time of every trial

disp 'reading MESc headers'

d = dir(fullfile(datafolder,'*.mesc'));
[dum ord] = sort([d.datenum]); %#ok<ASGLU>
fmesc = {d(ord).name};
if ~isscalar(fmesc)
    answer = questdlg(['Please check the order of MESc files:', fmesc], '', ...
        'OK', 'Change', 'OK');
    if isempty(answer), disp 'interrupted', return, end
    if strcmp(answer, 'Change')
        fmesc = fmesc(fn_list(fmesc));
    end
end

fmescfull = fn_map(@(s)fullfile(datafolder,s),fmesc);
sperday = 60*60*24;
%%
nfmesc = length(fmesc);
hmesc = cell(1,nfmesc);
for i=1:nfmesc
    [h desc] = mesc_header(fmescfull{i});
    units = [h.sessions.Units];
    trials = (1:length(units));

    if isempty(fieldnames(units))
        fprintf('file ''%s'': no trial could be read correctly\n',fmescfull{i}) 
        bad = true(1,length(units));
    else
        t = {units.MeasurementDatePosix};
        bad = fn_isemptyc(t);
    end
    units(bad) = [];
    desc(bad,:) = [];
    trials(bad) = [];
    nunit = length(units);

    if ~all(bad)
        t = double([units.MeasurementDatePosix])/sperday + datenum('01Jan1970 02:00:00');
    else
        t = zeros(1,0);
    end
    hmesc{i} = struct( ...
        'file',     i, ...
        'trial',    num2cell(trials), ...
        'time',     num2cell(t), ...
        'timestr',  cellstr(datestr(t))', ...
        'desc',     cellstr(desc)');
end
hmesc = [hmesc{:}];

%% Scan for DAT files - get time of every trial

disp 'reading Elphy files'

d = dir(fullfile(analysisfolder,'*.DAT'));
[dum ord] = sort([d.datenum]);
fdat = {d(ord).name};
fdatfull = fn_map(@(s)fullfile(analysisfolder,s),fdat);

nfdat = length(fdat);
hdat = cell(1,nfdat);
for i=1:nfdat
    [r vectors mpar xpar] = elphy_read(fdatfull{i}); 
    if isempty(fieldnames(xpar)), xpar = struct('fix',mpar); end
    ntrial = size(r,2);
    nstimperep = xpar.fix.NStimPerEp;
    try
        trecord = reshape(vectors.trecord(1:ntrial*nstimperep),nstimperep,ntrial);
    catch
        try
            trecord = reshape(vectors.iCond(1:ntrial*nstimperep),nstimperep,ntrial);
        catch
            trecord = reshape(vectors.tprebuild(1:ntrial*nstimperep),nstimperep,ntrial);
        end
    end
    conds = cell(1,ntrial);
    for k=1:ntrial, conds{k}=trecord(:,k)'; end
    hdat{i} = struct( ...
        'file',     i, ...
        'filename', fdat{i}, ...
        'trial',    num2cell(1:ntrial), ...
        'timestr',  vectors.dates, ...
        'time',     num2cell(datenum(vectors.dates)'), ...
        'cond',     conds, ...
        'fullpar',  xpar, ...
        'sound',[],'image',[]);
    if isfield(xpar,'soundlist') && ~isempty(xpar.soundlist)
        try
            for k=1:ntrial, hdat{i}(k).sound = xpar.soundlist(conds{k}); end
        catch
            disp 'something went wrong with the sound list'
        end
    end
    if isfield(xpar,'imagelist') && ~isempty(xpar.imagelist)
        try
            for k=1:ntrial, hdat{i}(k).image = xpar.imagelist(conds{k}); end
        catch
            disp 'something went wrong with the image list'
        end
    end
end
hdat = [hdat{:}];
if isempty(hdat)
    answer = questdlg('No Elphy data found. Was there no stimulation?','DataPreprocessScript','Indeed','Cancel','Indeed');
    if ~strcmp(answer,'Indeed'), disp 'interrupted', return, end
    xpar = [];
end


%% Manual alignment of trials

seq = [];
hmesc0 = hmesc; hdat0 = hdat;
A = cell(1,length(hmesc));
for i=1:length(hmesc)
    hi = hmesc(i);
    A{i} = sprintf('%s - %s - %s',fmesc{hi.file},hi.desc,hi.timestr);
end
B = cell(1,length(hdat));
for i=1:length(hdat)
    hi = hdat(i);
    B{i} = sprintf('%s - trial %i - %s',fdat{hi.file},hi.trial,hi.timestr);
end

if isempty(hdat)
    
    % Align mesc trial with binary stimulation sequence (Anthony's protocol).   
    seq_file = fullfile(datafolder, [data_name, '_sequence.bin']);
    
    if ~isfile(seq_file)
        answer = questdlg('No binary stimulation sequence found (Anthony). Was there no stimulation?','Stimulation sequence','Indeed','Cancel','Indeed');
        if ~strcmp(answer,'Indeed'), disp 'interrupted', return, end
    
        disp 'select trials'
        deltaA = ([hmesc0.time] - [hmesc0([1 1:end-1]).time])*sperday;
        A1 = A; for i=1:length(A), A1{i}=[A{i} sprintf(' (+%.0fs)',deltaA(i))]; end
        indices = fn_listedit(A1);
        hmesc = hmesc0(indices);
    else
        binary_sequence = true;
        f = fopen(seq_file);
        seq = fread(f, 'uint8');
        fclose(f);
        
        disp 'alignment between mesc and binary sequence'
        deltaA = ([hmesc0.time] - [hmesc0([1 1:end-1]).time])*sperday;
        A1 = A; for i=1:length(A), A1{i}=[A{i} sprintf(' (+%.0fs)',deltaA(i))]; end
        
        indices = fn_listedit(A1, seq);
        hmesc = hmesc0(indices(:,1));
        seq = seq(indices(:,2));
    end
    
else    

    disp 'alignment between MESc and Elphy'

    % Load correspondences / Try to guess
    fcorres = fullfile(analysisfolder,'MESc-Elphy correspondences.mat');
    if exist(fcorres,'file')
        indices = fn_loadvar(fcorres);
        sperday = 60*60*24;
        deltat = ([hmesc0(indices(:,1)).time]-[hdat0(indices(:,2)).time])*sperday;
        fprintf('(correspondence file already exists; time differences between\nMESc and Elphy are comprised between %.0f and %.0fs)\n',min(deltat),max(deltat));
    else
        if length(hmesc0)==length(hdat0)
            % try to guess
            indices = repmat((1:length(hmesc0))',1,2);
        else
            indices = [];
        end
        okcorr = false;
        while ~okcorr
            if ~isempty(indices)
                deltat = ([hmesc0(indices(:,1)).time]-[hdat0(indices(:,2)).time])*sperday;
                fulllist = [char(A(indices(:,1))) repmat(' - ',size(indices,1),1) char(B(indices(:,2))) num2str(deltat',' - diff=%.0fs')];
                figure(1), clf
                uicontrol('units','normalized','pos',[0 0 1 1],'style','listbox','string',fulllist);
                sperday = 60*60*24;
                hf = msgbox( ...
                    sprintf('Time differences between MESc and Elphy are comprised between %.0f and %.0fs. Accept? (close figure for not accepting)',min(deltat),max(deltat)));
                hu = findall(hf,'type','uicontrol');
                set(hu,'callback',@(u,e)delete(u));
                waitfor(hu)
                okcorr = ishandle(hf);
                if okcorr
                    delete(hf)
                end
                if ishandle(1), close(1), end
                if ~okcorr
                    answer = questdlg('What do you want to do?', '', ...
                        'Align MESc/Elphy again', 'Reorder MESc files', 'Interrupt', 'Align MESc/Elphy again');
                    switch answer
                        case {'', 'Interrupt'}
                            disp interrupted
                            return
                        case 'Reorder MESc files'
                            disp 'restarting DataPreprocessScript'
                            DataPreprocessScript_restart_ = true;
                            DataPreprocessScript
                            return
                    end
                end
            end
            if ~okcorr
                deltaA = ([hmesc0.time] - [hmesc0([1 1:end-1]).time])*sperday;
                A1 = A; for i=1:length(A), A1{i}=[A{i} sprintf(' (+%.0fs)',deltaA(i))]; end
                deltaB = ([hdat0.time] - [hdat0([1 1:end-1]).time])*sperday;
                B1 = B; for i=1:length(B), B1{i}=[B{i} sprintf(' (+%.0fs)',deltaB(i))]; end
                indices = fn_listedit(A1,B1);
            end
        end
        fn_savevar(fcorres,indices)
    end

    hmesc = hmesc0(indices(:,1));
    hdat = hdat0(indices(:,2));
end

%% Same protocol for every Elphy file or not?

if nfdat>1

    % Get protocol number for each Elphy file
    nconds = zeros(1,nfdat);
    s = struct;
    spec = struct('lab',{'label','Please indicate a protocol number for each Elphy file'});    
    for k=1:nfdat
        idx = find([hdat.file]==k);
        if isempty(idx), continue, end % data from this file is not used
        % number of conditions for this Elphy file
        nconds(k) = length(hdat(idx(1)).fullpar.table);
        if any([hdat(idx).cond]>nconds(k)), error 'condition number exceeds number of conditions!', end
        % prepare structure that user will edit
        f = ['E' num2str(k)];
        s.(f) = 1;
        spec(1).(f) = 'stepper 1 1 Inf';
        spec(2).(f) = fdat{k};
    end
    s = fn_structedit(s,spec);
    if isempty(s), disp 'interrupted', return, end
    protocol = zeros(1,nfdat);
    for k=1:nfdat
        if ~any([hdat.file]==k), continue, end % data from this file is not used
        f = ['E' num2str(k)];
        protocol(k) = s.(f);
    end

    % Assign global condition numbers
    nprotocol = max(protocol);
    protocoloffset = zeros(1,nprotocol);
    for kp=1:nprotocol-1
        kdat = find(protocol==kp);
        if isempty(kdat)
            answer = inputdlg( ...
                ['In order to set global condition numbers, please indicate the number of conditions of missing protocol ' num2str(kp)], ...
                'DataPreprocessScript',1,{'0'});
            protocoloffset(kp+1) = protocoloffset(kp) + str2double(answer{1});
        elseif any(diff(nconds(kdat)))
            error(['Files ' fn_strcat(fdat(kdat),', ') ' have different numbers of conditions, therefore cannot belong to the same protocol'])
        else
            protocoloffset(kp+1) = protocoloffset(kp) + nconds(kdat(1));
        end
    end
    for k=1:length(hdat)
        hdat(k).protocol = protocol(hdat(k).file);
        hdat(k).condglob = protocoloffset(protocol(hdat(k).file)) + hdat(k).cond;
    end
elseif ~isempty(hdat)

    [hdat.protocol] = deal(1);
    [hdat.condglob] = deal(hdat.cond);

end


%% Make list of keys to read from mesc files and retrieve dt and zdim.

nloops = length(hmesc);
key_count = 1;
keys = {};

for iloop=1:nloops
    % Solve session and unit indices.
    session = hmesc(iloop).desc;
    session = session(8);
    unit = strsplit(hmesc(iloop).desc);
    unit = unit{1};
    unit = unit(14:end);
    keys{key_count} = sprintf('/MSession_%s/MUnit_%s', session, unit);
    % Retrieve dt and zdim for each loop.
    hmesc(iloop).zdim = h5readatt(fmescfull{hmesc(iloop).file}, keys{key_count}, 'ZDim');
    hmesc(iloop).dt = h5readatt(fmescfull{hmesc(iloop).file}, keys{key_count}, 'ZAxisConversionConversionLinearScale');
    key_count = key_count + 1;
end


%% Save important info for suite2p

disp 'Saving info for suite2p.'
[filepath, name, ext] = fileparts(finfo);
finfo_suite2p = [name, '_suite2p', ext];
finfo_suite2p = fullfile(filepath, finfo_suite2p);

fn_savevar(finfo_suite2p, analysisfolder, fmescfull, hmesc, hdat, keys, xpar, seq)

disp('Done, you can run suite2p!')
