function avn_createMissingRDA(rdaFile1,rdaFile2)
% function developDifference
% first_rdaFile = 'orig.rda';
% second_rdaFile = 'supp.rda';
% first_rdaFile = mot_c07_day2_scan1_orig.rda;
% second_rdaFile = mot_c07_day2_scan1_supp.rda;

% new_rdaFile = 'mot_diff.rda';
%%



if nargin<1 || isempty(rdaFile1) || isempty(rdaFile2)
    
    [FILENAME,PATHNAME] = uigetfile('*.rda','Select rda files', 'MultiSelect','on');
    %     [FILENAME, PATHNAME] = uigetfile('*.rda');
end
if iscell(FILENAME)
    nrFiles = numel(FILENAME);
    
    for iFile = 1:nrFiles
        if iscell(PATHNAME)
            rdaFiles{iFile} = fullfile(PATHNAME{iFile},FILENAME{iFile});
        else
            rdaFiles{iFile} = fullfile(PATHNAME,FILENAME{iFile});
        end
    end
else
    rdaFiles{1} = fullfile(PATHNAME,FILENAME);
    nrFiles = 1;
end

for iFile = 1:nrFiles
    rdaFile = rdaFiles{iFile};
    fid = fopen(rdaFile);
    
    
    if isempty(fid)
        return
    end
    
    tline = fgets(fid);
    
    head_start_text = '>>> Begin of header <<<';
    head_end_text   = '>>> End of header <<<';
    iCount = 1;
    while (isempty(strfind(tline , head_end_text)))
        
        tline = fgets(fid);
        if ( isempty(strfind (tline , head_start_text)) + isempty(strfind (tline , head_end_text )) == 2)
            
            % Store this data in the appropriate format
            storeTline1{iCount} = tline;
            occurence_of_colon = findstr(':',tline);
            variable = tline(1:occurence_of_colon-1) ;
            value    = tline(occurence_of_colon+1 : length(tline)) ;
            
            switch variable
                case {'InstanceComments'  }
                    idxON = regexp(value,'On');
                    idxOFF = regexp(value,'Off');
                    idxDIFF = regexp(value,'Diff');
                    
                    if ~isempty(idxON)
                        ONOFF = 1;
                    elseif ~isempty(idxOFF)
                        ONOFF = 2;
                    elseif ~isempty(idxDIFF)
                        ONOFF = 3;
                    else
                        ONOFF = 0;
                    end
                    
                case {'VectorSize' }
                    %Integers
                    eval(['rda.' , variable , ' = str2num(value); ']);
                case {'CSIMatrixSize[0]' }
                    rda.CSIMatrix_Size(1) = str2num(value);
                case {'CSIMatrixSize[1]' }
                    rda.CSIMatrix_Size(2) = str2num(value);
                case {'CSIMatrixSize[2]' }
                    rda.CSIMatrix_Size(3) = str2num(value);
                    
                otherwise
                    % just store the value otherwise
            end
            
        else
            % Don't bother storing this bit of the output
        end
        
        iCount = iCount+1;
    end
    
    % Siemens documentation suggests that the data should be in a double complex format (8bytes for real, and 8 for imaginary?)
    
    complex_data = fread(fid , rda.CSIMatrix_Size(1) * rda.CSIMatrix_Size(1) *rda.CSIMatrix_Size(1) *rda.VectorSize * 2 , 'double');
    hmm = reshape(complex_data,  2 , rda.VectorSize , rda.CSIMatrix_Size(1) ,  rda.CSIMatrix_Size(2) ,  rda.CSIMatrix_Size(3) );
    hmm_complex = complex(hmm(1,:,:,:,:),hmm(2,:,:,:,:));
    
    fileData.complex{ONOFF} = complex_data;
    fileData.hmm{ONOFF} = hmm;
    fileData.hmm_complex{ONOFF} = hmm_complex;
    checkONOFF(iFile) = ONOFF;
    
    fclose(fid);
end


scanOptions = {'ONres','OFFres','difference'};
scanOptionsNamer = {'onrs','offr','difr'};
scanOptionsSiemens = {'On Resonance','Off Resonance','Difference'};

missingData = setdiff([1,2,3],checkONOFF);

keyboard;
%% only one file
if numel(missingData)>1
    fprintf('only one File to display, no files restored \n')
    
    figure; hold on
    for iFile = 1;
        plottable_data{iFile} = real(fft(fileData.hmm_complex{checkONOFF(iFile)}));
        plot(plottable_data{iFile})
        xlim([15 250])
    end
else
    
    %% two or three files
    if ~isempty(missingData)
        switch missingData
            case 1
                fileData.complex{missingData} = fileData.complex{2}-fileData.complex{3};
            case 2
                fileData.complex{missingData} = fileData.complex{1}-fileData.complex{3};
            case 3
                fileData.complex{missingData} = fileData.complex{1}-fileData.complex{2};
        end
        fileData.hmm{missingData} = reshape(fileData.complex{missingData},  2 , rda.VectorSize , rda.CSIMatrix_Size(1) ,  rda.CSIMatrix_Size(2) ,  rda.CSIMatrix_Size(3) );
        fileData.hmm_complex{missingData}  = complex(fileData.hmm{missingData}(1,:,:,:,:),fileData.hmm{missingData}(2,:,:,:,:));
        
        %% write missing .rda file
        newRDAfile = [rdaFiles{1}(1:end-8) sprintf('%s.rda',scanOptionsNamer{missingData})];
        
        fileID = fopen(newRDAfile, 'w');
        
        % write start header
        tlines = sprintf('%s\r',head_start_text);
        fwrite(fileID,tlines);
        
        % write other header points, replace relevant ones
        for iPart = 1:numel(storeTline1)
            
            tlines = storeTline1{iPart};
            occurence_of_colon = findstr(':',tlines);
            variable = tlines(1:occurence_of_colon-1) ;
            value    = tlines(occurence_of_colon+1 : length(tlines)) ;
            
            switch variable
                case {'InstanceNumber' }
                    idx = regexp(tlines,value);
                    tlines(idx+1) = '3';
                case {'InstanceComments' }
                    idx = regexp(tlines,'On Resonance|Difference|Off Resonance');
                    tlines = [tlines(1:idx-1) sprintf('%s',scanOptionsSiemens{missingData}) tlines(end)];
            end
            fwrite(fileID,tlines);
        end
        
        % write end header
        tlines = head_end_text;
        fwrite(fileID,tline);
        
        
        % write new complex numbers
        fwrite(fileID,fileData.complex{missingData},'double');
        fclose(fileID);
        
    end
    
    FigHandle=figure; hold on
    set(FigHandle, 'Position', [100, 100, 1400  ,500 ]);
    for iFile = 1:3;
        
        plottable_data{iFile} = real(fft(fileData.hmm_complex{iFile}));
        subplot(1,3,iFile);
        plot(plottable_data{iFile})
        xlim([15 250])
        if iFile==missingData
            title(sprintf('%s (RESTORED)',scanOptions{iFile}))
        else
            title(sprintf('%s',scanOptions{iFile}))
        end
    end
    
end






