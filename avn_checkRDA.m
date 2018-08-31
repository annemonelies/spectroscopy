function avn_checkRDA(rdaFile1,rdaFile2)
% function developDifference
% first_rdaFile = 'orig.rda';
% second_rdaFile = 'supp.rda';
% first_rdaFile = mot_c07_day2_scan1_orig.rda;
% second_rdaFile = mot_c07_day2_scan1_supp.rda;

% new_rdaFile = 'mot_diff.rda';
%%

noFile = false;
fid = [];
if nargin<1 || isempty(rdaFile1) || isempty(rdaFile2)
    
    [FILENAME,PATHNAME] = uigetfile('*.rda','Select rda files', 'MultiSelect','on');
    %     [FILENAME, PATHNAME] = uigetfile('*.rda');
    if ~iscell(FILENAME)
        if ischar(FILENAME)
            rdaFiles{1} = fullfile(PATHNAME,FILENAME);
            nrFiles = 1;
        else
            if isempty(FILENAME)||FILENAME
                noFile = true;
            end
        end
    end
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

if ~noFile

    for iFile = 1:nrFiles
        rdaFile = rdaFiles{iFile};
        
        [~,allFileNames{iFile},~] =fileparts(rdaFile);
        
        fid = fopen(rdaFile);
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
        
        fileData.complex{iFile} = complex_data;
        fileData.hmm{iFile} = hmm;
        fileData.hmm_complex{iFile} = hmm_complex;
        checkONOFF(iFile) = ONOFF;
        fclose(fid);
        
        timestamp = [storeTline1{22}(15:16),':',storeTline1{22}(17:18),':',storeTline1{22}(19:20)];
        newDisript = sprintf('%s %s %s',storeTline1{1}(end-10:end-2),timestamp,storeTline1{13}(15:26));
        filediscript{iFile} = strrep(newDisript,'_','-');
        
    end
    
    scanOptions = {'ONres','OFFres','difference'};
    FigHandle=figure; hold on
    figWidth = 500*nrFiles;
    set(FigHandle, 'Position', [100, 100, figWidth  ,500 ]);
    
    for iFile = 1: nrFiles
        if nrFiles > 1
            subplot(1,nrFiles,iFile);
        end
        plottable_data{iFile} = real(fft(fileData.hmm_complex{iFile}));
        plot(plottable_data{iFile})
        xlim([15 250])
        insertTitle = strrep(allFileNames{iFile},'_','-');
        title(sprintf('%s (%s)',insertTitle,scanOptions{checkONOFF(iFile)}))
        xlabel(filediscript{iFile})
    end
  
else
    fprintf('no file given, process aborted \n\n')
    
end


