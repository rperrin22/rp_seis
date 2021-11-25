function rp_convert_seg2(filename)

%% read in file
%filename = 'Rec_00156.seg2';
[filehdr,filetxt,trctxt,data] = seg2read(filename, ...
        'want','filehdr,filetxt,trctxt,data');
    
%% create header vectors for segy output
% FFID
% Receiver Station
% Shot Station

% initialize vectors
FFID = zeros(length(trctxt),1);
Rec_stn = zeros(length(trctxt),1);
Shot_stn = zeros(length(trctxt),1);

% populate vectors
for count = 1:length(trctxt)
   FFID(count) = str2double(trctxt(count).shot_sequence_number);
   Rec_stn(count) = str2double(trctxt(count).receiver_station_number);
   Shot_stn(count) = str2double(trctxt(count).source_station_number);
end

%% write a segy
segyfilename = sprintf('%s.sgy',filename(1:end-5));

WriteSegy(segyfilename,data,'dt',.000125,'FieldRecord',FFID','SourceX',Shot_stn','GroupX',Rec_stn');

WriteSegyTraceHeaderValue(segyfilename,FFID,'key','FieldRecord');
WriteSegyTraceHeaderValue(segyfilename,Shot_stn,'key','SourceX');
WriteSegyTraceHeaderValue(segyfilename,Rec_stn,'key','GroupX');