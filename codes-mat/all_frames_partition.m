% % % Brief summary of this function.
% % % Detailed explanation of this function.
% % % Detector independence - add to hist_all frames
% % % 
% % % imporvisation ideas:
% % % generate all the \Delta T bins and calculate everything in one go. (keep 2 levels though : fine and coarse)
% % % save as individual files for all the bins

function  [partition] = all_frames_partition(name_str,window_lb,window_ub,sourceflag,data_path,file_name,load_dat,current_path)

% tic
cd(data_path)

switch sourceflag 
    case 'monte'
        files = [file_name];
     
        gate_len = load_dat.export_vars.global_pars.TCSPC.n_photon_channels;
        MicroTime_Resolution = (4096*4/gate_len)*1E-3;
        
    case 'data'
        files = ls('*.csv');
        gate_len = 4096;
        MicroTime_Resolution = (4096*4/gate_len)*1E-3;
        
end


pathname = [pwd '\'];


MacroTime_Resolution = 25.6792;


MicroTime = (1:gate_len)*MicroTime_Resolution;

% IRF zero adjustment ?
% separate into detectors
% adjust IRF zero time shift 
% 
%% Define histogram parameters
loop_len = length(files(:,1));

TwoD_d1d1 = zeros(gate_len,gate_len);
TwoD_d2d2 = zeros(gate_len,gate_len);
TwoD_d1d2 = zeros(gate_len,gate_len);
TwoD_d2d1 = zeros(gate_len,gate_len);

OneD_dec1 = zeros(1,gate_len);
OneD_dec2 = zeros(1,gate_len);
dec1_all_traces = zeros(loop_len,gate_len);
dec2_all_traces = zeros(loop_len,gate_len);



for iter_files = 1:loop_len
    
    filename = files(iter_files,:); 
    spc = strfind(filename,'.csv');
    filename(spc+4:end) = [];
    
    fn = [pathname filename]
   
    n_gates  = gate_len;
%      wtl = waitbar(iter_files/loop_len);
%% time bins
    % units in seconds 
    
    del_lb = window_lb;
    del_ub = window_ub;

    fname = [name_str '_' num2str(iter_files) '_' num2str(del_lb) '_' num2str(del_ub)  '.mat' ];
    
%% 2DFD construction
switch sourceflag
    case 'monte'
        d1_macro_time = (load_dat.export_vars.measurement.D1.MacroTime);
        d1_gate       = (load_dat.export_vars.measurement.D1.Gate');
        d2_macro_time = (load_dat.export_vars.measurement.D2.MacroTime);
        d2_gate       = (load_dat.export_vars.measurement.D2.Gate');
        
    case 'data'
        A = readtable(fn);
        
        MacroTime = table2array(A(:,6));
        MacroTime = (MacroTime)*1E-9; % convert to seconds
        gate = table2array(A(:,9));
        gate = gate+1; % zero offset
        channel = table2array(A(:,3));  
        
        d1_macro_time = MacroTime(channel==1);
        d1_gate = gate(channel==1);
        d2_macro_time = MacroTime(channel==2);
        d2_gate = gate(channel==2);
    
 end
    
    
    cd(current_path)
    mapd1d1 = partition_2D_core(n_gates,del_lb,del_ub,d1_macro_time,d1_gate,d1_macro_time,d1_gate);
    disp '11'
    mapd1d2 = partition_2D_core(n_gates,del_lb,del_ub,d1_macro_time,d1_gate,d2_macro_time,d2_gate);
    disp '12'
    mapd2d1 = partition_2D_core(n_gates,del_lb,del_ub,d2_macro_time,d2_gate,d1_macro_time,d1_gate);
    disp '21'
    mapd2d2 = partition_2D_core(n_gates,del_lb,del_ub,d2_macro_time,d2_gate,d2_macro_time,d2_gate);
    disp '22'

    TwoD_d1d1 = TwoD_d1d1 + mapd1d1;
    TwoD_d1d2 = TwoD_d1d2 + mapd1d2;
    TwoD_d2d2 = TwoD_d2d2 + mapd2d2;
    TwoD_d2d1 = TwoD_d2d1 + mapd2d1;
%     toc


% % % % extract long dec
    bin_edges = 1:gate_len+1; 
    tdec1 = histcounts(d1_gate,bin_edges);
    dec1_all_traces(iter_files,:) = tdec1;
    OneD_dec1 = OneD_dec1 + tdec1;
    
    tdec2 = histcounts(d2_gate,bin_edges);
    dec2_all_traces(iter_files,:) = tdec2;
    OneD_dec2 = OneD_dec2 + tdec2;

end

%% wrapup
%      close(wtl)
     time_stamp = datestr(now);
     
%      save(fname,'TwoD_d1d1','TwoD_d2d2',"TwoD_d2d1","TwoD_d1d2",...
%          'OneD_dec1','OneD_dec2', 'dec2_all_traces','dec1_all_traces')
%      

     partition.fname = fname;
     partition.gate_len = gate_len;
     partition.TwoD_d1d1 = TwoD_d1d1;
     partition.TwoD_d2d2=TwoD_d2d2;
     partition.TwoD_d2d1=TwoD_d2d1;
     partition.TwoD_d1d2=TwoD_d1d2;
     partition.OneD_dec1=OneD_dec1;
     partition.OneD_dec2=OneD_dec2;
     partition.dec2_all_traces=dec2_all_traces;
     partition.dec1_all_traces=dec1_all_traces;
     
    
%     timeeee = toc
    
end 



