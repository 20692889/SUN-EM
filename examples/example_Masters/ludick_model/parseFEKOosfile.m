function [FEKO_currents] = parseFEKOosfile (Const, os_filename)

    %parseFEKOosfile
    %   Date: 2021.03.07
    %   Usage:
    %       [FEKO_currents] = parseFEKOosfile (Const, os_filename)
    %
    %   Input Arguments:
    %       Const: A global struct containing general program flow settings.
    %       os_filename:
    %              The particular *.os/*.ol filename that will be used    
    %
    %   Output Arguments:
    %       FEKO_currents:
    %           Struct containing the currents sampling points and values.
    %
    %   Description:
    %       Reads in a FEKO *.os/*ol file and extracts the surface current values,
    %       as well as sampling points.
    %
    %   =======================
    %   Written by Willem de la Bat on 2021.03.07
    %   Stellenbosch University
    %   Email: 20692889.sun.ac.za

    narginchk(2,2);

    fid = fopen(os_filename,'r');

    if fid == -1
        message_fc(Const,sprintf('Error reading FEKO *.os file: %s',os_filename));
        error(['Error reading FEKO *.os file: %s' os_filename]);
    end

    message_fc(Const,' ');
    message_fc(Const,...
        '------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Parsing the FEKO *.os file'));
    message_fc(Const,sprintf('  *.os file: %s',os_filename));

    % Initialise the return values.
    FEKO_currents = [];

    %FEKO file version number supported -- from sunem_init
    
        Const.FEKO_os_file_format = 6; % FEKO *.os file format

    % ========================
    % Read the file type
    % ========================
    line=fgetl(fid);
    os_line_data = strsplit(line);
    % Make sure we are working with the correct data (Currents)
    if (~strcmp(os_line_data{3},'Currents'))
        message_fc(Const,sprintf('Expecting current values'));
        error(['Expecting current values']);
    end%if

    % ========================
    % Read the file version
    % ========================
    line=fgetl(fid);
    os_line_data = strsplit(line);
    % Make sure we are working with the correct data (Currents)
    if (str2num(os_line_data{3}) ~= Const.FEKO_os_file_format)
        message_fc(Const,sprintf('Unsupported file format for *.os file'));
        error(['Unsupported file format for *.os file']);
    end%if
    
    % Skip a few lines until we see the "#Frequency line" below.
    frequency_line_found = false;
    while(~frequency_line_found)
        line=fgetl(fid);
        os_line_data = strsplit(line);
        if(strcmp(os_line_data{1},'#Frequency:'))
            frequency_line_found = true;
        end%if
    end%while (~frequency_line_found)
    
    % ========================
    % Read the Frequency
    % ========================
    % We have already split the frequency line into tokens
    %line=fgetl(fid);
    %efe_line_data = strsplit(line);
    FEKO_currents.frequency = str2double(os_line_data{2});

    % ========================
   % Read No. of Electric Current Triangle Samples
    % ========================
   
     line=fgetl(fid);
     os_line_data = strsplit(line);
     FEKO_currents.number_Ecurr_triangles = str2num(os_line_data{7});

    % ========================
    % Read the number of R samples
    % ========================
%     line=fgetl(fid);
%     os_line_data = strsplit(line);
%     FEKO_currents.number_r_samples = str2num(os_line_data{5});

    % ========================
    % Read the number of Theta samples
    % ========================
%     line=fgetl(fid);
%     os_line_data = strsplit(line);
%     FEKO_currents.number_theta_samples = str2num(os_line_data{5});
% 
%     % ========================
    % Read the number of Phi samples
    % ========================
%     line=fgetl(fid);
%     os_line_data = strsplit(line);
%     FEKO_currents.number_phi_samples = str2num(os_line_data{5});

    % Skip a few lines again
    line=fgetl(fid);
    line=fgetl(fid);
    
    %number_field_points = FEKO_currents.number_r_samples*FEKO_currents.number_theta_samples*FEKO_currents.number_phi_samples;

    % Initialise the field-point values
    % Note: The frequency axis has not yet been defined.
    number_Ecurr_triangles = FEKO_currents.number_Ecurr_triangles;
    % Cartesian co-ordinate (x, y, z)
    FEKO_currents.x_samples_m = zeros(number_Ecurr_triangles,1);
    FEKO_currents.z_samples_m = zeros(number_Ecurr_triangles,1);
    FEKO_currents.y_samples_m = zeros(number_Ecurr_triangles,1);
        
    % Current value (Jx, Jy, Jz)
    FEKO_currents.Jx = complex(zeros(number_Ecurr_triangles,1));
    FEKO_currents.Jy = complex(zeros(number_Ecurr_triangles,1));
    FEKO_currents.Jz = complex(zeros(number_Ecurr_triangles,1));

    for current_indx = 1:number_Ecurr_triangles
        line=fgetl(fid);
        os_line_data = strsplit(line);
        
        FEKO_currents.x_samples_m(current_indx,1) = str2double(os_line_data{2});
        FEKO_currents.y_samples_m(current_indx,1) = str2double(os_line_data{3});
        FEKO_currents.z_samples_m(current_indx,1) = str2double(os_line_data{4});
        
        FEKO_currents.Jx(current_indx,1) = str2double(os_line_data{5}) + 1i*str2double(os_line_data{6});
        FEKO_currents.Jy(current_indx,1) = str2double(os_line_data{7}) + 1i*str2double(os_line_data{8});
        FEKO_currents.Jz(current_indx,1) = str2double(os_line_data{9}) + 1i*str2double(os_line_data{10});

    end%while end_flag == 0
    
    