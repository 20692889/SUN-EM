function [FEKO_Hfield] = parseFEKOhfefile(Const, hfe_filename)
    %parseFEKOefefile
    %   Date: 2018.05.24
    %   Usage:
    %       [Const, FEKO_data] = parseFEKOefefile(Const)
    %
    %   Input Arguments:
    %       Const: A global struct containing general program flow settings.
    %       efe_filename:
    %              The particular *.efe filename that will be used    
    %
    %   Output Arguments:
    %       FEKO_Efield:
    %           Struct containing the E-field sampling points and values.
    %
    %   Description:
    %       Reads in a FEKO *.efe file and extracts the electric field values,
    %       as well as sampling points.
    %
    %   =======================
    %   Written by Danie Ludick on 2018.05.24
    %   Stellenbosch University
    %   Email: dludick.sun.ac.za

    narginchk(2,2);

    fid = fopen(hfe_filename,'r');

    if fid == -1
        message_fc(Const,sprintf('Error reading FEKO *.hfe file: %s',hfe_filename));
        error(['Error reading FEKO *.hfe file: %s' hfe_filename]);
    end

    message_fc(Const,' ');
    message_fc(Const,...
        '------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Parsing the FEKO *.hfe file'));
    message_fc(Const,sprintf('  *.hfe file: %s',hfe_filename));

    % Initialise the return values.
    FEKO_Hfield = [];

    % ========================
    % Read the file type
    % ========================
    line=fgetl(fid);
    hfe_line_data = strsplit(line);
    % Make sure we are working with the correct data (Electric field)
    if (~strcmp(hfe_line_data{3},'Magnetic'))
        message_fc(Const,sprintf('Expecting magnetic field values'));
        error(['Expecting magnetic field values']);
    end%if

    % ========================
    % Read the file version
    % ========================
    line=fgetl(fid);
    hfe_line_data = strsplit(line);
    % Make sure we are working with the correct data (Magnetic field)
    if (str2num(hfe_line_data{3}) ~= Const.FEKO_efe_file_format)
        message_fc(Const,sprintf('Unsupported file format for *.hfe file'));
        error(['Unsupported file format for *.hfe file']);
    end%if
    
    % Skip a few lines until we see the "#Frequency line" below.
    frequency_line_found = false;
    while(~frequency_line_found)
        line=fgetl(fid);
        hfe_line_data = strsplit(line);
        if(strcmp(hfe_line_data{1},'#Frequency:'))
            frequency_line_found = true;
        end%if
    end%while (~frequency_line_found)
    
    % ========================
    % Read the Frequency
    % ========================
    % We have already split the frequency line into tokens
    %line=fgetl(fid);
    %efe_line_data = strsplit(line);
    FEKO_Hfield.frequency = str2double(hfe_line_data{2});

    % ========================
    % Read the co-ordinate system
    % ========================
    line=fgetl(fid);
    hfe_line_data = strsplit(line);
    FEKO_Hfield.coordinate_system = hfe_line_data{3};

    % ========================
    % Read the number of x samples
    % ========================

    line=fgetl(fid);
    efe_line_data = strsplit(line);
    FEKO_Hfield.number_x_samples = str2num(efe_line_data{5});

    % ========================
    % Read the number of y samples
    % ========================
    line=fgetl(fid);
    efe_line_data = strsplit(line);
    FEKO_Hfield.number_y_samples = str2num(efe_line_data{5});
    % ========================
    % Read the number of z samples
    % ========================
    line=fgetl(fid);
    efe_line_data = strsplit(line);
    FEKO_Hfield.number_z_samples = str2num(efe_line_data{5});

    % Skip a few lines again
    line=fgetl(fid);
    line=fgetl(fid);
    line=fgetl(fid);
    
    number_field_points = FEKO_Hfield.number_x_samples*FEKO_Hfield.number_y_samples*FEKO_Hfield.number_z_samples;

    % Initialise the field-point values
    % Note: The frequency axis has not yet been defined.

    % Cartesian co-ordinate (x, y, z)
    FEKO_Hfield.x_samples_m = zeros(number_field_points,1);
    FEKO_Hfield.y_samples_m = zeros(number_field_points,1);
    FEKO_Hfield.z_samples_m = zeros(number_field_points,1);

    % Spherical co-ordinate (r, theta, phi)
   % FEKO_Hfield.r_samples_m = zeros(number_field_points,1);
   % FEKO_Hfield.theta_samples_deg = zeros(number_field_points,1);
   % FEKO_Hfield.phi_samples_deg = zeros(number_field_points,1);
        
    % H-field value (Er, Etheta, Ephi)
    FEKO_Hfield.Hx = complex(zeros(number_field_points,1));
    FEKO_Hfield.Hy = complex(zeros(number_field_points,1));
    FEKO_Hfield.Hz = complex(zeros(number_field_points,1));

    for hfield_indx = 1:number_field_points
        line=fgetl(fid);
        hfe_line_data = strsplit(line);
        
        FEKO_Hfield.x_samples_m(hfield_indx,1) = str2double(hfe_line_data{2});
        FEKO_Hfield.y_samples_m(hfield_indx,1) = str2double(hfe_line_data{3});
        FEKO_Hfield.z_samples_m(hfield_indx,1) = str2double(hfe_line_data{4});
        
        FEKO_Hfield.Hx(hfield_indx,1) = str2double(hfe_line_data{5}) + 1i*str2double(hfe_line_data{6});
        FEKO_Hfield.Hy(hfield_indx,1) = str2double(hfe_line_data{7}) + 1i*str2double(hfe_line_data{8});
        FEKO_Hfield.Hz(hfield_indx,1) = str2double(hfe_line_data{9}) + 1i*str2double(hfe_line_data{10});

    end%while end_flag == 0