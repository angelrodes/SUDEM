%% Shielding Under a Digital Elevation Model (SUDEM)
% SUDEM calculates the attenuation of cosmogenic radiation 
% at a set of locations (samples.xyz) 
% under a Digital Elevation Model (dem.xyz)

% Angel Rodes, 2024

clear
close all hidden
tic

%% Define parametres

% number of points for radial dems
n_azhimuts=100; % higher number will produce more detailed radial DEMs
n_distances=300; % higher number will produce more detailed radial DEMs
min_distance=0.01; % min distance in metres

randomized_models=1; % Number of randomized model used to estimate uncertainties

% Attenuationg lengths for cosgenic productions
L=[160 850 5000 500]; % g/cm2. See Rodés, Á. (2021) https://doi.org/10.3390/geosciences11090362

rho=2.65; % Density of the rock

display_3D_dem=1; % Set to 0 to save time when using big DEMs

%% Import Digital Elevation Model (all in metres)
% Load XYZ data from the file
data = importdata('dem.xyz');
% Extract x, y, and z coordinates
x = data(:, 1);
y = data(:, 2);
z = data(:, 3);


%% Import sample locations (all in metres)
% Load XYZ data from the file
data = importdata('samples.xyz');
% Extract x, y, and z coordinates
x_samples = data(:, 1);
y_samples = data(:, 2);
z_samples = data(:, 3);

% Choose how sample elevations are defined
% relative_sample_elevations is 0 for absolute (metres), and 1 for relative to ground (negative elevations in metres)
if max(z_samples)>0 % if there are positive elevations
    relative_sample_elevations=0; % sample elevations must be absolute (no samples floating in the air)
else % if elevations are 0 or below
    if min(z)>100 % and the base of the DEM is at least 100 m above sea level
        relative_sample_elevations=1; % sample elevations must be relative (change previous line if you want to calculate shieldings at depths below 100 m and below sea level usign absolute elevations)
    else % if not, ask user...
        % Create a yes/no dialog
        choice = questdlg('How are sample elevations defined?', ...
                  'Confirmation', ...
                  'Relative to ground','Absolute','Relative to ground');
              % Process the user's choice
              switch choice
                  case 'Relative to ground'
                      relative_sample_elevations=1;
                  case 'Absolute'
                      relative_sample_elevations=0; 
                  otherwise
                      relative_sample_elevations=0;
              end
    end
end
       

if relative_sample_elevations==1 % if relative elevations
    % convert to absolute elevations
    z_samples=z_samples+griddata(x,y,z,x_samples,y_samples);
end
kk
if sum(x_samples<min(x) | x_samples>max(x) | y_samples<min(y) | y_samples>max(y))>0
   warning([num2str(sum(x_samples<min(x) | x_samples>max(x) | y_samples<min(y) | y_samples>may(y))) ' samples out of DEM limits!']) 
end

%% Plot original DEM and samples
if display_3D_dem==1
    h=figure('units','normalized','outerposition',[0 0 0.5 0.5],'Name','DEM & samples','NumberTitle','off');
    
    hold on
    plot3(x,y,z,'.b')
    plot3(x_samples,y_samples,z_samples,'pr','MarkerFaceColor','r')
    for si=1:numel(x_samples)
        plot3(x_samples(si)*[1,1],y_samples(si)*[1,1],[z_samples(si) max(z(:))+10],'-r')
        plot3(x_samples(si),y_samples(si),max(z(:))+10,'vr')
    end
    xlabel('Lon.');
    ylabel('Lat.');
    zlabel('Elv.')
    grid on; box on;
    axis equal;
    view(-37.5+180,30)
    xlim([min(x) max(x)])
    ylim([min(y) max(y)])
    
    drawnow
end

%% Calculate shielding factors

SF_restults=zeros(numel(x_samples),numel(L)); % init. table
dSF_restults=zeros(numel(x_samples),numel(L)); % init. table of uncertainties

randomized_models=max(0,round(randomized_models)); % correct nonsense values of randomized_models

if randomized_models<2
    if randomized_models>0
        warning('Less than 2 randomized_models. Uncertainties will not be estimated properly.')
    else
        warning('No randomized_models. Uncertainties will not be estimated.')
    end
end

for sampleindex=1:numel(x_samples)
     disp(' ')
    disp(['Sample #' num2str(sampleindex)])
    disp(['----------'])
    x_sample=x_samples(sampleindex);
    y_sample=y_samples(sampleindex);
    z_sample=z_samples(sampleindex);
    
    SF_i=zeros(randomized_models+1,numel(L));
    
    for run_index=1:randomized_models+1 % run randomized_models+1 models
        
        % compute angles on actual DEM
        dist=abs(((x-x_sample).^2+(y-y_sample).^2).^0.5);
        z_dem=max(0,z-z_sample);
        angles=atan(z_dem./dist);
        
        % calculate distances from min_distance to max distance on DEM
        dd=logspace(log10(min_distance),log10(max(dist(angles>0))),n_distances);
        
        % azimuths
        azaz = linspace(0,1,n_azhimuts)*2*pi;
        azaz=azaz(1:end-1); % do not repeat azimuths (0 and 2*pi)
        
        if run_index>1 %randomize points
            dd=interp1(dd,[rand*min_distance,[1:numel(dd)-1]+rand(1,numel(dd)-1)]);
            azaz=interp1(azaz,[1:numel(azaz)-1]+rand(1,numel(azaz)-1));
        end
        
        % create radial DEM
        [az,dist] = meshgrid(azaz,dd);
        xq = x_sample+dist.*cos(az);
        yq = y_sample+dist.*sin(az);
        %     disp(['Total points:' num2str(numel(xq)) ' (similar to a ' num2str(ceil(numel(xq)^0.5)) 'x' num2str(ceil(numel(xq)^0.5)) ' grid)'])
        zq = griddata(x, y, z, xq, yq, 'natural'); % Interpolate the data onto the grid
        zq_extrap=griddata(x, y, z, xq, yq, 'nearest'); % Nearest Neighbour for extrapolations
        zq(isnan(zq))=zq_extrap(isnan(zq));
        DEM=zq;
        
        
        SF=0.*L; % init Shielding Factors
        
        for L_i=1:numel(L) % repeat for all attenuation lenghts
            
            % calculate Shielding Factor
            sample_elv=z_sample;
            z_dem=DEM-sample_elv;
            angles=atan((DEM-sample_elv)./dist);
            depths=(z_dem.^2+dist.^2).^0.5.*(angles>0); % compute no shielding for negative angles
            wehight1=sin(abs(angles)).^2.3; % weight for radiation angle
            angles_reference=linspace(atan(0/1),atan(1/0),100);
            [counts,angles_ref] = hist(abs(angles(angles>0)),angles_reference);
            wehight2=1./interp1(angles_ref,counts,abs(angles),'nearest'); % weight for angle distribution
            wehights=wehight1.*wehight2;
            attenuation_factors=(1-exp(-depths*100*rho/L(L_i)));
            shielding=1-sum(attenuation_factors(angles>0).*wehights(angles>0))/sum(wehights(angles>0));
            SF_i(run_index,L_i)=shielding;
            
            if run_index==1 && L_i==1 % plot only the first model, and after calculating the first SF
                % Plot distances vs. elevations
                h=figure('units','normalized',...
                    'outerposition',[0+0.4*sampleindex/numel(x_samples) 0.5-0.4*sampleindex/numel(x_samples) 0.5 0.5],...
                    'Name',['Sample ' num2str(sampleindex)],'NumberTitle','off');
                subplot(1,2,1)
                min_angle=20;
                hold on
                for az_i=unique(az(:))'
                    xplot=dist(az==az_i);
                    zplot=zq(az==az_i);
                    hue_value = (1 + cosd(mod(rad2deg(az_i), 360))) / 2; % Map azimuths to hue values 
                    color_az = hsv2rgb([hue_value, ones(size(hue_value)), ones(size(hue_value))]);
                    plot(xplot,zplot,'-','LineWidth',1,'Color',color_az)
                end
                plot(0,z_sample,'pk','MarkerFaceColor','w','MarkerSize',12,'LineWidth',2)
                xlim([0 max(dist(angles>min_angle/360*2*pi))])
                z_range=max(zq(angles>min_angle/360*2*pi))-z_sample;
                ylim([z_sample-z_range*0.1 max(zq(angles>min_angle/360*2*pi))+z_range*0.2])
                xlabel('Distance (m)')
                ylabel('Elevation (m)')
                title(['Sample #' num2str(sampleindex)])
                grid on
                box on
                
                subplot(1,2,2)
                min_angle=0;
                hold on
                for az_i=unique(az(:))'
                    xplot=dist(az==az_i);
                    zplot=zq(az==az_i);
                    hue_value = (1 + cosd(mod(rad2deg(az_i), 360))) / 2; % Map azimuths to hue values 
                    color_az = hsv2rgb([hue_value, ones(size(hue_value)), ones(size(hue_value))]);
                    plot(xplot,zplot,'-','LineWidth',1,'Color',color_az)
                end
                plot(0,z_sample,'pk','MarkerFaceColor','w','MarkerSize',12,'LineWidth',2)
                xlim([0 max(dist(angles>min_angle/360*2*pi))])
                z_range=max(zq(angles>min_angle/360*2*pi))-z_sample;
                ylim([z_sample-z_range*0.1 max(zq(angles>min_angle/360*2*pi))+z_range*0.2])
                xlabel('Distance (m)')
                ylabel('Elevation (m)')
                title(['SF_{(' num2str(L(L_i)) ' g/cm^{2})}=' num2str(shielding*100,3) '%'])
                grid on
                box on
                
                drawnow
                
            end
        end
    end
    
    SF=SF_i(1,:); % Shielding Factors from not randomized model
    if size(SF_i,1)>1
        dSF=std(SF_i,1); % Uncetainty based on randomized models
    else
       dSF=0*SF; % no uncertainty in case of one model only
    end
    

    % Display shielding factors in a table
    string_L=arrayfun(@(x) ['L_' num2str(x)], L, 'UniformOutput', false);
    T = array2table([SF;dSF;dSF./SF*100],'VariableNames',string_L,'RowName',{'Shielding Factors','Uncertainties','% uncertainties'});
    disp(T) 
    
    % Store sample data
    SF_restults(sampleindex,:)=SF;
    dSF_restults(sampleindex,:)=dSF;
    
    
end

%% Display Summary
disp(' ')
disp('Shielding Factors - summary')
disp('---------------------------')

if randomized_models>0
    SF_restults_table=arrayfun(@(m,d) [...
        num2str(100*round(m*10^(-floor(log10(d))+1))/10^(-floor(log10(d))+1),['%.' num2str(-floor(log10(d))-1) 'f'])...
        ' ± '...
        num2str(100*round(d*10^(-floor(log10(d))+1))/10^(-floor(log10(d))+1),['%.' num2str(-floor(log10(d))-1) 'f'])...
        ' %'],...
        SF_restults,dSF_restults, 'UniformOutput', false);
else
        SF_restults_table=arrayfun(@(m,d) [...
        num2str(100*m,4)...
        ' %'],...
        SF_restults,dSF_restults, 'UniformOutput', false);
end

% Display shielding factors in a table
string_L=arrayfun(@(x) ['L_' num2str(x)], L, 'UniformOutput', false);
string_samples=arrayfun(@(x) ['Sample #' num2str(x)], 1:numel(x_samples)', 'UniformOutput', false);

T = array2table(categorical(SF_restults_table),'VariableNames',string_L,'RowName',string_samples);
disp(T)

% Display time
disp(' ')
toc