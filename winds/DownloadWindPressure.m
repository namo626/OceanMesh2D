%% Download wind and pressure files for specific months

% clear all
% close all
% clc

storms = [];
% storms{end+1}.name = 'Katrina2005';
% storms{end}.days = '08/23/2005-08/31/2005';
% 
% storms{end+1}.name = 'Rita2005';
% storms{end}.days = '09/18/2005-09/26/2005';
% 
storms{end+1}.name = 'Ike2008';
storms{end}.days = '08/06/2008-09/5/2008';

% storms{end+1}.name = 'Gustav2008';
% storms{end}.days = '08/25/2008-09/05/2008';
% storms{end}.days = '08/01/2008-09/05/2008';

% storms{end+1}.name = 'Irene2011';
% storms{end}.days = '08/21/2011-08/30/2011';
% storms{end}.days = '08/01/2011-09/30/2011';
% 
% storms{end+1}.name = 'Isaac2012';
% storms{end}.days = '08/21/2012-09/03/2012';
% 
% storms{end+1}.name = 'Sandy2012';
% storms{end}.days = '10/22/2012-11/02/2012';
% storms{end}.days = '10/01/2012-11/02/2012';
% 
% storms{end+1}.name = 'Matthew2016';
% storms{end}.days = '09/28/2016-10/10/2016';
% 
% storms{end+1}.name = 'Irma2017';
% storms{end}.days = '08/30/2017-09/12/2017';
% storms{end}.days = '08/10/2017-09/12/2017';
% 
% storms{end+1}.name = 'Maria2017';
% storms{end}.days = '09/16/2017-10/02/2017';
% 
% storms{end+1}.name = 'Harvey2017';
% storms{end}.days = '08/17/2017-09/03/2017';

% storms{end+1}.name = 'Matthew2016';
% storms{end}.days = '09/13/2016-10/11/2016';

% storms{end+1}.name = 'Harvey2017';
% storms{end}.days = '07/26/2017-09/15/2017';

% storms{end+1}.name = 'Florence2018';
% storms{end}.days = '08/23/2018-09/26/2018';

% storms{end+1}.name = 'Nate2017';
% storms{end}.days = '09/16/2017-10/08/2017';

% storms{end+1}.name = 'Ian2022';
% storms{end}.days = '09/19/2022-10/02/2022';


fileID = fopen('download_wind_pressure_cfsv2_tides.sh','w');
fprintf(fileID,'%s \n','start=`date +%s`');
fprintf(fileID,'%s \n','shopt -s nullglob');

for i = 1:length(storms)
    
    sep = strfind(storms{i}.days,'-');
    day1 = datetime(storms{i}.days(1:sep-1));
    day2 = datetime(storms{i}.days(sep+1:end));
    
    mkdir_str = ['mkdir ',storms{i}.name];
    fprintf(fileID,'%s \n',mkdir_str);
    cddir_str = ['cd ',storms{i}.name];
    fprintf(fileID,'%s \n',cddir_str);
    
%     fileID2 = fopen([storms{i}.name,'.txt'],'w');
%     fprintf(fileID2,'%s \n',storms{i}.name);
%     fprintf(fileID2,'%s \n',datetime(day1-7,'format','yyyy-MM-dd HH:mm:ss'));
%     fprintf(fileID2,'%s \n',datetime(day2,'format','yyyy-MM-dd HH:mm:ss'));
%     fclose(fileID2);
    
    if day1<datetime('04/01/2011') % Use Reanalysis
        
        %link1 = 'https://www.ncei.noaa.gov/data/climate-forecast-system/access/reanalysis/time-series';        
        link1 = 'https://www.ncei.noaa.gov/oa/prod-cfs-reanalysis/time-series';
        
        yr = year(day1);
        month1 = month(day1);
        month2 = month(day2);
        months = unique([month1,month2]);
        if length(months)==1
            days = [num2str(yr),num2str(months(1),'%02.0f')];
            path1 = [link1,'/',days,'/'];
            filename = ['wnd10m.gdas.',days,'.grb2'];
            fprintf(fileID,'%s \n',['wget ',path1,filename]);
            fprintf(fileID,'%s \n',['mv ',filename,' fort.222.grb2']);
            
            path1 = [link1,'/',days,'/'];
            filename = ['prmsl.gdas.',days,'.grb2'];
            fprintf(fileID,'%s \n',['wget ',path1,filename]);
            fprintf(fileID,'%s \n',['mv ',filename,' fort.221.grb2']);
        else
            for j = 1:length(months)
                if j ==1
                    
                    days = [num2str(yr),num2str(months(j),'%02.0f')];
                    path1 = [link1,'/',days,'/'];
                    filename = ['wnd10m.gdas.',days,'.grb2'];
                    fprintf(fileID,'%s \n',['wget ',path1,filename]);
                    fprintf(fileID,'%s \n',['mv ',filename,' fort.222.grb2']);
                    
                    path1 = [link1,'/',days,'/'];
                    filename = ['prmsl.gdas.',days,'.grb2'];
                    fprintf(fileID,'%s \n',['wget ',path1,filename]);
                    fprintf(fileID,'%s \n',['mv ',filename,' fort.221.grb2']);
                    
                else
                    days = [num2str(yr),num2str(months(j),'%02.0f')];
                    path1 = [link1,'/',days,'/'];
                    filename = ['wnd10m.gdas.',days,'.grb2'];
                    fprintf(fileID,'%s \n',['wget ',path1,filename]);
                    fprintf(fileID,'%s \n',['cat fort.222.grb2 ',filename,'  > fort.222.grb2']);
                    
                    path1 = [link1,'/',days,'/'];
                    filename = ['prmsl.gdas.',days,'.grb2'];
                    fprintf(fileID,'%s \n',['wget ',path1,filename]);
                    fprintf(fileID,'%s \n',['cat fort.221.grb2 ',filename,'  > fort.221.grb2']);
                    
                end
            end
        end
    else % Use CFSv2 Analysis
        link1 = 'https://www.ncei.noaa.gov/data/climate-forecast-system/access/operational-analysis/time-series';
        yr = year(day1);
        month1 = month(day1);
        month2 = month(day2);
        months = unique([month1,month2]);
        if length(months) ==1
            days = [num2str(yr),num2str(months(1),'%02.0f')];
            path1 = [link1,'/',num2str(yr),'/',days,'/'];
            filename = ['wnd10m.gdas.',days,'.grib2'];
            fprintf(fileID,'%s \n',['wget',path1,filename]);
            %fprintf(fileID,'%s \n',['mv ',filename,' fort.222.grb2']);
            
            path1 = [link1,'/',num2str(yr),'/',days,'/'];
            filename = ['prmsl.gdas.',days,'.grib2'];
            fprintf(fileID,'%s \n',['wget ',path1,filename]);
            %fprintf(fileID,'%s \n',['mv ',filename,' fort.221.grb2']);
            
        else
            for j = 1:length(months)
                if j==1
                    days = [num2str(yr),num2str(months(j),'%02.0f')];
                    path1 = [link1,'/',num2str(yr),'/',days,'/'];
                    filename = ['wnd10m.gdas.',days,'.grib2'];
                    fprintf(fileID,'%s \n',['wget ',path1,filename]);
                    %fprintf(fileID,'%s \n',['mv ',filename,' fort.222.grb2']);
                    
                    path1 = [link1,'/',num2str(yr),'/',days,'/'];
                    filename = ['prmsl.gdas.',days,'.grib2'];
                    fprintf(fileID,'%s \n',['wget ',path1,filename]);
                    %fprintf(fileID,'%s \n',['mv ',filename,' fort.221.grb2']);
                else
                    days = [num2str(yr),num2str(months(j),'%02.0f')];
                    path1 = [link1,'/',num2str(yr),'/',days,'/'];
                    filename = ['wnd10m.gdas.',days,'.grib2'];
                    fprintf(fileID,'%s \n',['wget ',path1,filename]);
                    fprintf(fileID,'%s \n',['cat fort.222.grb2 ',filename,' > fort.222.grb2']);
                    %fprintf(fileID,'%s \n',['rm ',filename]);
                    
                    path1 = [link1,'/',num2str(yr),'/',days,'/'];
                    filename = ['prmsl.gdas.',days,'.grib2'];
                    fprintf(fileID,'%s \n',['wget ',path1,filename]);
                    fprintf(fileID,'%s \n',['cat fort.221.grb2 ',filename,' > fort.221.grb2']);
                    %fprintf(fileID,'%s \n',['rm ',filename]);
                end
            end
        end
    end
    
%     fprintf(fileID,'%s \n',['mv ./../',storms{i}.name,'.txt ./']);
%     fprintf(fileID,'%s \n','cd ..');
    
    
end


fprintf(fileID,'%s \n','end=`date +%s`');
fprintf(fileID,'%s \n','runtime=$((end-start))');
fprintf(fileID,'%s \n','echo $runtime');
fclose(fileID);


system('bash download_wind_pressure_cfsv2_tides.sh')



