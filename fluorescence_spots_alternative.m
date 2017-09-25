%script to analyze fluorescence pictures of droplet arrays (time-lapse)
clc; %clear screen
close all; %close figures
%% definitions
%array positions are from left to right and up to down in increasing number
%10x magnification was used with 250 um spot size and 500 um spot distance
%(borders have to be changed otherwise)
%distance between picture rows and columns is always the same
%Youscope recorded time and filenames

%% load and extract data
directory=uigetdir('/Users/nnuti/Documents/Nicola (local)/Projects');
files = dir([directory,'\*.jpg']);
positions=zeros(size(files,1),2);
time=zeros(size(files,1),1);
image_properties = readtable([directory,'\images.csv'],'Delimiter','semi');

%extract images and properties
data_pic=cell(size(files,1),4); %position(picture),time,image, results image
for j =1:size(files,1)
    data_pic(j,1)={[str2double(image_properties.ImageFile{j}(10:11)),str2double(image_properties.ImageFile{j}(12:13))]}; %extract position of image
    data_pic(j,2)={str2double(image_properties.MeasurementTime_ms_{j})}; %extract time
    if j==1
        starttime= data_pic{1,2};
    end;
    data_pic(j,2) = {data_pic{j,2}-starttime};          %reference time to the beginning
    data_pic(j,3)= {imread([directory,'\',image_properties.ImageFile{j}])};     %extract images
end;

%% image processing
borders=[50,200,400;200,400,0]; %borders for spots in image
for j=1:size(files,1)
    %read image
    img=data_pic{j,3};
    
    %mask generation
    img_bw=im2bw(img,graythresh(img));      %imbinarize
    se = strel('disk',8);
    img_bw_2=imclose(img_bw,se);
    mask=imerode(img_bw_2,se); %mask for spot detection
    
    %detect droplets
    %img_grayscale = rgb2gray(img);
    regionprops = regionprops(mask,img,'EquivDiameter','WeightedCentroid','PixelValues');
    centers = zeros(size(regionprops,1),2);
    fluorescence_median= zeros(size(regionprops,1),1);
    fluorescence_mean= zeros(size(regionprops,1),1);
    radii = zeros(size(regionprops,1),1);
    for i = 1:size(regionprops,1)
        centers(i,:)=regionprops(i).WeightedCentroid;  %access centers
        fluorescence_median(i) = median(regionprops(i).PixelValues);
        fluorescence_mean(i) = mean(regionprops(i).PixelValues);
        radii(i,:)=regionprops(i).EquivDiameter/2;
    end;
    clear regionprops;
   %identify row and column of droplet in picture, centers (x,y,col,row),
   %col=3 means out of boundaries
    data_position = zeros(size(centers,2),2);
    for i=1:size(centers,1)
        if centers(i,1)>borders(1,1) && centers(i,1)<borders(1,2)
            data_position(i,1)=1;
        elseif centers(i,1)>borders(1,2) && centers(i,1)<borders(1,3)
            data_position(i,1)=2;
        else
            data_position(i,1)=3;
        end;
        if centers(i,2)<borders(2,1) 
            data_position(i,2)=1;
        elseif centers(i,2)>borders(2,1) && centers(i,2)<borders(2,2)
            data_position(i,2)=2;
        else
            data_position(i,2)=3;
        end;
    end;
    data = [centers,data_position(:,[2,1]),fluorescence_median,fluorescence_mean,radii];
    data=sortrows(data,[3,4]);
    data_pic(j,4)= {data};  %results image: x,y,row,col,median int.,mean int.
end;

%% quality control of image analysis
for j=1:size(files,1)
    imshow(data_pic{j,3});
    title(num2str(j));
    viscircles(data_pic{j,4}(:,1:2), data_pic{j,4}(:,7),'LineStyle','--');
    pause(0.2);
end;

%% Combine pictures, data sorting
data_spot=zeros(size(files,1)*6,4); % row(total),column(total),time,median fluorescence,  
k=1; %Count Datenpunkte     %attention row and column detection doesnt work
for j=1:size(files,1)
    data_spot(k,1)= (data_pic{j,1}(1)-1)*3+data_pic{j,4}(1,3); %row(total)
    data_spot(k,2)= (data_pic{j,1}(2)-1)*2+data_pic{j,4}(1,4); %row(total)
    data_spot(k,3)= data_pic{j,2}; %time
    data_spot(k,4)= data_pic{j,4}(1,5);
    k=k+1;
    for i=2: size(data_pic{j,4},1) % number of spot in picture
        if data_pic{j,4}(i,4)<3
            if or(ne(data_pic{j,4}(i,4),data_pic{j,4}(i-1,4)), ne(data_pic{j,4}(i,3),data_pic{j,4}(i-1,3)))
                data_spot(k,1)= (data_pic{j,1}(1)-1)*3+data_pic{j,4}(i,3); %row(total)
                data_spot(k,2)= (data_pic{j,1}(2)-1)*2+data_pic{j,4}(i,4); %row(total)
                data_spot(k,3)= data_pic{j,2};
                data_spot(k,4)= data_pic{j,4}(i,5);
                k=k+1;
            end;
        end;
        
    end;
end;
display([num2str(size(data_spot,1)-k+1),' spots were not detected']);
data_spot = data_spot(1:k-1,:); %delete empty rows

%% Kinetics
%sort kinetics
data_spot_kinetics=cell(max(data_spot(:,1)),max(data_spot(:,2)));
for k=1:size(data_spot,1)
    data_spot_kinetics{data_spot(k,1),data_spot(k,2)} = [data_spot_kinetics{data_spot(k,1),data_spot(k,2)};[data_spot(k,3),data_spot(k,4)]];
end;
% fit and calculate activity values
data_activity=cell(size(data_spot_kinetics));
for i=1:size(data_spot_kinetics,1) %row
  for j=1:size(data_spot_kinetics,2)  %column
      if size(data_spot_kinetics{i,j},1)>= 3
        [fitdata,gof]=fit(data_spot_kinetics{i,j}(1:17,1)/1000,data_spot_kinetics{i,j}(1:17,2),'poly1');
        data_activity{i,j}(1,:) = coeffvalues(fitdata);% activity
        data_activity{i,j}(2,:) = [gof.rsquare,gof.rmse];% fit exactness
        %plot kinetics
        figure;
        plot(fitdata,data_spot_kinetics{i,j}(1:17,1)/1000,data_spot_kinetics{i,j}(1:17,2)); %change timepoints for new data
        title(['row: ',num2str(i),', column: ',num2str(j)]);
        pause(0.2);
      end;
  end;
end;

% plot activities
for i=1:size(data_spot_kinetics,1) %row
  for j=1:size(data_spot_kinetics,2)  %column
      if size(data_activity{i,j},1)>= 0
        plot([i,j],data_activity{i,j}(1,1));
      end;
  end;
end;