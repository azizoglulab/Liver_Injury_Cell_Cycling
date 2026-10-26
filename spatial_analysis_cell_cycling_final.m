%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%             Spatial Proliferation Analysis    %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Sherry Li Zheng (slzheng@stanford.edu)      Updated: 10-01-2026 %
% Updates by: Shannon Rhoads (srhoads@unc.edu)        Updated: 10-16-2026 %


%% INSTRUCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.Update the parameters specific to your data in the "USER INPUT"
%   section below.
% 2.Run script and provide user defined regions for the central vein (CV)
%   centroid and boarder when instructed. To define the CV regions:
    % - Wait for a pop up window to appear displaying your first fluorescently
    %   labeled lobule image. 
    % - In the first tab that appears, DOUBLE CLICK on the center of the
    %   central vein.
    % - In the second tab that appears, CLICK AND HOLD to draw a region of
    %   interest around the boarder of the central vein. 
    % - Continue this process for all images in the dataset. Each image
    %   will be display consecutively.
% 3.After manual inputs, quantification will run automatically and be
%   saved in the specified output folders. Read below for more information
%   on the expected inputs and outputs.




%% OVERVIEW OF CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: 
    % This script is used to analyze the spatial distribution of fluorescently
    % marked nuclei of interest in a liver lobule.

% INPUTS: 
    % 1.Single channel 2D images of fluorescently-labeled liver lobules 
    %   (called "Procr_signal" below). Expecting: ImageN_X.tiff
    % 2.Single channel 2D instance segmentation or fluorescence intensity 
    %   image of marked nuclei of interest. Example: ERG/Ki67 double
    %   positive nuclei. Expecting: ImageN_X.tif
    
    % Naming convensions:
        % - N = image number
        % - X = lobule marker name 

% PIPELINE COMPONENTS: 
    % USER INPUT
        % 1. Sample number
        % 2. Number of images
        % 3. Brightness of liver lobule signal
        % 4. Paths to the input and output folders
    
    % SECTION 1
        % 1. Create structural elements/initialize data cell arrays for
        %    downstream processing
        % 2. Collect manual inputs for Procr line mask (CV centroid and boundary)
        % 3. Create mask for Procr line
    
    % SECTION 2
        % 4. Create (or import) masks for marked nuclei of interest
        % 5. Analysis of marked nuclei of interest
            % a. Exclude nuclei outside Procr line
            % b. Draw CV-PV and CV-centroid lines
            % c. Draw CV-clone edge and PV-clone edge lines
            % d. Run calculations and fill in single cell data cell array
            
    % SECTION 3   
        % 6. Export data (including Procr calculations) in csv files

% OUTPUT:
    % 1.One .csv file containing the percent distance of each marked nuclei
    %   of interest to the central vein (CV) and peripheral
    %   vein (PV) and the area in pixels of the clone. This file is saved to
    %   the "quant_outpath" location specified below.
    % 2.Four .jpg images displaying the line used to calculate the
    %   distance of the clone from CV to PV. These files are saved to
    %   the "quant_outpath" location specified below.
    % 3.One .bmp image displaying the lobules intensity image with the
    %   central vein area removed. This file is saved to the "quant_outpath"
    %   location specified below. 
    % 4.One .jpg image displaying all of the clones labeled (outline
    %   color) by their classification: inside or outside of the Procr
    %   line. This file is saved to the "mask_outpath" location specified
    %   below. 
    % 5.One .jpg image displaying only the marked nuclei of interest.
    %   This file is saved to the "mask_outpath" location specified
    %   below.



%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EDIT THIS SECTION BEFORE RUNNING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set sample and number of images included in analysis
Exp_sample_ID = 0;      % adjust to reflect correct sample number (ExpXX); for naming purposes only
first_image = 1;        % Image number (N) of the first image to process
last_image = 3;         % Image number (N) of the LAST image to process -- all images between first_image and last_image will be included in the analysis

% Indicate how bright your 
Procr_signal = 2;       % enter 1 if Procr_signal is bright; enter 2 if Procr_signal is faint

% INPUT file paths:
Procr_inpath = 'C:/Users/Folder/EPCR';        % Path to the fluorescently labeled lobules images
marked_nuclei_mask_inpath = 'C:/Users/Folder/Overlap mask';     % Path to the mask images of marked nuclei of interest

% OUTPUT file paths: 
quant_outpath = 'C:/Users/Folder/EPCR';       % Path to save the quantitative data and related central to peripheral vein .jpg
mask_outpath = 'C:/Users/Folder/Overlap mask';        % Path to save the masks of marked nuclei of interest only
    


%% SECTION 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic %starts timer

% Create structural elements
disk_4 = strel('disk',4); % small disk
disk_10 = strel('disk',10); % big disk

% Preallocate CV coordinate matrix
CV_point = zeros(last_image,2);

% Constants for filling data arrays; do not adjust 
clone_sample_number = 0; 
singlecell_sample_number = 0; 

% Initialize data cell arrays
clone_data_cell_array = {}; 
singlecell_data_cell_array = {}; 
Procr_data_cell_array = {}; 

% If output file locations don't exist already, create them
if ~isfolder(quant_outpath)
    mkdir(quant_outpath);
    fprintf('Created folder: %s\n', quant_outpath);
end

if ~isfolder(mask_outpath)
    mkdir(mask_outpath);
    fprintf('Created folder: %s\n', mask_outpath);
end


%% Manual inputs for Procr line mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for image_number = first_image:last_image 
    
    % Load Procr image
    unaltered_Procr_image=imread(strcat(Procr_inpath,'/Image',num2str(image_number),'_Procr.tif'));
    
    % Collect CV center point and mask
    % display instructions for first image only
    if image_number == first_image
        % Select CV point by double clicking
        figure; imshow(unaltered_Procr_image);
        uiwait(msgbox('Double-click: Select the center of the central vein.','Instructions','modal'));
        [cv_x, cv_y] = getpts; %double click

        % Create a binary mask around the central vein area by clicking and dragging in a circle; release click to finish circle
        figure; imshow(unaltered_Procr_image, []);
        uiwait(msgbox('Freehand ROI: Click and drag to draw a mask around the central vein. Release the click when finished.','Instructions','modal'));
        subtract_CV = imfreehand(); %click and drag in a circle, then release

        if image_number ~= last_image
            uiwait(msgbox('Repeat the double-click, then the freehand ROI for all remaining images in the dataset. The images will continue to pop up until they have all been displayed.','Instructions','modal'));
        end
    
    else
        % Select CV point by double clicking
        figure; imshow(unaltered_Procr_image);
        [cv_x, cv_y] = getpts; %double click

        % Create a binary mask around the central vein area by clicking and dragging in a circle; release click to finish circle
        figure; imshow(unaltered_Procr_image, []);
        subtract_CV = imfreehand(); %click and drag in a circle, then release
    end

    % Add CV centroid to image
    CV_point(image_number,1) = cv_x;
    CV_point(image_number,2) = cv_y;

    % Create mask of area outside CV based on selection
    deleted_CV_area = subtract_CV.createMask();

    % Remove CV area from lobule intensity image
    adjusted_Procr = unaltered_Procr_image;
    adjusted_Procr(deleted_CV_area) = 0;
    
    % Save as new Procr image 
    file_location = fullfile(quant_outpath,sprintf(strcat('Image',num2str(image_number),'_Procr_Exp',num2str(Exp_sample_ID),'.bmp')));
    imwrite(adjusted_Procr,file_location); %use this image for next section

    if image_number == last_image
        uiwait(msgbox('All images have been displayed. Now proceeding with the quantification.','Instructions','modal'));
    end

end
     

%% Create Procr line mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for image_number = first_image:last_image 
    % Load new Procr images
    original_Procr_image = imread(strcat(Procr_inpath,'/Image', num2str(image_number), '_Procr.tif')); %for later overlays
    Procr_image = imread(strcat(quant_outpath,'/Image',num2str(image_number),'_Procr_Exp',num2str(Exp_sample_ID),'.bmp'));
    Procr_img_size = size(Procr_image);
    
    if Procr_signal == 1 %if Procr_signal is bright
        
        % Adjust contrast and sort intensity values
        Procr_contrast = adapthisteq(Procr_image);
        Procr_contrast_sort = sort(Procr_contrast(:));
    
        % Use background subtraction (100 AU) and intensity thresholding (select top 15%)
        Procr_subtract = Procr_contrast-100; 
        Procr_premask = Procr_subtract>(max(0.15*Procr_contrast_sort));
        
        % Filter to remove small debris (10 pxl or less); blur removes more
        Procr_filter = bwareaopen(Procr_premask,10);
        Procr_blur = imfilter(Procr_filter,fspecial('disk',3));
    end
    
    if Procr_signal == 2 %if Procr_signal is faint
        
        % Skip contrast adjustment and background subtraction; filter and blur
        Procr_filter = bwareaopen(Procr_image,10);
        Procr_blur = imfilter(Procr_filter,fspecial('disk',3));
    end
    
    % Use dilation to connect the Procr mask fragments; some images may need extra dilation if stain is more faint or patchy
    Procr_dilate1 = imdilate(Procr_blur,disk_10);
    Procr_dilate2 = imdilate(Procr_dilate1,disk_10);
    Procr_dilate3 = imdilate(Procr_dilate2,disk_10);
    Procr_enclose=imclose(Procr_dilate3,true(50)); 
    Procr_linemask=bwskel(Procr_enclose);
    pre_fatline = imdilate(Procr_linemask,disk_10);
    fatline = bwareaopen(pre_fatline,2000);
    
    % Determine whether mask needs to be dilated further by testing whether line is continuous (whether inside and outside areas can be connected in any way)
    flipped_fatline = ~fatline; %fatline = 0, other areas = 1
    bw_flipped_fatline = bwlabel(flipped_fatline); %counts number of other areas
    max_bw_flipped_fatline = max(bw_flipped_fatline(:)); %if >1, fatline is continuous; if =1, fatline is broken and needs to be dilated more
    
    % Additional code for further dilation, if needed
    if max_bw_flipped_fatline == 1
        Procr_dilate4 = imdilate(Procr_dilate3,disk_10);
    	Procr_dilate5 = imdilate(Procr_dilate4,disk_10);
    	Procr_dilate6 = imdilate(Procr_dilate5,disk_10);
        Procr_dilate7 = imdilate(Procr_dilate6,disk_10);
        Procr_dilate8 = imdilate(Procr_dilate7,disk_10);
    	Procr_enclose=imclose(Procr_dilate8,true(50)); 
    	Procr_linemask=bwskel(Procr_enclose);
    	fatline = imdilate(Procr_linemask,disk_10);
        flipped_fatline = ~fatline;
        bw_flipped_fatline = bwlabel(flipped_fatline);
        max_bw_flipped_fatline = max(bw_flipped_fatline(:));
        if max_bw_flipped_fatline == 1 %more additional code, if needed
            Procr_dilate9 = imdilate(Procr_dilate8,disk_10);
            Procr_dilate10 = imdilate(Procr_dilate9,disk_10);
            Procr_enclose=imclose(Procr_dilate10,true(50)); 
            Procr_linemask=bwskel(Procr_enclose);
            fatline = imdilate(Procr_linemask,disk_10);
            flipped_fatline = ~fatline;
            bw_flipped_fatline = bwlabel(flipped_fatline);
            max_bw_flipped_fatline = max(bw_flipped_fatline(:));
            if max_bw_flipped_fatline == 1 %more additional code, if needed
                Procr_dilate11 = imdilate(Procr_dilate10,disk_10);
                Procr_dilate12 = imdilate(Procr_dilate11,disk_10);
                Procr_enclose=imclose(Procr_dilate12,true(50)); 
                Procr_linemask=bwskel(Procr_enclose);
                fatline = imdilate(Procr_linemask,disk_10);
            end
        end
    end

  
    % Export Procr line as binary to check for continuity 
    file_location = fullfile(quant_outpath,sprintf(strcat(num2str(image_number),'_Procr_line_Exp',num2str(Exp_sample_ID),'.jpg')));
    imwrite(fatline,file_location);

    % Export Procr line as an overlay on the CV removed image
    [Procr_image_height,Procr_image_width]=size(Procr_image);
    Procr_check=zeros(Procr_image_height,Procr_image_width,3);
    Procr_check(:,:,1)=bwperim(Procr_linemask);
    Procr_check(:,:,3)=imadjust(mat2gray(Procr_image));
    file_location = fullfile(quant_outpath,sprintf(strcat(num2str(image_number),'_Procr_overlay1_Exp',num2str(Exp_sample_ID),'.jpg')));
    imwrite(Procr_check,file_location);

    % Export Procr line as an overlay on the original image to check accuracy
    full_Procr_check=zeros(Procr_image_height,Procr_image_width,3);
    full_Procr_check(:,:,2)=bwperim(Procr_linemask);
    full_Procr_check(:,:,3)=imadjust(mat2gray(original_Procr_image));
    file_location = fullfile(quant_outpath,sprintf(strcat(num2str(image_number),'_Procr_overlay2_Exp',num2str(Exp_sample_ID),'.jpg')));
    imwrite(full_Procr_check,file_location);

    
    % Get pixel values for Procr linemask border (used for minimum distance measurements later on)
    Procr_region_props = regionprops(Procr_linemask,'PixelList','Area','Perimeter');
    Procr_pixel_list = Procr_region_props.PixelList;
    Procr_pixel_list_x = round(Procr_pixel_list(:,1));
    Procr_pixel_list_y = round(Procr_pixel_list(:,2));
    

    
    %% SECTION 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Analysis of marked nuclei of interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Load images of marked nuclei of interest; sort contrast
    clone_image = imread(strcat(marked_nuclei_mask_inpath,'/Image',num2str(image_number),'_X.tiff'));
    clone_contrast = imadjust(clone_image);
    clone_contrast_sort = sort(clone_contrast(:));
    
    % Use background subtraction and intensity thresholding
    clone_subtract = clone_contrast; 
    clone_mask = clone_subtract>(max(0.01*clone_contrast_sort));
    
    % Filter out single cells
    clone_mask = bwareaopen(clone_mask,0); 
    original_clone_mask = clone_mask; %for later overlays
    
    % Export clone masks overlayed on original clone images to check for fit
    [clone_image_height,clone_image_width]=size(clone_image);
    clone_check=zeros(clone_image_height,clone_image_width,3);
    clone_check(:,:,1)=imadjust(mat2gray(clone_image));
    clone_check(:,:,2)=bwperim(clone_mask);
    file_location = fullfile(mask_outpath,sprintf(strcat(num2str(image_number),'_clone_masks_Exp',num2str(Exp_sample_ID),'_X.jpg')));
    imwrite(clone_check,file_location); %check to make sure masks are precise and exclude single cells
    
    
    %% Exclude clones outside Procr line %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Find clone centroid values and label clone masks
    CC = bwconncomp(clone_mask, 8);
    clone_region_props1 = regionprops(CC,'centroid','ConvexHull'); 
    clone_centroid_matrix1 = cat(1,clone_region_props1.Centroid);
    clone_length1 = length(clone_centroid_matrix1);
    
    % If there are no clones in the lobule, fill data tab with nothing; analysis ends here
    if clone_length1 == 0
        clone_sample_number = clone_sample_number + 1;
       	identifier = strcat('image_',num2str(image_number),'_X');
       	clone_out_row = {identifier, 'NA', 'NA', 'NA'};
      	clone_data_cell_array(clone_sample_number,:) = clone_out_row;
        
    % Otherwise, move on    
    else
    
        % Fill Procr mask and test whether each centroid is inside filled mask
        % If yes, keep centroid; if no, remove the clone before downstream quantification
        % If centroid is zero, multiply clone mask to (zeroed) convexhull (masking the clone mask) to zero out that clone
        ProcrMask = imfill(Procr_linemask,'holes');

        % CC = bwconncomp(clone_mask, 8); % count number of clones
        % props = regionprops(CC, 'Centroid');
        keep = false(1,CC.NumObjects); % create array to flag objects to keep
        
        % loop through array and determine which clones are inside Procr mask
        for k = 1:CC.NumObjects
            c = clone_region_props1(k).Centroid;
            cx = round(c(1)); cy = round(c(2));
            % skip out-of-bounds centroids
            if cx < 1 || cx > size(ProcrMask,2) || cy < 1 || cy > size(ProcrMask,1)
                keep(k) = false;
            else
                keep(k) = ProcrMask(cy, cx); % true if centroid inside Procr
            end
        end
        
        % build cleaned mask: keep only components with centroid inside ProcrMask
        L = labelmatrix(CC);
        clean_clone_mask = ismember(L, find(keep));

        % Remove any debris left over (due to rounding, etc) and check overlay
        clone_mask_final = bwareaopen(clean_clone_mask,10);
        clone_mask_check=zeros(clone_image_height,clone_image_width,3);
        clone_mask_check(:,:,1)=bwperim(clone_mask);
        clone_mask_check(:,:,2)=bwperim(clone_mask_final);
        clone_mask_check(:,:,3)=bwperim(Procr_linemask);
        file_location = fullfile(mask_outpath,sprintf(strcat(num2str(image_number),'_clone_exclusion_Exp',num2str(Exp_sample_ID),'_X.jpg')));
        imwrite(clone_mask_check,file_location); %red are excluded, yellow are included

        % Find centroid and area measurements for all included clones
        clone_region_props2 = regionprops(clone_mask_final,'area','centroid','PixelList','ConvexHull'); 
        clone_centroid_matrix2 = cat(1,clone_region_props2.Centroid);
        clone_length2 = length(clone_centroid_matrix2);
        clone_area_matrix = cat(1,clone_region_props2.Area); 
        filled_clone_mask = imfill(clone_mask_final, 'holes');

        % If there are no clones in the lobule, analysis ends here
        if clone_length2 == 0
            clone_sample_number = clone_sample_number + 1;
            identifier = strcat('image_',num2str(image_number),'_X');
            clone_out_row = {identifier, 'NA', 'NA', 'NA'};
            clone_data_cell_array(clone_sample_number,:) = clone_out_row;

        % Otherwise, move on    
        else
            %% Draw lines from CV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for clone_number = 1:length(clone_centroid_matrix2(:,1))
                %% Draw line from CV to image edge %%%%%%%%%%%%%%%%%%%%%%%%

                % Find equation of line from CV through centroid to image edge
                clone_centroid_x=clone_centroid_matrix2(clone_number,1);
                clone_centroid_y=clone_centroid_matrix2(clone_number,2); 
                i_x=[(CV_point(image_number,1)) clone_centroid_x];
                i_y=[(CV_point(image_number,2)) clone_centroid_y];
                i_intercept_slope = [[1; 1]  i_x(:)]\i_y(:);
                i_intercept = i_intercept_slope(1);
                i_slope = i_intercept_slope(2);

                    % If clone is left of the CV, extend line to the left edge
                    % If clone is right of the CV, extend line to the right edge
                    if clone_centroid_x < (CV_point(image_number,1)) 
                        end_x = 0;
                        end_y = i_intercept; 
                    else 
                        end_x = Procr_img_size(2);
                        end_y = end_x*i_slope+i_intercept;
                    end

                % Draw line from CV to image edge
                i_nPoints = ceil(sqrt((end_x - (CV_point(image_number,1))).^2));
                i_xvalues = round(linspace((CV_point(image_number,1)), end_x, i_nPoints));
                i_yvalues = round(linspace((CV_point(image_number,2)), end_y, i_nPoints));
                CV_edge_linemask = zeros(Procr_img_size);
                    for i_linepoint = 1:length(i_xvalues)
                        if i_xvalues(i_linepoint) > 0 && i_xvalues(i_linepoint) < Procr_img_size(2)
                           if i_yvalues(i_linepoint) > 0 && i_yvalues(i_linepoint) < Procr_img_size(1)
                             CV_edge_linemask(sub2ind(Procr_img_size, i_yvalues(i_linepoint), i_xvalues(i_linepoint))) = 1;
                           end
                        end
                    end


                %% Draw line from CV to Procr edge %%%%%%%%%%%%%%%%%%%%%%%%

                % Find intersection between CV-image edge and Procr mask lines
                dilated_linemask = imdilate(CV_edge_linemask, disk_10); 
                multiplied_mask = dilated_linemask .* Procr_linemask;
                mult_check = max(max(multiplied_mask)); %if =0, there is no intersection because the (dotted) CV_edge_linemask has missed the Procr_linemask

                % Since line is dotted, CV_edge_linemask may require (a lot of) extra dilation
                if mult_check == 0
                    dilated_linemask = imdilate(dilated_linemask, disk_4); 
                    multiplied_mask = dilated_linemask .* Procr_linemask;
                    mult_check = max(max(multiplied_mask));
                    if mult_check == 0 %further dilation, if needed
                        dilated_linemask = imdilate(dilated_linemask, disk_4); 
                        multiplied_mask = dilated_linemask .* Procr_linemask;
                        mult_check = max(max(multiplied_mask));
                        if mult_check == 0 %further dilation, if needed
                            dilated_linemask = imdilate(dilated_linemask, disk_10); 
                            multiplied_mask = dilated_linemask .* Procr_linemask;
                            mult_check = max(max(multiplied_mask));
                        end
                    end
                end
                
                % If there is still no intercept, analysis ends here with an error
                if mult_check == 0
                    clone_sample_number = clone_sample_number + 1;
                    identifier = strcat('image_',num2str(image_number),'_','clone_',num2str(clone_number),'_X');
                    clone_out_row = {identifier, 'ER', 'ER', 'ER'};
                    clone_data_cell_array(clone_sample_number,:) = clone_out_row;

                % Otherwise, move on    
                else

                    % Get intersection points (where multiplied mask produces a non-zero number)
                    [cvProcr_y,cvProcr_x] = find (multiplied_mask > 0);

                    % Draw line from CV to Procr edge using the median intersection point
                    j_nPoints = ceil(sqrt(((median(cvProcr_x)) - (CV_point(image_number,1))).^2 + ((median(cvProcr_y)) - (CV_point(image_number,2))).^2)); 
                    j_xvalues = round(linspace((CV_point(image_number,1)), (median(cvProcr_x)), j_nPoints)); 
                    j_yvalues = round(linspace((CV_point(image_number,2)), (median(cvProcr_y)), j_nPoints));
                    CV_Procr_linemask = zeros(Procr_img_size); 
                    CV_Procr_linemask(sub2ind(Procr_img_size, j_yvalues, j_xvalues)) = 1; 
                    CV_Procr_linemask = imdilate(CV_Procr_linemask, disk_4);

                    %% Draw line from CV to clone centroid %%%%%%%%%%%%%%%%

                    % Truncate CV-Procr line at the centroid
                    k_nPoints = ceil(sqrt((clone_centroid_x - (CV_point(image_number,1))).^2 + (clone_centroid_y - (CV_point(image_number,2))).^2)); 
                    k_xvalues = round(linspace((CV_point(image_number,1)), clone_centroid_x, k_nPoints));
                    k_yvalues = round(linspace((CV_point(image_number,2)), clone_centroid_y, k_nPoints));
                    CV_centroid_linemask = zeros(Procr_img_size);
                    CV_centroid_linemask(sub2ind(Procr_img_size, k_yvalues, k_xvalues)) = 1;
                    CV_centroid_linemask = imdilate(CV_centroid_linemask, disk_4);

                    %% Draw line from CV and PV to closest clone edge %%%%%

                    % Obtain list of pixel coordinates at the perimeter of the clone mask
                    clone_pixel_list = clone_region_props2(clone_number).PixelList;
                    clone_pixel_list_x = round(clone_pixel_list(:,1));
                    clone_pixel_list_y = round(clone_pixel_list(:,2));

                    % Loop through to find shortest distance between clone mask edge and CV
                    CV_clone_minimum_dist = Inf;
                    PV_clone_minimum_dist = Inf;
                    for clone_pixel_number = 1:length(clone_pixel_list_x)
                        CV_coordinate_pairs = [clone_pixel_list_x(clone_pixel_number) clone_pixel_list_y(clone_pixel_number); (CV_point(image_number,1)) (CV_point(image_number,2))];
                        CV_dist = pdist(CV_coordinate_pairs);
                        if CV_dist < CV_clone_minimum_dist
                            CV_clone_minimum_dist = CV_dist;
                            clone_CV_edge_x = clone_pixel_list_x(clone_pixel_number);
                            clone_CV_edge_y = clone_pixel_list_y(clone_pixel_number);
                        end

                        % Loop through to find shortest distance between clone mask edge and PV
                        for Procr_pixel_number = 1:length(Procr_pixel_list(:,1))
                            PV_coordinate_pairs = [clone_pixel_list_x(clone_pixel_number) clone_pixel_list_y(clone_pixel_number); Procr_pixel_list_x(Procr_pixel_number) Procr_pixel_list_y(Procr_pixel_number)];
                            PV_dist = pdist(PV_coordinate_pairs);
                            if PV_dist < PV_clone_minimum_dist
                                PV_clone_minimum_dist = PV_dist;
                                clone_PV_edge_x = clone_pixel_list_x(clone_pixel_number);
                                clone_PV_edge_y = clone_pixel_list_y(clone_pixel_number);
                                min_Procr_x = Procr_pixel_list_x(Procr_pixel_number);
                                min_Procr_y = Procr_pixel_list_y(Procr_pixel_number);
                            end
                        end
                    end


                    % Draw shortest line between CV and clone edge
                    l_nPoints = ceil(sqrt((clone_CV_edge_x - (CV_point(image_number,1))).^2 + (clone_CV_edge_y - (CV_point(image_number,2))).^2)); 
                    l_xvalues = round(linspace((CV_point(image_number,1)), clone_CV_edge_x, l_nPoints));
                    l_yvalues = round(linspace((CV_point(image_number,2)), clone_CV_edge_y, l_nPoints));
                    CV_clone_linemask = zeros(Procr_img_size);
                    CV_clone_linemask(sub2ind(Procr_img_size, l_yvalues, l_xvalues)) = 1;
                    CV_clone_linemask = imdilate(CV_clone_linemask, disk_4);

                    % Draw shortest line between PV and clone edge
                    m_nPoints = ceil(sqrt((clone_PV_edge_x - min_Procr_x).^2 + (clone_CV_edge_y - min_Procr_y).^2)); 
                    m_xvalues = round(linspace(min_Procr_x, clone_PV_edge_x, m_nPoints));
                    m_yvalues = round(linspace(min_Procr_y, clone_PV_edge_y, m_nPoints));
                    PV_clone_linemask = zeros(Procr_img_size);
                    PV_clone_linemask(sub2ind(Procr_img_size, m_yvalues, m_xvalues)) = 1;
                    PV_clone_linemask = imdilate(PV_clone_linemask, disk_4);

                    % Check CV-PV line
                    j_mask_check=zeros(clone_image_height,clone_image_width,3);
                    j_mask_check(:,:,1)=bwperim(clone_mask);
                    j_mask_check(:,:,2)=bwperim(CV_Procr_linemask);
                    j_mask_check(:,:,3)=imadjust(rgb2gray(full_Procr_check));
                    file_location = fullfile(quant_outpath,sprintf(strcat(num2str(image_number),'_', num2str(clone_number),'_CV_PV_nuclei_Exp',num2str(Exp_sample_ID),'_X.jpg')));
                    imwrite(j_mask_check,file_location);

                    % Check CV-centroid line
                    k_mask_check=zeros(clone_image_height,clone_image_width,3);
                    k_mask_check(:,:,1)=bwperim(clone_mask);
                    k_mask_check(:,:,2)=bwperim(CV_centroid_linemask);
                    k_mask_check(:,:,3)=imadjust(rgb2gray(full_Procr_check));
                    file_location = fullfile(quant_outpath,sprintf(strcat(num2str(image_number),'_', num2str(clone_number),'_CV_centroid_nuclei_Exp',num2str(Exp_sample_ID),'_X.jpg')));
                    imwrite(k_mask_check,file_location);

                    % Check CV-clone line
                    l_mask_check=zeros(clone_image_height,clone_image_width,3);
                    l_mask_check(:,:,1)=bwperim(clone_mask);
                    l_mask_check(:,:,2)=bwperim(CV_clone_linemask);
                    l_mask_check(:,:,3)=imadjust(rgb2gray(full_Procr_check));
                    file_location = fullfile(quant_outpath,sprintf(strcat(num2str(image_number),'_', num2str(clone_number),'_CV_clone_nuclei_Exp',num2str(Exp_sample_ID),'_X.jpg')));
                    imwrite(l_mask_check,file_location);

                    % Check PV-clone line
                    m_mask_check=zeros(clone_image_height,clone_image_width,3);
                    m_mask_check(:,:,1)=bwperim(clone_mask);
                    m_mask_check(:,:,2)=bwperim(PV_clone_linemask);
                    m_mask_check(:,:,3)=imadjust(rgb2gray(full_Procr_check));
                    file_location = fullfile(quant_outpath,sprintf(strcat(num2str(image_number),'_', num2str(clone_number),'_PV_clone_nuclei_Exp',num2str(Exp_sample_ID),'_X.jpg')));
                    imwrite(m_mask_check,file_location);


                    %% Clone data calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    clone_sample_number = clone_sample_number + 1;
                    identifier = strcat('image_',num2str(image_number),'_','clone_',num2str(clone_number),'_X');
                    CV_centroid_distance = sqrt((clone_centroid_x - cv_x).^2 + (clone_centroid_y - cv_y).^2);
                    PV_centroid_distance = sqrt((clone_PV_edge_x - min_Procr_x).^2 + (clone_CV_edge_y - min_Procr_y).^2);
                    CV_centroid_percent = CV_centroid_distance/(CV_centroid_distance+ PV_centroid_distance)*100;
                    PV_centroid_percent = PV_centroid_distance/(CV_centroid_distance+ PV_centroid_distance)*100;
                    clone_area_pixels = clone_area_matrix(clone_number,:);
                    clone_out_row = {identifier, CV_centroid_percent, PV_centroid_percent, clone_area_pixels};
                    clone_data_cell_array(clone_sample_number,:) = clone_out_row;
                end
            end
        end 
    end
end



%% SECTION 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Export data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% csv file for clone data
clone_header = {'identifier', 'CV_centroid_percent_distance', 'PV_centroid_percent_distance', 'clone_area_pixels'};
clone_data_array = [clone_header; clone_data_cell_array];
writecell(clone_data_array,strcat(quant_outpath,'/Exp',num2str(Exp_sample_ID),'_raw_clone_data.csv'));

toc %ends timer and returns elapsed time
uiwait(msgbox('Your quantification has finished. All outputs have been saved to your output folder(s).','Instructions','modal'));



%% END OF CODE & ACKNOWLEDGEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ellen Rim and Berfin Azizoglu (Nusse Lab, Stanford University School of Medicine)
% Tanner Jensen and Colin Unger (Stanford University SOM and CS) 
% Uwe Schmidt (Max Planck Institute) and Martin Weigert (University of Chicago) for the StarDist FIJI plugin
    % downloaded October 2020 from FIJI update --> plugin menu
    % works cited: U Schmidt, M Weigert, C Broaddus, and G Myers. Cell Detection with Star-convex Polygons. International Conference on Medical Image Computing and Computer-Assisted Intervention (MICCAI), Granada, Spain, September 2018.
% Dylan Muir (SynSense/aiCTX) for the ReadImageJROI function
    % accessed October 2020 from https://www.mathworks.com/matlabcentral/fileexchange/33479-readimagejroi
    % works cited: DR Muir and BM Kampa. 2014. FocusStack and StimServer: A new open source MATLAB toolchain for visual stimulation and analysis of two-photon calcium neuronal imaging data, Frontiers in Neuroinformatics.