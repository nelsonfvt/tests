%% Downsampling bvecs
% This function downsamples the angular resolution of a given bvec text
% file given a target angular resolution.
% The inputs are the following:
% raw_bvec = the full path to the raw bvec file;
% raw_bval = the full path to the raw bval file;
% output_txt = full path to the output file;
% output_index = full path to the file with the indexes of the new bvecs;
% target_angular_resolution = desired target angular resolution;
%%
function downsample_angular(raw_bvec, raw_bval, target_angular_resolution, output_bvec, output_bval)

    %check raw_bvec and raw_bval
    if ~exist(raw_bvec,'file')
       error(' Missing original bvec txt file') 
    end
    if ~exist(raw_bval,'file')
       error(' Missing original bval txt file') 
    end
 
    if isstring(target_angular_resolution)
        target_angular_resolution = str2double(target_angular_resolution);
    end

    % load raw_bvec and raw_bval
    raw = load(raw_bvec);
    if size(raw,1) == 3
        raw = raw'; % make sure raw is an Nx3 matrix
    end
    
    rawbval = load(raw_bval);
    if size(rawbval,1) == 1
        rawbval = rawbval'; % make sure raw is an Nx3 matrix
    end

    % check the number of b0s and the number of DWIs
    % It will only find the bzeros less than 20. This is hardcoded.
    num_b0 = size(find(rawbval <= 20)); 
    disp([' Your bval contains ',num2str(num_b0),' bzeros.'])
    num_bx = size(find(rawbval > 20));
    disp([' Your bval contains ',num2str(num_bx),' DWIs.'])
    
    % check logic issue
    if target_angular_resolution >= num_bx
       error(' The number of DWIs in the raw bvec file is smaller than your target, so no need to downsample.')  
    end

    %remove b0 from raw
    b0_index = find(rawbval <= 20);
    bx_index = find(rawbval > 20);
    b0 = raw(b0_index,:);
    bx = raw(bx_index,:);
    bval0 = rawbval(b0_index);
    bvalx = rawbval(bx_index);
    new_subset = core(bx,target_angular_resolution);

    newbvecs = [b0; new_subset(:,1:3)];  %vertical concatenation
    newbvals = [bval0; bvalx(new_subset(:,4))];
    
    % These are the original indices from the input DWI.
    index_original = bx_index(new_subset(:,4));
    index_original_zeroindex = index_original - 1;
    % Concatenating the original indices (1-index and 0-index)
    % and the new indices of the DWI bvecs without bzeros.
    index_all = horzcat(index_original, index_original_zeroindex, new_subset(:,4));
    
    % Defining the output name for the indices
    [outpath, n, e] = fileparts(output_bvec);
    output_index = strcat(outpath, '/new_bvecs_indices');
    
    dlmwrite(output_bvec, newbvecs', 'delimiter', '\t', 'precision',12)
    dlmwrite(output_bval, newbvals', 'delimiter', '\t', 'precision',12)
    dlmwrite(output_index, index_all', 'delimiter', '\t', 'precision',12)
    
end

function new_subset = core(original_subset,target_num)
    % created by Liang Zhan @ Aug 8, 2008
    % updated by Liang Zhan @ March 28, 2010

    % new_subset is num*4 format, 
    % the first three columns are coordinators (x,y,z) 
    % the last column is the index in the original vector

    % keep original_subset to be num*3 format
    x=size(original_subset);
    if x(2)>3
        original_subset=original_subset';
    end
    x=size(original_subset);
    num=x(1);
    original_q=zeros(num,4);
    for i=1:num
        original_q(i,1:3)=original_subset(i,:);
        original_q(i,4)=i;
    end

    index=1;
    new_subset(index,:)=original_q(index,:);
    original_q=setdiff(original_q,new_subset,'rows');

    while index<target_num
        dim1=size(new_subset);
        dim1=dim1(1);
        dim2=size(original_q);
        dim2=dim2(1);

        distance_sum=zeros(dim2,1);
        for i=1:dim2
            for j=1:dim1
                possible_vector=original_q(i,1:3); 
                vector=new_subset(j,1:3);
                distance_sum(i)=distance_sum(i)+cal_dis(possible_vector, vector);
            end        
        end
        object=find(distance_sum==max(distance_sum));
        if length(object)>1
            object=object(1);
        end
        index=index+1;    
        new_subset(index,:)=original_q(object,:);
        original_q=setdiff(original_q,new_subset,'rows');    
    end

    % new_subset(:,4)=[];
end


function dist = cal_dis(dir1, dir2)
    %adjust dir1 and dir2 to be 3*1 format
    dim=size(dir1);
    if dim(2)>1
        dir1=dir1';
    end

    dim=size(dir2);
    if dim(2)>1
        dir2=dir2';
    end

    [lon1,lat1,R] = cart2sph(dir1(1),dir1(2),dir1(3));
    [lon2,lat2,R] = cart2sph(dir2(1),dir2(2),dir2(3));
    [lon3,lat3,R] = cart2sph(-dir2(1),-dir2(2),-dir2(3));

    dist1 = greatCircleDistance(lat1,lon1,lat2,lon2,1);
    dist2 = greatCircleDistance(lat1,lon1,lat3,lon3,1);
    dist=1/(dist1^2+dist2^2);
end

function d = greatCircleDistance(phi_s, lambda_s, phi_f, lambda_f, r)
    % compute the great circle distance given lat and long for two points
    % optionally, a fifth parameter (r) can be specified. If this paramter
    % isn't specified it's assumed to be the mean radius of the earth. The
    % calculation is done using the Vincenty formula.
    %
    % INPUTS:
    % phi_s    = latitude of the standpoint (base) [rad]
    % lambda_s = longitude of the standpoint (base) [rad]
    % phi_f    = latitude of the forepoint (destination) [rad]
    % lambda_f = longitude of the forepoint (destination) [rad]
    % r        = radius of the sphere [units determine units of d]
    %
    % OUTPUT:
    % d        = great circle distance from standpoint to forepoint
    %
    % See http://en.wikipedia.org/wiki/Great-circle_distance

    % If no arguments, bail out
    if nargin < 4
        fprintf('Usage: greatCircleDistance(phi_s, lambda_s, phi_f, lambda_f, r)\n')
        return
    end

    % If no radius supplied, assume the mean radius of the earth in km
    if nargin < 5
        r = 6371.01; % km
    end

    % convert from degrees minutes seconds to radians as needed
    if isstruct(phi_s) || (length(phi_s) > 1 && ~isstruct(phi_s))
        phi_s = dms2r(phi_s);
    end
    if isstruct(lambda_s) || (length(lambda_s) > 1 && ~isstruct(lambda_s))
        lambda_s = dms2r(lambda_s);
    end
    if isstruct(phi_f) || (length(phi_f) > 1 && ~isstruct(phi_f))
        phi_f = dms2r(phi_f);
    end
    if isstruct(lambda_f) || (length(lambda_f) > 1 && ~isstruct(lambda_f))
        lambda_f = dms2r(lambda_f);
    end

    % Compute Delta lambda (delta longitude)
    Delta_lambda = lambda_f - lambda_s;

    % Compute Delta sigma (central angle)
    Delta_sigma = atan2(sqrt((cos(phi_f)*sin(Delta_lambda))^2 + (cos(phi_s)*sin(phi_f) - sin(phi_s)*cos(phi_f)*cos(Delta_lambda))^2), ...
        sin(phi_s)*sin(phi_f) + cos(phi_s)*cos(phi_f)*cos(Delta_lambda));

    d = r*Delta_sigma;
end

function r = dms2r(dms)
    if isstruct(dms)
        r = sign(dms.deg)*(abs(dms.deg) + (dms.min + dms.sec/60)/60)*pi/180;
    elseif length(dms) == 3
        r = sign(dms(1))*(abs(dms(1)) + (dms(2) + dms(3)/60)/60)*pi/180;
    elseif length(dms) == 2
        r = sign(dms(1))*(abs(dms(1)) + dms(2)/60)*pi/180;
    else
        r = nan;
    end
end
