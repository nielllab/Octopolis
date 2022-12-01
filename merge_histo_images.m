%%% folder for each channel
folders = {'DAPI','FMRFA-620','FMRFR-570','PROTPRQFV-520','PRQFV-690'};
%%% base of name to store output files
name_base = 'merge\';
%%% amount to downsample to save memory
downsample = 0.125;

%%% get list of filenames from first channel folder
files = dir([folders{1} '\*.tif']);
for f = 1:length(files)
    fnames{f} = files(f).name;
end
fnames = sort(fnames)

%%% collect size of each individual image from first channel
for f = 1:length (fnames);
    f
    info = imfinfo(fullfile(folders{1},fnames{f}));
    rows(f) = info(1).Height;
    cols(f) = info(1).Width;
end

%%% set minimum crop size
nrows = min(rows);
ncols = min(cols);

%%% loop over all files
for f = 10:14  %length(files);
    fnames{f}
   tic
   for ch = 1:length(folders);
        %fullfile(folders{ch},files(f).name)
        im = imread(fullfile(folders{ch},fnames{f}));  %%% read file
        im_flat = max(im,[],3);                             %%% flatten the 3 RGB channels
        im_crop = imcrop(im_flat, [0 0 ncols nrows]);       %%% crop to miminum size
        im_small = imresize(im_crop,downsample);                 %%% downsize to make files manageable
        if ch==1
            imwrite(im_small,[name_base fnames{f}]);
        else
            imwrite(im_small,[name_base fnames{f}],'WriteMode','append');
        end
   end
    toc
end



