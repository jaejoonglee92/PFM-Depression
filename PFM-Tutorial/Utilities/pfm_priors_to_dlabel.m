function pfm_priors_to_dlabel(PriorsSpatial,Priors,Threshold,OutFile,OutDir,WorkbenchBinary,NetworkIdx,MinSize,MidthickSurfs,CorrectedAreas,Fill)
% Threshold the spatial priors and write a single dlabel of the networks.
% Label keys follow Priors.NetworkLabels order (1 = Default_Parietal, ...,
% 20 = SomatoCognitiveAction) and colours come from Priors.NetworkColors, so the
% result is consistent with pfm_parse_manual_decisions_filter.
%
%   PriorsSpatial   path to priors_Spatial.dscalar.nii
%   Priors          the Priors struct (from priors.mat), or path to it
%   Threshold       inclusive: grayordinates with prior >= Threshold are kept;
%                   0 labels every grayordinate with any positive prior
%   OutFile         output basename (no extension)
%   OutDir          output directory
%   WorkbenchBinary path to wb_command
%   NetworkIdx      (optional) which prior columns to include, default all
%   MinSize         (optional) remove clusters of a network smaller than this
%                   (mm^2), per network, as pfm_spatial_filtering does for the
%                   individual PFM labels (MPCONE: 50); default 0 = no filtering
%   MidthickSurfs   {left, right} midthickness surfaces (needed if MinSize > 0);
%                   when given, network borders <OutFile>.L/R.border are also
%                   written (-cifti-label-to-border)
%   CorrectedAreas  (optional) {left, right} vertex-area metrics for the size
%                   threshold; for group-average surfaces use the HCP
%                   *_midthickness_va maps, so MinSize means individual cortex
%   Fill            (optional) true: removed grayordinates take the nearest
%                   surviving network (as in the PFM pipeline); false: they are
%                   left unlabelled. Default: true when Threshold == 0 (every
%                   cortical grayordinate labelled), false otherwise (removed
%                   islands are mostly isolated in unlabelled cortex).
%                   Grayordinates unlabelled before filtering stay unlabelled.
%
% The priors sum to <= 1 at every grayordinate, so from 0.5 up at most one
% network can pass and the labelling is unambiguous. Lower thresholds do let a
% grayordinate pass several (205 at 0.45, 2689 at 0.4); those go to the
% highest-scoring network, and the count is reported.

if ischar(Priors)
    tmp = load(Priors);
    Priors = tmp.Priors;
end

S = ft_read_cifti_mod(PriorsSpatial);
if nargin < 7 || isempty(NetworkIdx)
    NetworkIdx = 1:size(S.data,2);
end
NetworkIdx = NetworkIdx(:)';

D = S.data(:,NetworkIdx);
[Best,Which] = max(D,[],2);              % winner takes all if several pass
% a grayordinate with no prior at all (subcortex and cerebellum, which the
% priors do not cover) would otherwise pass a threshold of 0 and take column 1
% from the all-zero tie in max(), so a positive prior is always required
Keep = Best >= Threshold & Best > 0;

O = S;
O.data = zeros(size(D,1),1);
O.data(Keep) = NetworkIdx(Which(Keep));  % key = the network's column in the priors
O.mapname = {OutFile};
O.dimord = 'pos_scalar';

fprintf('threshold %.3g: %d of %d grayordinates labelled (%.1f%%)\n', ...
    Threshold, sum(Keep), size(D,1), 100*sum(Keep)/size(D,1));
nTie = sum(sum(D >= Threshold & D > 0,2) > 1);
if nTie > 0
    fprintf('  %d grayordinate(s) passed several networks; assigned to the highest prior\n', nTie);
end

% optional: remove small clusters of each network, as pfm_spatial_filtering
% does for the individual PFM labels (per network: wb_command
% -cifti-find-clusters with MinSize; removed grayordinates then take the
% nearest network, -cifti-dilate -nearest). Unlike pfm_spatial_filtering,
% only the removed grayordinates are filled (and only if Fill), so cortex the
% threshold left unlabelled stays unlabelled.
if nargin < 8 || isempty(MinSize); MinSize = 0; end
if nargin < 10; CorrectedAreas = []; end
if nargin < 11 || isempty(Fill); Fill = Threshold == 0; end
if MinSize > 0
    if nargin < 9 || numel(MidthickSurfs) ~= 2
        error('MinSize > 0 needs MidthickSurfs = {left, right} midthickness surfaces');
    end
    if isempty(CorrectedAreas)
        SurfArgs = [' -left-surface ' MidthickSurfs{1} ' -right-surface ' MidthickSurfs{2}];
    else
        SurfArgs = [' -left-surface ' MidthickSurfs{1} ' -corrected-areas ' CorrectedAreas{1} ...
            ' -right-surface ' MidthickSurfs{2} ' -corrected-areas ' CorrectedAreas{2}];
    end
    Pre = O.data;
    Kept = false(size(Pre));
    B = O; B.dimord = 'pos_time'; B.time = 0;
    for i = unique(nonzeros(Pre))'
        B.data = double(Pre == i);
        ft_write_cifti_mod([OutDir '/Tmp_clu'],B);
        Status = system([WorkbenchBinary ' -cifti-find-clusters ' OutDir '/Tmp_clu.dtseries.nii 0 ' num2str(MinSize) ...
            ' 0 ' num2str(MinSize) ' COLUMN ' OutDir '/Tmp_clu.dtseries.nii' SurfArgs ' -merged-volume']);
        if Status ~= 0; error('wb_command -cifti-find-clusters failed (status %d)',Status); end
        C = ft_read_cifti_mod([OutDir '/Tmp_clu.dtseries.nii']);
        Kept = Kept | C.data ~= 0;
    end
    Removed = Pre > 0 & ~Kept;
    Post = Pre; Post(Removed) = 0;
    if Fill && any(Removed)
        B.data = Post;
        ft_write_cifti_mod([OutDir '/Tmp_clu'],B);
        Status = system([WorkbenchBinary ' -cifti-dilate ' OutDir '/Tmp_clu.dtseries.nii COLUMN 50 50 -left-surface ' ...
            MidthickSurfs{1} ' -right-surface ' MidthickSurfs{2} ' ' OutDir '/Tmp_clu.dtseries.nii -nearest']);
        if Status ~= 0; error('wb_command -cifti-dilate failed (status %d)',Status); end
        C = ft_read_cifti_mod([OutDir '/Tmp_clu.dtseries.nii']);
        Post(Removed) = C.data(Removed);
    end
    delete([OutDir '/Tmp_clu.dtseries.nii']);
    fprintf('  clusters < %g mm^2 removed: %d grayordinates; %d refilled from the nearest network, %d unlabelled\n', ...
        MinSize, sum(Removed), sum(Removed & Post > 0), sum(Removed & Post == 0));
    for i = NetworkIdx
        if any(Pre == i & Removed) || any(Post == i & Removed)
            fprintf('    %-32s %5d -> %5d\n', Priors.NetworkLabels{i}, sum(Pre == i), sum(Post == i));
        end
    end
    O.data = Post;
end

% label list: name, then "key R G B A", for every network actually present
fid = fopen([OutDir '/LabelList.txt'],'w');
for i = NetworkIdx
    if any(O.data == i)
        c = round(Priors.NetworkColors(i,:) * 255);
        fprintf(fid,'%s\n%d %d %d %d 255\n',Priors.NetworkLabels{i},i,c(1),c(2),c(3));
        fprintf('  %2d  %-32s %4d %4d %4d  n = %d\n', ...
            i, Priors.NetworkLabels{i}, c(1), c(2), c(3), sum(O.data==i));
    end
end
fclose(fid);

ft_write_cifti_mod([OutDir '/Tmp_lab.dscalar.nii'],O);
Status = system([WorkbenchBinary ' -cifti-label-import ' OutDir '/Tmp_lab.dscalar.nii ' ...
    OutDir '/LabelList.txt ' OutDir '/' OutFile '.dlabel.nii -discard-others']);
delete([OutDir '/Tmp_lab.dscalar.nii'],[OutDir '/LabelList.txt']);
if Status ~= 0
    error('wb_command -cifti-label-import failed (status %d)',Status);
end

fprintf('wrote %s/%s.dlabel.nii\n',OutDir,OutFile);

% network borders (as pfm_parse_manual_decisions_filter does for the
% individual PFM labels), when surfaces are given; borders are tied to the
% 32k vertices, so they can be displayed on any 32k_fs_LR surface
if nargin >= 9 && numel(MidthickSurfs) == 2
    h = {'L','R'};
    for s = 1:2
        Status = system([WorkbenchBinary ' -cifti-label-to-border ' OutDir '/' OutFile '.dlabel.nii -border ' ...
            MidthickSurfs{s} ' ' OutDir '/' OutFile '.' h{s} '.border']);
        if Status ~= 0
            error('wb_command -cifti-label-to-border failed (status %d)',Status);
        end
    end
    fprintf('wrote %s/%s.{L,R}.border\n',OutDir,OutFile);
end

end