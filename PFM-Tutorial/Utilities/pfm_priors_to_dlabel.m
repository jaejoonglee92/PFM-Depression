function pfm_priors_to_dlabel(PriorsSpatial,Priors,Threshold,OutFile,OutDir,WorkbenchBinary,NetworkIdx)
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

end